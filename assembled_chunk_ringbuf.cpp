#include <iostream>
#include "ch_frb_io_internals.hpp"
#include "chlog.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif

// This is a lightweight scoped lock
typedef std::lock_guard<std::mutex> guard_t;
// This is also a scoped lock that supports use of a condition variable.
typedef std::unique_lock<std::mutex> ulock_t;

assembled_chunk_ringbuf::assembled_chunk_ringbuf(const intensity_network_stream::initializer &ini_params_, int beam_id_, int stream_id_) :
    max_fpga_flushed(0),
    max_fpga_retrieved(0),
    first_fpgacount(0),
    first_packet_received(false),
    ini_params(ini_params_),
    stream_id(stream_id_),
    beam_id(beam_id_),
    output_devices(ini_params.output_devices)
{
    if ((beam_id < 0) || (beam_id > constants::max_allowed_beam_id))
	throw runtime_error("ch_frb_io: bad beam_id passed to assembled_chunk_ringbuf constructor");
    
    if (ini_params.assembled_ringbuf_capacity <= 0)
	throw runtime_error("ch_frb_io: assembled_chunk_ringbuf constructor: assembled_ringbuf_capacity must be > 0");

    if ((ini_params.nt_align < 0) || (ini_params.nt_align % constants::nt_per_assembled_chunk))
	throw runtime_error("ch_frb_io: 'nt_align' must be a multiple of nt_per_assembled_chunk(=" + to_string(constants::nt_per_assembled_chunk) + ")");

    for (int n: ini_params.telescoping_ringbuf_capacity) {
	if (n < 2)
	    throw runtime_error("ch_frb_io: assembled_chunk_ringbuf constructor: all telescoping_ringbuf_capacities must be >= 2");
    }

#ifndef __AVX2__
    if (ini_params.force_fast_kernels)
	throw runtime_error("ch_frb_io: the 'force_fast_kernels' flag was set, but this machine does not have the AVX2 instruction set");
#endif

    this->num_downsampling_levels = max(ini_params.telescoping_ringbuf_capacity.size(), 1UL);
    this->ringbuf_pos.resize(num_downsampling_levels, 0);
    this->ringbuf_size.resize(num_downsampling_levels, 0);
    this->ringbuf_capacity.resize(num_downsampling_levels, 0);
    this->ringbuf.resize(num_downsampling_levels);

    // Note that ringbuf_capacity[0] is the sum of 'ini_params.assembled_ringbuf_capacity'
    // and 'ini_params.telescoping_ringbuf_capacity[0]'.

    this->ringbuf_capacity[0] = ini_params.assembled_ringbuf_capacity;
    for (unsigned int i = 0; i < ini_params.telescoping_ringbuf_capacity.size(); i++)
	this->ringbuf_capacity[i] += ini_params.telescoping_ringbuf_capacity[i];

    for (int ids = 0; ids < num_downsampling_levels; ids++)
	this->ringbuf[ids].resize(ringbuf_capacity[ids]);

    this->downstream_pos = 0;
    this->downstream_bufsize = ini_params.assembled_ringbuf_capacity;
    
    this->_check_invariants();
}

void assembled_chunk_ringbuf::print_state() 
{
    guard_t lock(mutx);

    cout << "Beam " << beam_id << "\n";

    cout << "  downstream: [";
    for (int ipos = downstream_pos; ipos < ringbuf_pos[0] + ringbuf_size[0]; ipos++)
	cout << " " << this->ringbuf_entry(0,ipos)->ichunk;
    cout << " ]\n";
    
    for (int ids = 0; ids < num_downsampling_levels; ids++) {
	int i0 = ringbuf_pos[ids];
	int i1 = (ids > 0) ? (ringbuf_pos[ids] + ringbuf_size[ids]) : downstream_pos;
	
	cout << "  binning " << ids << ": [";
	for (int ipos = i0; ipos < i1; ipos++)
	    cout << " " << this->ringbuf_entry(ids,ipos)->ichunk;
	cout << " ]\n";
    }
}

shared_ptr<assembled_chunk>
assembled_chunk_ringbuf::find_assembled_chunk(uint64_t fpga_counts, bool top_level_only)
{
    ulock_t lock(mutx);

    // Return an empty pointer iff stream has ended, and chunk is requested past end-of-stream.
    // (If anything else goes wrong, an exception will be thrown.)
    if (this->doneflag && (fpga_counts >= this->final_fpga))
	return shared_ptr<assembled_chunk> ();
    
    // Scan telescoping ring buffer
    int start_level = (top_level_only ? 0 : num_downsampling_levels-1);
    for (int lev = start_level; lev >= 0; lev--) {
	for (int ipos = ringbuf_pos[lev]; ipos < ringbuf_pos[lev] + ringbuf_size[lev]; ipos++) {
	    auto ch = this->ringbuf_entry(lev, ipos);
	    if (ch->fpga_begin == fpga_counts)
		return ch;
	}
    }

    throw runtime_error("ch_frb_io::assembled_chunk::find_assembled_chunk(): couldn't find chunk, maybe your ring buffer is too small?");
}

vector<pair<shared_ptr<assembled_chunk>, uint64_t>>
assembled_chunk_ringbuf::get_ringbuf_snapshot(uint64_t min_fpga_counts, uint64_t max_fpga_counts)
{
    // Preallocate vector, before acquiring lock.
    vector<pair<shared_ptr<assembled_chunk>, uint64_t>> ret;
    ret.reserve(sum(ringbuf_capacity));

    guard_t lock(mutx);

    // Scan telescoping ring buffer, in a time-ordered way.
    for (int ids = num_downsampling_levels-1; ids >= 0; ids--) {
	for (int ipos = ringbuf_pos[ids]; ipos < ringbuf_pos[ids] + ringbuf_size[ids]; ipos++) {
	    auto chunk = this->ringbuf_entry(ids, ipos);

	    if (min_fpga_counts && (chunk->fpga_end <= min_fpga_counts))
		continue;   // no overlap
	    if (max_fpga_counts && (chunk->fpga_begin > max_fpga_counts))
		continue;   // no overlap

	    uint64_t where = 1 << (ids+1);   // Note: works since l1_ringbuf_level::L1RB_LEVELn == 2^n.
	    if ((ids == 0) && (ipos >= downstream_pos))
		where = l1_ringbuf_level::L1RB_DOWNSTREAM;

	    ret.push_back({ chunk, where });
	}
    }
    return ret;
}


// Returns stats about the ring buffer, for the get_statistics RPC.
void assembled_chunk_ringbuf::get_ringbuf_size(uint64_t *ringbuf_fpga_next,
                                               uint64_t *ringbuf_n_ready,
                                               uint64_t *ringbuf_total_capacity,
                                               uint64_t *ringbuf_nelements,
                                               uint64_t *ringbuf_fpga_min,
                                               uint64_t *ringbuf_fpga_max,
                                               int level) 
{
    guard_t lock(mutx);

    if (ringbuf_fpga_next && (level == 0)) {
	*ringbuf_fpga_next = 0;

	if (downstream_pos < ringbuf_pos[0] + ringbuf_size[0]) {
	    // Use initial FPGA count of first chunk which has been assembled,
	    // but not yet processed by "downstream" thread.
	    *ringbuf_fpga_next = this->ringbuf_entry(0, downstream_pos)->fpga_begin;
	}
	else if (ringbuf_size[0] > 0) {
	    // All chunks have been processed by "downstream" thread.
	    // Use final FPGA count of last chunk in buffer.
	    *ringbuf_fpga_next = this->ringbuf_entry(0, ringbuf_pos[0] + ringbuf_size[0] - 1)->fpga_end;
	}
    }

    if (ringbuf_n_ready && (level == 0)) {
	// Number of chunks which have been assembled, but not yet processed by "downstream" thread.
        *ringbuf_n_ready = ringbuf_pos[0] + ringbuf_size[0] - downstream_pos;
    }

    if (ringbuf_total_capacity) {
        if (level == 0) {
            *ringbuf_total_capacity = sum(ringbuf_capacity);
        } else if (level > num_downsampling_levels) {
            *ringbuf_total_capacity = 0;
        } else {
            *ringbuf_total_capacity = ringbuf_capacity[level-1];
        }
    }
    
    if (ringbuf_nelements) {
        if (level == 0) {
            *ringbuf_nelements = sum(ringbuf_size);
        } else if (level > num_downsampling_levels) {
            *ringbuf_nelements = 0;
        } else {
            *ringbuf_nelements = ringbuf_size[level-1];
        }
    }

    if (ringbuf_fpga_min) {
	*ringbuf_fpga_min = 0;
        if (level == 0) {
            for (int lev = num_downsampling_levels-1; lev >= 0; lev--) {
                if (ringbuf_size[lev] > 0) {
                    int ipos = ringbuf_pos[lev];
                    *ringbuf_fpga_min = this->ringbuf_entry(lev,ipos)->fpga_begin;
                    break;
                }
            }
        } else if (level <= num_downsampling_levels) {
            if (ringbuf_size[level-1] > 0) {
                int ipos = ringbuf_pos[level-1];
                *ringbuf_fpga_min = this->ringbuf_entry(level-1,ipos)->fpga_begin;
            }
        }
    }

    if (ringbuf_fpga_max) {
	*ringbuf_fpga_max = 0;
        if (level == 0) {
            for (int ids = 0; ids < num_downsampling_levels; ids++) {
                if (ringbuf_size[ids] > 0) {
                    int ipos = ringbuf_pos[ids] + ringbuf_size[ids] - 1;
                    *ringbuf_fpga_max = this->ringbuf_entry(ids,ipos)->fpga_end;
                    break;
                }
            }
        } else if (level <= num_downsampling_levels) {
            if (ringbuf_size[level-1] > 0) {
                int ipos = ringbuf_pos[level-1];
                *ringbuf_fpga_max = this->ringbuf_entry(level-1,ipos)->fpga_end;
            }
        }
    }
}


void assembled_chunk_ringbuf::stream_to_files(const string &filename_pattern, int priority, bool need_rfi)
{
    guard_t lock(mutx);
    this->stream_pattern = filename_pattern;
    this->stream_priority = priority;
    this->stream_rfi_mask = need_rfi;
    this->stream_chunks_written = 0;
    this->stream_bytes_written = 0;
}


// In assembled_chunk_ringbuf::put_unassembled_packet(), it's OK to modify 'event_counts' 
// without acquiring any locks.  This is because the assembler thread passes an event_subcounts 
// array which is updated on a per-packet basis, and accumulated into the global event_counts 
// on a per-udp_packet_list basis (with locks acquired!).

void assembled_chunk_ringbuf::put_unassembled_packet(const intensity_packet &packet, int64_t *event_counts)
{
    uint64_t packet_t0 = packet.fpga_count / packet.fpga_counts_per_sample;
    uint64_t packet_ichunk = packet_t0 / constants::nt_per_assembled_chunk;

    if (!first_packet_received) {
	uint64_t first_ichunk = packet_ichunk;

	if (ini_params.nt_align > 0) {
	    uint64_t chunk_align = ini_params.nt_align / constants::nt_per_assembled_chunk;
	    first_ichunk = ((first_ichunk + chunk_align - 1) / chunk_align) * chunk_align;
	}
	
	this->active_chunk0 = this->_make_assembled_chunk(first_ichunk, 1);
	this->active_chunk1 = this->_make_assembled_chunk(first_ichunk+1, 1);
	this->first_packet_received = true;

	// We initialize 'first_fpgacount' to the FPGA count of the first assembled_chunk.
	// (Note that this can be either earlier or later than the FPGA count of the packet.)
	// This makes sense because 'first_fpgacount' is used to convert between FPGA counts and
	// time sample indices in rf_pipelines/bonsai.
	
        this->first_fpgacount = first_ichunk * constants::nt_per_assembled_chunk * ini_params.fpga_counts_per_sample;
    }

    // We test these pointers instead of 'doneflag' so that we don't need to acquire the lock in every call.
    if (_unlikely(!active_chunk0 || !active_chunk1))
	throw runtime_error("ch_frb_io: internal error: assembled_chunk_ringbuf::put_unassembled_packet() called after end_stream()");

    if (packet_ichunk >= active_chunk0->ichunk + 2) {
	//
	// If we receive a packet whose timestamps extend past the range of our current
	// assembly buffer, then we advance the buffer and send an assembled_chunk to the
	// "downstream" thread.
	//
	// A design decision here: for a packet which is far in the future, we advance the 
	// buffer by one assembled_chunk, rather than advancing all the way to the packet
	// timestamp.  This is to avoid a situation where a single rogue packet timestamped
	// in the far future effectively kills the L1 node.
	//

        /*
         chlog("Got packet with ichunk = " << packet_ichunk << ", vs active " << active_chunk1->ichunk << ", " << active_chunk0->ichunk
         << " (" << active_chunk1->packets_received << ", " << active_chunk0->packets_received
         << " packets received, " << active_chunk0->packets_missed << " missed) -- sender " << ip_to_string(packet.sender));
         */

	this->_put_assembled_chunk(active_chunk0, event_counts);

        // After _put_assembled_chunk(), active_chunk0 has been reset to a null pointer.
        active_chunk0.swap(active_chunk1);

        // Note that we've just swapped active_chunk1 down to active_chunk0, so active_chunk1's ichunk is (active0 + 1).
	active_chunk1 = this->_make_assembled_chunk(active_chunk0->ichunk + 1, 1);
    }

    if (packet_ichunk == active_chunk0->ichunk) {
	event_counts[intensity_network_stream::event_type::assembler_hit]++;
	active_chunk0->add_packet(packet);
    }
    else if (packet_ichunk == active_chunk1->ichunk) {
	event_counts[intensity_network_stream::event_type::assembler_hit]++;
	active_chunk1->add_packet(packet);
    }
    else {
      //chlog("Assembler miss (packet_ichunk " << packet_ichunk << " vs active " << active_chunk0->ichunk << " and " << active_chunk1->ichunk
      //<< "), from " << ip_to_string(packet.sender));
	event_counts[intensity_network_stream::event_type::assembler_miss]++;
	active_chunk0->packets_missed++;
	if (_unlikely(ini_params.throw_exception_on_assembler_miss))
	    throw runtime_error("ch_frb_io: assembler miss occurred, and this stream was constructed with the 'throw_exception_on_assembler_miss' flag");
    }
}

struct streaming_write_chunk_request : public write_chunk_request {
    weak_ptr<assembled_chunk_ringbuf> assembler;
    int udelay;
    virtual void status_changed(bool finished, bool success,
                                const std::string &state,
                                const std::string &error_message) override {
        if (udelay) {
            usleep(udelay);
        }
        if (finished && success) {
            // "lock" our weak pointer to the assembler; this fails if
            // it has been deleted already (in which case we do nothing).
            shared_ptr<assembled_chunk_ringbuf> realpointer = assembler.lock();
            if (realpointer)
                realpointer->chunk_streamed(filename);
	    else
                cout << "Assembled_chunk_ringbuffer: write chunk finished, but assembler has been deleted. No problem!" << endl;
        }
    }
    virtual ~streaming_write_chunk_request() { }
};

void assembled_chunk_ringbuf::chunk_streamed(const std::string &filename) {
    //chlog("Assembled_chunk streamed: " << filename);
    struct stat st;
    int err = stat(filename.c_str(), &st);

    if (err < 0) {
	chlog("warning: failed to stat file " + filename + " that was just streamed: " + strerror(errno));
	return;
    }

    size_t len = st.st_size;
    guard_t lock(mutx);
    this->stream_chunks_written ++;
    this->stream_bytes_written += len;
}

void assembled_chunk_ringbuf::get_streamed_chunks(int &achunks, size_t &abytes) {
    guard_t lock(mutx);
    achunks = stream_chunks_written;
    abytes = stream_bytes_written;
}

// Helper function called assembler thread, to add a new assembled_chunk to the ring buffer.
// Resets 'chunk' to a null pointer.
// Warning: only safe to call from assembler thread.
bool assembled_chunk_ringbuf::_put_assembled_chunk(unique_ptr<assembled_chunk> &chunk, int64_t *event_counts)
{
    if (!chunk)
	throw runtime_error("ch_frb_io: internal error: empty pointer passed to assembled_chunk_ringbuf::_put_unassembled_packet()");
    if (chunk->has_rfi_mask)
	throw runtime_error("ch_frb_io: internal error: chunk passed to assembled_chunk_ringbuf::_put_unassembled_packet() has rfi_mask flag set");

    //chlog("Assembled chunk " << chunk->ichunk << " beam " << beam_id << ": received " << chunk->packets_received << " packets");

    // Step 1: prepare all data needed to modify the ring buffer.  In this step, we do all of our
    // buffer allocation and downsampling, without the lock held.  In step 2, we will acquire the
    // lock and modify the ring buffer (without expensive operations like allocation/downsampling).
    //
    // It is very important to note that we can read (but not modify) the ring buffer without 
    // acquiring the lock!  This is because _put_assembled_chunk() is called from the assembler
    // thread, and only the assembler thread modifies the ring buffer ("single-producer").  

    int nds = this->num_downsampling_levels;

    // List of chunks to be pushed and popped at each level of the ring buffer (in step 2!)
    vector<shared_ptr<assembled_chunk>> pushlist(nds);
    vector<shared_ptr<assembled_chunk>> poplist(2*nds);

    uint64_t chunk_fpga_end = chunk->fpga_end;
    
    // Converts unique_ptr -> shared_ptr, and resets 'chunk' to a null pointer.
    pushlist[0] = shared_ptr<assembled_chunk> (chunk.release());

    // Without lock held...
    for (int ids = 0; ids < nds; ids++) {
	// At top of loop, we want to add the chunk pushlist[ids] at level 'ids' of
	// the telescoping ring buffer.  Is there space available...?

	if (ringbuf_size[ids] < ringbuf_capacity[ids])
	    break;  // ... Yes, no problem.

	// ... No space available!  Need to pop chunks.
	// If we're at the bottom level of the buffer, just pop a single chunk...

	if (ids == nds-1) {
	    poplist[2*ids] = this->ringbuf_entry(ids, ringbuf_pos[ids]);

	    // This assert and its counterpart below ensure that a chunk never leaves the telescoping
	    // ring buffer before its RFI mask is filled.  (If this could happen, we might hang on to
	    // the reference forever in output_device::_awaiting_rfi and get a memory leak.)

	    if ((ini_params.nrfifreq > 0) && !poplist[2*ids]->has_rfi_mask)
		throw runtime_error("ch_frb_io: _put_assembled_chunk(): rfimask not initialized as expected, maybe your ring buffer is too small?");

	    break;
	}

	// ... Otherwise, pop two chunks, downsample, and push the downsampled chunk
	// to the next level of the telescoping ring buffer.

	poplist[2*ids] = this->ringbuf_entry(ids, ringbuf_pos[ids]);
	poplist[2*ids+1] = this->ringbuf_entry(ids, ringbuf_pos[ids]+1);
	
	if ((ini_params.nrfifreq > 0) && (!poplist[2*ids]->has_rfi_mask || !poplist[2*ids+1]->has_rfi_mask))
	    throw runtime_error("ch_frb_io: _put_assembled_chunk(): rfimask not initialized as expected, maybe your ring buffer is too small?");

	pushlist[ids+1] = _make_assembled_chunk(poplist[2*ids]->ichunk, 1 << (ids+1));

	// Note: this test is currently superfluous, since _make_assembled_chunk() throws
	// an exception (rather than returning NULL) if the allocation fails.  It's 
	// just a placeholder to remind myself that the return value of this function
	// is supposed to indicate success/failure, and that more thought needs to
	// be put into assembled_chunk memory management.

	if (!pushlist[ids+1])
	    return false;

	pushlist[ids+1]->downsample(poplist[2*ids].get(), poplist[2*ids+1].get());
    }

    // Step 2: acquire lock and modify the ring buffer.  We have already computed the chunks to
    // be added/removed at each level (pushlist/poplist), so we don't malloc/free/downsample with
    // the lock held.

    ulock_t lock(mutx);
    
    if (this->doneflag)
	throw runtime_error("ch_frb_io: internal error: assembled_chunk_ringbuf::put_unassembled_packet() called after end_stream()");

    for (int ids = 0; ids < nds; ids++) {
	// Number of chunks to be removed from level 'ids' of the telescoping ring buffer.
	int npop = 0;
	if (poplist[2*ids]) npop++;
	if (poplist[2*ids+1]) npop++;
	
	// Remove chunks from ring buffer, by resetting shared_ptrs
	// Note that we are still holding references to these chunks in poplist[].
	// This ensures that assembled_chunk destructors are called without the lock held (at the end of this function).

	for (int p = 0; p < npop; p++)
	    this->ringbuf_entry(ids, ringbuf_pos[ids]+p) = shared_ptr<assembled_chunk> ();

	ringbuf_pos[ids] += npop;
	ringbuf_size[ids] -= npop;

	// Add chunk to level 'ids' of the telescoping ring buffer.

	if (pushlist[ids]) {
	    this->ringbuf_entry(ids, ringbuf_pos[ids] + ringbuf_size[ids]) = pushlist[ids];
	    ringbuf_size[ids]++;
	}
    }

    // Last step while holding lock: handle case where downstream thread is running
    // slow, and chunks were dropped.

    int num_assembled_chunks_dropped = 0;
    int max_allowed_downstream_pos = ringbuf_pos[0] + ringbuf_size[0] - downstream_bufsize;

    if (downstream_pos < max_allowed_downstream_pos) {
	num_assembled_chunks_dropped = max_allowed_downstream_pos - downstream_pos;
	downstream_pos = max_allowed_downstream_pos;
    }

    // Make thread-local copies with lock held.
    string loc_stream_pattern = this->stream_pattern;
    int loc_stream_priority = this->stream_priority;
    bool loc_stream_rfi_mask = this->stream_rfi_mask;
    
    this->cond_assembled_chunks_added.notify_all();
    lock.unlock();

    // Stream new chunk to disk (if 'stream_pattern' is a nonempty string).
    // It's better to do this processing without the lock held, we just need to use
    // 'loc_stream_pattern' and 'loc_stream_priority' here, for thread-safety.

    if (loc_stream_pattern.size() > 0) {
	shared_ptr<streaming_write_chunk_request> wreq = make_shared<streaming_write_chunk_request> ();
	wreq->filename = pushlist[0]->format_filename(loc_stream_pattern);
	wreq->priority = loc_stream_priority;
        wreq->need_wait = loc_stream_rfi_mask;
	// DEBUG
	if (wreq->priority == -1000)
            wreq->udelay = 1000000;
	wreq->chunk = pushlist[0];	
        wreq->assembler = shared_from_this();
        
	// return value from enqueue_write_request() is ignored.
	output_devices.enqueue_write_request(wreq);
    }

    // This call to _check_invariants() is a good test during debugging, but
    // shouldn't be enabled in production.
    //
    // FIXME!!  Make sure this line gets commented out eventually.
    //this->_check_invariants();

    // For even more debugging, uncomment this line!
    // this->print_state();

    if (event_counts) {
	event_counts[intensity_network_stream::event_type::assembled_chunk_queued]++;
	event_counts[intensity_network_stream::event_type::assembled_chunk_dropped] += num_assembled_chunks_dropped;
    }

    assert(chunk_fpga_end > this->max_fpga_flushed);
    this->max_fpga_flushed = chunk_fpga_end;

    if (ini_params.emit_warning_on_buffer_drop && (num_assembled_chunks_dropped > 0))
	cout << "ch_frb_io: warning: processing thread is running too slow, dropping assembled_chunk" << endl;
    if (ini_params.throw_exception_on_buffer_drop && (num_assembled_chunks_dropped > 0))
	throw runtime_error("ch_frb_io: assembled_chunk was dropped and stream was constructed with 'throw_exception_on_buffer_drop' flag");

    // Note: when this function returns, stray references in poplist[*][*] are dropped, and assembled_chunk destructors get called.
    return true;
}


void assembled_chunk_ringbuf::_check_invariants()
{
    // It's OK to access the ringbuf_* fields read-only without acquiring the lock,
    // since _check_invariants() is only called from the assembler thread.
    //
    // Some checks in this function are redundant with checks elsewhere, but that's OK!

    ch_assert(num_downsampling_levels > 0);
    ch_assert(ringbuf_pos.size() == (unsigned) num_downsampling_levels);
    ch_assert(ringbuf_size.size() == (unsigned) num_downsampling_levels);
    ch_assert(ringbuf_capacity.size() == (unsigned) num_downsampling_levels);
    ch_assert(ringbuf.size() == (unsigned) num_downsampling_levels);

    for (int ids = 0; ids < num_downsampling_levels; ids++) {
	ch_assert(ringbuf_pos[ids] >= 0);
	ch_assert(ringbuf_size[ids] >= 0);
	ch_assert(ringbuf_capacity[ids] >= 2);
	ch_assert(ringbuf_size[ids] <= ringbuf_capacity[ids]);
	ch_assert(ringbuf[ids].size() == (unsigned) ringbuf_capacity[ids]);

	for (int ipos = ringbuf_pos[ids]; ipos < ringbuf_pos[ids] + ringbuf_capacity[ids]; ipos++) {
	    shared_ptr<assembled_chunk> chunk = this->ringbuf_entry(ids, ipos);

	    // These entries of the ring buffer should be empty.
	    if (ipos >= ringbuf_pos[ids] + ringbuf_size[ids]) {
		ch_assert(!chunk);
		continue;
	    }
	    
	    // Nonempty entries...
	    ch_assert(chunk);
	    ch_assert(chunk->beam_id == this->beam_id);
	    ch_assert(chunk->nupfreq == this->ini_params.nupfreq);
	    ch_assert(chunk->nt_per_packet == this->ini_params.nt_per_packet);
	    ch_assert(chunk->fpga_counts_per_sample == this->ini_params.fpga_counts_per_sample);
	    ch_assert(chunk->binning == (1 << ids));
	    ch_assert(chunk->isample == chunk->ichunk * constants::nt_per_assembled_chunk);

	    if ((ini_params.nrfifreq > 0) && (ids > 0))
		ch_assert(chunk->has_rfi_mask);

	    // Now check logical contiguousness of the telescoping ring buffer, by
	    // checking that 'chunk' is contiguous with the next chunk in the buffer.

	    shared_ptr<assembled_chunk> next;

	    if (ipos < ringbuf_pos[ids] + ringbuf_size[ids] - 1) {
		// Next chunk is in same level of telescoping ring buffer.
		next = this->ringbuf_entry(ids, ipos+1);
	    }
	    else if (ids > 0) {
		// Next chunk is in a different level of the telescoping ring buffer.
		ch_assert(ringbuf_size[ids-1] > 0);
		next = this->ringbuf_entry(ids-1, ringbuf_pos[ids-1]);
	    }
	    else
		continue;   // Last chunk in buffer, there is no 'next'

	    ch_assert(next);
	    ch_assert(next->ichunk == chunk->ichunk + chunk->binning);
	}
    }

    // We do need to acquire the lock to access 'downstream_pos', since it's modified
    // by the downstream thread.

    int dpos;
    {
        guard_t lock(mutx);
        dpos = this->downstream_pos;
    }

    ch_assert(downstream_bufsize > 0);
    ch_assert(downstream_bufsize <= ringbuf_capacity[0]);

    ch_assert(dpos >= ringbuf_pos[0]);
    ch_assert(dpos <= ringbuf_pos[0] + ringbuf_size[0]);
    ch_assert(dpos >= ringbuf_pos[0] + ringbuf_size[0] - downstream_bufsize);
}


bool assembled_chunk_ringbuf::inject_assembled_chunk(assembled_chunk* chunk) 
{
    uint64_t ich = chunk->ichunk;
    unique_ptr<assembled_chunk> uch(chunk);
    bool worked = _put_assembled_chunk(uch, NULL);
    // Danger: monkey with the active_chunk0, active_chunk1 variables,
    // which are not lock-protected and only supposed to be accessed
    // by the assembler thread.
    active_chunk0 = this->_make_assembled_chunk(ich + 1, 1);
    active_chunk1 = this->_make_assembled_chunk(ich + 2, 1);
    return worked;
}


shared_ptr<assembled_chunk> assembled_chunk_ringbuf::get_assembled_chunk(bool wait)
{
    shared_ptr<assembled_chunk> chunk;
    ulock_t lock(mutx);

    for (;;) {
	if (downstream_pos < ringbuf_pos[0] + ringbuf_size[0]) {
	    chunk = this->ringbuf_entry(0, downstream_pos);
	    downstream_pos++;
	    break;
	}

        if (!wait)
            break;

	if (this->doneflag)
	    break;   // Ring buffer is empty and end_stream() has been called

	// Wait for chunks to be added to the ring buffer.
        this->cond_assembled_chunks_added.wait(lock);
    }

    if (chunk) {
        assert(chunk->fpga_end > this->max_fpga_retrieved);
        this->max_fpga_retrieved = chunk->fpga_end;
    }
    return chunk;
}


// Called by the assembler thread, when it exits.
void assembled_chunk_ringbuf::end_stream(int64_t *event_counts)
{
    if (!active_chunk0 || !active_chunk1)
	throw runtime_error("ch_frb_io: internal error: empty pointers in assembled_chunk_ringbuf::end_stream(), this can happen if end_stream() is called twice");

    // Local variable (will shortly assign to this->final_fpga, after acquiring lock).
    uint64_t loc_final_fpga = (active_chunk0->ichunk + 2) * uint64_t(constants::nt_per_assembled_chunk * active_chunk0->fpga_counts_per_sample);

    // After these calls, 'active_chunk0' and 'active_chunk1' will be reset to null pointers.
    this->_put_assembled_chunk(active_chunk0, event_counts);
    this->_put_assembled_chunk(active_chunk1, event_counts);

    {
        ulock_t lock(mutx);
        if (doneflag)
            throw runtime_error("ch_frb_io: internal error: doneflag already set in assembled_chunk_ringbuf::end_stream()");

        // With lock held
        this->doneflag = true;
        this->final_fpga = loc_final_fpga;
        // Wake up processing thread, if it is waiting for data
        this->cond_assembled_chunks_added.notify_all();
    }
}


std::unique_ptr<assembled_chunk> assembled_chunk_ringbuf::_make_assembled_chunk(uint64_t ichunk, int binning, bool zero)
{
    struct assembled_chunk::initializer chunk_params;

    chunk_params.beam_id = this->beam_id;
    chunk_params.nupfreq = this->ini_params.nupfreq;
    chunk_params.nrfifreq = this->ini_params.nrfifreq;
    chunk_params.nt_per_packet = this->ini_params.nt_per_packet;
    chunk_params.fpga_counts_per_sample = this->ini_params.fpga_counts_per_sample;
    chunk_params.force_reference = this->ini_params.force_reference_kernels;
    chunk_params.force_fast = this->ini_params.force_fast_kernels;
    chunk_params.stream_id = this->stream_id;
    chunk_params.binning = binning;
    chunk_params.ichunk = ichunk;

    if (ini_params.memory_pool) {
	chunk_params.pool = ini_params.memory_pool;
	chunk_params.slab = ini_params.memory_pool->get_slab(zero);

	if (!chunk_params.slab)
	    throw runtime_error("**** Too much memory pressure for this poor L1 node to survive!  Blowing up now... ****");
    }

    return assembled_chunk::make(chunk_params);
}


}  // namespace ch_frb_io
