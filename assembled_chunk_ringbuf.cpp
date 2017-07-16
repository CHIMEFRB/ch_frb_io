#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


assembled_chunk_ringbuf::assembled_chunk_ringbuf(const intensity_network_stream::initializer &ini_params_, int beam_id_, int nupfreq_,
						 int nt_per_packet_, uint64_t fpga_counts_per_sample_, uint64_t fpga_count0) :
    ini_params(ini_params_),
    beam_id(beam_id_),
    nupfreq(nupfreq_),
    nt_per_packet(nt_per_packet_),
    fpga_counts_per_sample(fpga_counts_per_sample_)
{
    if ((beam_id < 0) || (beam_id > constants::max_allowed_beam_id))
	throw runtime_error("ch_frb_io: bad beam_id passed to assembled_chunk_ringbuf constructor");
    if ((nupfreq < 0) || (nupfreq > constants::max_allowed_nupfreq))
	throw runtime_error("ch_frb_io: bad nupfreq value passed to assembled_chunk_ringbuf constructor");
    if ((nt_per_packet <= 0) || (constants::nt_per_assembled_chunk % nt_per_packet != 0))
	throw runtime_error("ch_frb_io: internal error: assembled_chunk_ringbuf::nt_per_packet must be a divisor of constants::nt_per_assembled_chunk");
    if ((fpga_counts_per_sample <= 0) || (fpga_counts_per_sample > constants::max_allowed_fpga_counts_per_sample))
	throw runtime_error("ch_frb_io: bad fpga_counts_per_sample value passed to assembled_chunk_ringbuf constructor");
    if (fpga_count0 % fpga_counts_per_sample != 0)
	throw runtime_error("ch_frb_io: assembled_chunk_ringbuf constructor: fpga_count0 was not a multiple of fpga_counts_per_sample");
    if (ini_params.assembled_ringbuf_capacity <= 0)
	throw runtime_error("ch_frb_io: assembled_chunk_ringbuf constructor: assembled_ringbuf_capacity must be > 0");

    for (int n: ini_params.telescoping_ringbuf_capacity) {
	if (n < 2)
	    throw runtime_error("ch_frb_io: assembled_chunk_ringbuf constructor: all telescoping_ringbuf_capacities must be >= 2");
    }

#ifndef __AVX2__
    if (ini_params.force_fast_kernels)
	throw runtime_error("ch_frb_io: the 'force_fast_kernels' flag was set, but this machine does not have the AVX2 instruction set");
#endif

    uint64_t packet_t0 = fpga_count0 / fpga_counts_per_sample;
    uint64_t ichunk = packet_t0 / constants::nt_per_assembled_chunk;

    this->active_chunk0 = this->_make_assembled_chunk(ichunk, 1);
    this->active_chunk1 = this->_make_assembled_chunk(ichunk+1, 1);

    pthread_mutex_init(&this->lock, NULL);
    pthread_cond_init(&this->cond_assembled_chunks_added, NULL);

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


assembled_chunk_ringbuf::~assembled_chunk_ringbuf()
{
    pthread_cond_destroy(&this->cond_assembled_chunks_added);
    pthread_mutex_destroy(&this->lock);
}


void assembled_chunk_ringbuf::print_state() 
{
    pthread_mutex_lock(&this->lock);

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

    pthread_mutex_unlock(&this->lock);
}


vector<pair<shared_ptr<assembled_chunk>, uint64_t>>
assembled_chunk_ringbuf::get_ringbuf_snapshot(uint64_t min_fpga_counts, uint64_t max_fpga_counts)
{
    // Preallocate vector, before acquiring lock.
    vector<pair<shared_ptr<assembled_chunk>, uint64_t>> ret;
    ret.reserve(sum(ringbuf_capacity));

    pthread_mutex_lock(&this->lock);

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

    pthread_mutex_unlock(&this->lock);
    return ret;
}


// Returns stats about the ring buffer, for the get_statistics RPC.
void assembled_chunk_ringbuf::get_ringbuf_size(uint64_t *ringbuf_fpga_next,
                                               uint64_t *ringbuf_n_ready,
                                               uint64_t *ringbuf_total_capacity,
                                               uint64_t *ringbuf_nelements,
                                               uint64_t *ringbuf_fpga_min,
                                               uint64_t *ringbuf_fpga_max) 
{
    pthread_mutex_lock(&this->lock);

    if (ringbuf_fpga_next) {
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

    if (ringbuf_n_ready) {
	// Number of chunks which have been assembled, but not yet processed by "downstream" thread.
        *ringbuf_n_ready = ringbuf_pos[0] + ringbuf_size[0] - downstream_pos;
    }

    if (ringbuf_total_capacity)
        *ringbuf_total_capacity = sum(ringbuf_capacity);

    if (ringbuf_nelements)
        *ringbuf_nelements = sum(ringbuf_size);

    if (ringbuf_fpga_min) {
	*ringbuf_fpga_min = 0;
	for (int ids = num_downsampling_levels-1; ids >= 0; ids--) {
	    if (ringbuf_size[ids] > 0) {
		int ipos = ringbuf_pos[ids];
		*ringbuf_fpga_min = this->ringbuf_entry(ids,ipos)->fpga_begin;
		break;
	    }
	}
    }

    if (ringbuf_fpga_max) {
	*ringbuf_fpga_max = 0;
	for (int ids = 0; ids < num_downsampling_levels; ids++) {
	    if (ringbuf_size[ids] > 0) {
		int ipos = ringbuf_pos[ids] + ringbuf_size[ids] - 1;
		*ringbuf_fpga_max = this->ringbuf_entry(ids,ipos)->fpga_end;
		break;
	    }
	}
    }

    pthread_mutex_unlock(&this->lock);
}


// In assembled_chunk_ringbuf::put_unassembled_packet(), it's OK to modify 'event_counts' 
// without acquiring any locks.  This is because the assembler thread passes an event_subcounts 
// array which is updated on a per-packet basis, and accumulated into the global event_counts 
// on a per-udp_packet_list basis (with locks acquired!).

void assembled_chunk_ringbuf::put_unassembled_packet(const intensity_packet &packet, int64_t *event_counts)
{
    // We test these pointers instead of 'doneflag' so that we don't need to acquire the lock in every call.
    if (_unlikely(!active_chunk0 || !active_chunk1))
	throw runtime_error("ch_frb_io: internal error: assembled_chunk_ringbuf::put_unassembled_packet() called after end_stream()");

    uint64_t packet_t0 = packet.fpga_count / packet.fpga_counts_per_sample;
    uint64_t packet_ichunk = packet_t0 / constants::nt_per_assembled_chunk;

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
	event_counts[intensity_network_stream::event_type::assembler_miss]++;
	if (_unlikely(ini_params.throw_exception_on_assembler_miss))
	    throw runtime_error("ch_frb_io: assembler miss occurred, and this stream was constructed with the 'throw_exception_on_assembler_miss' flag");
    }
}


// Helper function called assembler thread, to add a new assembled_chunk to the ring buffer.
// Resets 'chunk' to a null pointer.
// Warning: only safe to call from assembler thread.
bool assembled_chunk_ringbuf::_put_assembled_chunk(unique_ptr<assembled_chunk> &chunk, int64_t *event_counts)
{
    if (!chunk)
	throw runtime_error("ch_frb_io: internal error: empty pointer passed to assembled_chunk_ringbuf::_put_unassembled_packet()");

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
	    break;
	}

	// ... Otherwise, pop two chunks, downsample, and push the downsampled chunk
	// to the next level of the telescoping ring buffer.

	poplist[2*ids] = this->ringbuf_entry(ids, ringbuf_pos[ids]);
	poplist[2*ids+1] = this->ringbuf_entry(ids, ringbuf_pos[ids]+1);	

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

    // We process the stream_filename and chunk_callbacks here (without the lock held!)

    if (stream_filename_pattern.length()) {
        string fn = pushlist[0]->format_filename(stream_filename_pattern);
        // turn on compression, but revert the state of pushlist[0]->msgpack_bitshuffle after writing.
        bool bitpack = pushlist[0]->msgpack_bitshuffle;
        pushlist[0]->msgpack_bitshuffle = true;
        pushlist[0]->write_msgpack_file(fn);
        pushlist[0]->msgpack_bitshuffle = bitpack;
    }

    for (auto it=chunk_callbacks.begin(); it!=chunk_callbacks.end(); it++) {
	(*it)(pushlist[0]);
    }

    // Step 2: acquire lock and modify the ring buffer.  We have already computed the chunks to
    // be added/removed at each level (pushlist/poplist), so we don't malloc/free/downsample with
    // the lock held.

    pthread_mutex_lock(&this->lock);

    if (this->doneflag) {
	pthread_mutex_unlock(&this->lock);
	throw runtime_error("ch_frb_io: internal error: assembled_chunk_ringbuf::put_unassembled_packet() called after end_stream()");
    }

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

    pthread_cond_broadcast(&this->cond_assembled_chunks_added);
    pthread_mutex_unlock(&this->lock);

    // This call to _check_invariants() is a good test during debugging, but
    // shouldn't be enabled in production.
    //
    // FIXME!!  Make sure this line gets commented out eventually.
    this->_check_invariants();

    // For even more debugging, uncomment this line!
    // this->print_state();

    if (event_counts) {
	event_counts[intensity_network_stream::event_type::assembled_chunk_queued]++;
	event_counts[intensity_network_stream::event_type::assembled_chunk_dropped] += num_assembled_chunks_dropped;
    }

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
	    ch_assert(chunk->nupfreq == this->nupfreq);
	    ch_assert(chunk->nt_per_packet == this->nt_per_packet);
	    ch_assert(chunk->fpga_counts_per_sample == this->fpga_counts_per_sample);
	    ch_assert(chunk->binning == (1 << ids));
	    ch_assert(chunk->isample == chunk->ichunk * constants::nt_per_assembled_chunk);

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

    pthread_mutex_lock(&lock);
    int dpos = this->downstream_pos;
    pthread_mutex_unlock(&lock);

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
    pthread_mutex_lock(&this->lock);

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
	pthread_cond_wait(&this->cond_assembled_chunks_added, &this->lock);
    }

    pthread_mutex_unlock(&this->lock);
    return chunk;
}


void assembled_chunk_ringbuf::end_stream(int64_t *event_counts)
{
    if (!active_chunk0 || !active_chunk1)
	throw runtime_error("ch_frb_io: internal error: empty pointers in assembled_chunk_ringbuf::end_stream(), this can happen if end_stream() is called twice");

    if (active_chunk0)
        this->_put_assembled_chunk(active_chunk0, event_counts);

    if (active_chunk1)
        this->_put_assembled_chunk(active_chunk1, event_counts);

    pthread_mutex_lock(&this->lock);

    if (doneflag) {
	pthread_mutex_unlock(&this->lock);
	throw runtime_error("ch_frb_io: internal error: doneflag already set in assembled_chunk_ringbuf::end_stream()");
    }

    // Wake up processing thread, if it is waiting for data
    pthread_cond_broadcast(&this->cond_assembled_chunks_added);

    this->doneflag = true;
    pthread_mutex_unlock(&this->lock);
}


std::unique_ptr<assembled_chunk> assembled_chunk_ringbuf::_make_assembled_chunk(uint64_t ichunk, int binning, bool zero)
{
    struct assembled_chunk::initializer chunk_params;

    chunk_params.beam_id = this->beam_id;
    chunk_params.nupfreq = this->nupfreq;
    chunk_params.nt_per_packet = this->nt_per_packet;
    chunk_params.fpga_counts_per_sample = this->fpga_counts_per_sample;
    chunk_params.force_reference = this->ini_params.force_reference_kernels;
    chunk_params.force_fast = this->ini_params.force_fast_kernels;
    chunk_params.binning = binning;
    chunk_params.ichunk = ichunk;

    if (ini_params.memory_pool) {
	chunk_params.pool = ini_params.memory_pool;
	chunk_params.slab = ini_params.memory_pool->get_slab(zero);

	if (!chunk_params.slab)
	    throw runtime_error("ch_frb_io: assembled chunk allocation failed!  FIXME: currently treated as an error, should add code to recover gracefully");
    }

    return assembled_chunk::make(chunk_params);
}


}  // namespace ch_frb_io
