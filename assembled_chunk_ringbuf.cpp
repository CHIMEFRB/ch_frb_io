#include <iostream>
#include "ch_frb_io_internals.hpp"

#include "l1-ringbuf.hpp"

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

#ifndef __AVX2__
    if (ini_params.mandate_fast_kernels)
	throw runtime_error("ch_frb_io: the 'mandate_fast_kernels' flag was set, but this machine does not have the AVX2 instruction set");
#endif

    ringbuf = new L1Ringbuf(beam_id_, ini_params.ringbuf_n);

    uint64_t packet_t0 = fpga_count0 / fpga_counts_per_sample;
    uint64_t ichunk = packet_t0 / constants::nt_per_assembled_chunk;

    this->active_chunk0 = this->_make_assembled_chunk(ichunk);
    this->active_chunk1 = this->_make_assembled_chunk(ichunk+1);

    pthread_mutex_init(&this->lock, NULL);
    pthread_cond_init(&this->cond_assembled_chunks_added, NULL);
}


assembled_chunk_ringbuf::~assembled_chunk_ringbuf()
{
    pthread_cond_destroy(&this->cond_assembled_chunks_added);
    pthread_mutex_destroy(&this->lock);
    delete ringbuf;
}

void assembled_chunk_ringbuf::print_state() {
    cout << "Beam " << beam_id << endl;
    ringbuf->print();
}

vector<pair<shared_ptr<assembled_chunk>, uint64_t> >
assembled_chunk_ringbuf::get_ringbuf_snapshot(uint64_t min_fpga_counts,
                                              uint64_t max_fpga_counts)
{
    vector<pair<shared_ptr<assembled_chunk>, uint64_t> > chunks;
    pthread_mutex_lock(&this->lock);
    ringbuf->retrieve(min_fpga_counts, max_fpga_counts, chunks);
    pthread_mutex_unlock(&this->lock);
    return chunks;
}

void assembled_chunk_ringbuf::get_ringbuf_size(uint64_t* ringbuf_fpga_next,
                                               uint64_t* ringbuf_n_ready,
                                               uint64_t* ringbuf_capacity,
                                               uint64_t* ringbuf_nelements,
                                               uint64_t* ringbuf_fpga_min,
                                               uint64_t* ringbuf_fpga_max) {
    pthread_mutex_lock(&this->lock);

    if (ringbuf_fpga_next) {
        shared_ptr<assembled_chunk> nxt = ringbuf->peek();
        if (!nxt)
            *ringbuf_fpga_next = 0;
        else
            *ringbuf_fpga_next = nxt->fpgacounts_begin();
    }

    if (ringbuf_n_ready)
        *ringbuf_n_ready = ringbuf->n_ready();

    if (ringbuf_capacity)
        *ringbuf_capacity = ringbuf->total_capacity();

    if (ringbuf_nelements)
        *ringbuf_nelements = ringbuf->total_size();

    if (ringbuf_fpga_min || ringbuf_fpga_max)
        ringbuf->fpga_counts_range(ringbuf_fpga_min, ringbuf_fpga_max);

    pthread_mutex_unlock(&this->lock);
}

// It's OK to modify 'event_counts' without acquiring any locks.  This is because the assembler
// thread passes an event_subcounts array which is updated on a per-packet basis, and accumulated
// into the global event_counts on a per-udp_packet_list basis (with locks acquired!).

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
        // after _put_assembled_chunk(), active_chunk0 has been reset.
        active_chunk0.swap(active_chunk1);
        // note that we've just swapped active_chunk1 down to active_chunk0, so active_chunk1's ichunk is active0 + 1
	active_chunk1 = this->_make_assembled_chunk(active_chunk0->ichunk + 1);
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

bool assembled_chunk_ringbuf::_put_assembled_chunk(unique_ptr<assembled_chunk> &chunk, int64_t *event_counts)
{
    if (!chunk)
	throw runtime_error("ch_frb_io: internal error: empty pointer passed to assembled_chunk_ringbuf::_put_unassembled_packet()");

    pthread_mutex_lock(&this->lock);

    if (this->doneflag) {
	pthread_mutex_unlock(&this->lock);
        chunk.reset();
	throw runtime_error("ch_frb_io: internal error: assembled_chunk_ringbuf::put_unassembled_packet() called after end_stream()");
    }

    if (stream_filename_pattern.length()) {
        string fn = chunk->format_filename(stream_filename_pattern);
        // turn on compression, but revert the state of chunk->msgpack_bitshuffle after writing.
        bool bitpack = chunk->msgpack_bitshuffle;
        chunk->msgpack_bitshuffle = true;
        chunk->write_msgpack_file(fn);
        chunk->msgpack_bitshuffle = bitpack;
    }

    // Convert unique_ptr into bare pointer, and reset unique_ptr.
    assembled_chunk* ch = chunk.release();
    // Try to add to ringbuf...
    shared_ptr<assembled_chunk> sch = ringbuf->push(ch);

    if (sch) {
	pthread_cond_broadcast(&this->cond_assembled_chunks_added);
	pthread_mutex_unlock(&this->lock);

        // Call any callbacks for chunks being added.
        for (auto it=chunk_callbacks.begin(); it!=chunk_callbacks.end(); it++) {
            (*it)(sch);
        }

        if (event_counts)
            event_counts[intensity_network_stream::event_type::assembled_chunk_queued]++;
	return true;
    }

    // If we get here, the ring buffer was full.
    pthread_mutex_unlock(&this->lock);
    if (event_counts)
        event_counts[intensity_network_stream::event_type::assembled_chunk_dropped]++;
    delete ch;

    if (ini_params.emit_warning_on_buffer_drop)
	cerr << "ch_frb_io: warning: processing thread is running too slow, dropping assembled_chunk\n";
    if (ini_params.throw_exception_on_buffer_drop)
	throw runtime_error("ch_frb_io: assembled_chunk was dropped and stream was constructed with 'throw_exception_on_buffer_drop' flag");
    return false;
}

bool assembled_chunk_ringbuf::inject_assembled_chunk(assembled_chunk* chunk) {
    uint64_t ich = chunk->ichunk;
    unique_ptr<assembled_chunk> uch(chunk);
    bool worked = _put_assembled_chunk(uch, NULL);
    // Danger: monkey with the active_chunk0, active_chunk1 variables,
    // which are not lock-protected and only supposed to be accessed
    // by the assembler thread.
    active_chunk0 = this->_make_assembled_chunk(ich + 1);
    active_chunk1 = this->_make_assembled_chunk(ich + 2);
    return worked;
}

shared_ptr<assembled_chunk> assembled_chunk_ringbuf::get_assembled_chunk(bool wait)
{
    pthread_mutex_lock(&this->lock);
    shared_ptr<assembled_chunk> chunk;

    for (;;) {
        chunk = ringbuf->pop();
        if (chunk)
            break;
        if (!wait)
            break;
	if (this->doneflag)
	    // Ring buffer is empty and end_stream() has been called.
            // (chunk contains null pointer)
            break;
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


std::unique_ptr<assembled_chunk> assembled_chunk_ringbuf::_make_assembled_chunk(uint64_t ichunk)
{
    bool force_ref  = false;
    bool force_fast = false;
    if (ini_params.mandate_fast_kernels)
        force_fast = true;
    else if (ini_params.mandate_reference_kernels)
        force_ref = true;
    return assembled_chunk::make(beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk, force_ref, force_fast);
}


}  // namespace ch_frb_io
