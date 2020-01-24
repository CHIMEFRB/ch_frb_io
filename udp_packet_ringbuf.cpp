#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif

// This is a lightweight scoped lock
typedef std::lock_guard<std::mutex> guard_t;
// This is also a scoped lock that supports use of a condition variable.
typedef std::unique_lock<std::mutex> ulock_t;

udp_packet_ringbuf::udp_packet_ringbuf(int ringbuf_capacity_, int max_npackets_per_list_, int max_nbytes_per_list_)
    : ringbuf_capacity(ringbuf_capacity_), 
      max_npackets_per_list(max_npackets_per_list_),
      max_nbytes_per_list(max_nbytes_per_list_)
{
    if (ringbuf_capacity <= 0)
	throw runtime_error("udp_packet_ringbuf constructor: expected ringbuf_capacity > 0");

    this->ringbuf.resize(ringbuf_capacity);
    for (int i = 0; i < ringbuf_capacity; i++)
	ringbuf[i] = make_unique<udp_packet_list> (this->max_npackets_per_list, this->max_nbytes_per_list);
}


void udp_packet_ringbuf::get_size(int* currsize, int* maxsize) {
    guard_t lock(mutx);
    if (currsize)
        *currsize = ringbuf_size;
    if (maxsize)
        *maxsize = ringbuf_capacity;
}

bool udp_packet_ringbuf::put_packet_list(unique_ptr<udp_packet_list> &p, bool is_blocking)
{    
    if (!p)
	throw runtime_error("ch_frb_io: udp_packet_ringbuf::put_packet_list() was called with empty pointer");

    ulock_t lock(mutx);

    for (;;) {
	if (stream_ended)
	    throw runtime_error("ch_frb_io: internal error: udp_packet_ringbuf::put_packet_list() called after end of stream");

	if (ringbuf_size < ringbuf_capacity) {
	    int i = (ringbuf_pos + ringbuf_size) % ringbuf_capacity;
	    std::swap(this->ringbuf[i], p);
	    this->ringbuf_size++;

            cond_packets_added.notify_all();
	    p->reset();
	    return true;
	}

	if (!is_blocking) {
	    p->reset();
	    return false;
	}

        cond_packets_removed.wait(lock);
    }
}


bool udp_packet_ringbuf::get_packet_list(unique_ptr<udp_packet_list> &p)
{
    if (!p)
	throw runtime_error("ch_frb_io: udp_packet_ringbuf::get_packet_list() was called with empty pointer");

    p->reset();

    ulock_t lock(mutx);

    for (;;) {
	if (ringbuf_size > 0) {
	    int i = ringbuf_pos % ringbuf_capacity;
	    std::swap(this->ringbuf[i], p);
	    this->ringbuf_pos++;
	    this->ringbuf_size--;

            cond_packets_removed.notify_all();
	    return true;
	}

	if (stream_ended)
	    return false;

        cond_packets_added.wait(lock);
    }
}


void udp_packet_ringbuf::end_stream()
{
    ulock_t lock(mutx);
    this->stream_ended = true;
    cond_packets_added.notify_all();
    cond_packets_removed.notify_all();
}


bool udp_packet_ringbuf::is_alive()
{
    guard_t lock(mutx);
    bool ret = !this->stream_ended;
    return ret;
}


}  // namespace ch_frb_io
