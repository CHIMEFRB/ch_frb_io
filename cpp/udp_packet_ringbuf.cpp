#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


udp_packet_ringbuf::udp_packet_ringbuf(int ringbuf_capacity_, int max_npackets_per_list_, int max_nbytes_per_list_)
    : ringbuf_capacity(ringbuf_capacity_), 
      max_npackets_per_list(max_npackets_per_list_),
      max_nbytes_per_list(max_nbytes_per_list_)
{
    if (ringbuf_capacity <= 0)
	throw runtime_error("udp_packet_ringbuf constructor: expected ringbuf_capacity > 0");

    this->ringbuf.resize(ringbuf_capacity);

    for (int i = 0; i < ringbuf_capacity; i++) {
	udp_packet_list l(max_npackets_per_list, max_nbytes_per_list);
	std::swap(this->ringbuf[i], l);
    }

    pthread_mutex_init(&this->lock, NULL);
    pthread_cond_init(&this->cond_packets_added, NULL);
    pthread_cond_init(&this->cond_packets_removed, NULL);
}


udp_packet_ringbuf::~udp_packet_ringbuf()
{
    pthread_mutex_destroy(&lock);
    pthread_cond_destroy(&cond_packets_added);
    pthread_cond_destroy(&cond_packets_removed);
}


bool udp_packet_ringbuf::put_packet_list(udp_packet_list &packet_list, bool is_blocking)
{    
    pthread_mutex_lock(&this->lock);

    for (;;) {
	if (stream_ended) {
	    pthread_mutex_unlock(&this->lock);
	    throw runtime_error("ch_frb_io: internal error: udp_packet_ringbuf::put_packet_list() called after end of stream");
	}

	if (ringbuf_size < ringbuf_capacity) {
	    int i = (ringbuf_pos + ringbuf_size) % ringbuf_capacity;
	    std::swap(this->ringbuf[i], packet_list);
	    this->ringbuf_size++;
	
	    pthread_cond_broadcast(&this->cond_packets_added);
	    pthread_mutex_unlock(&this->lock);
	    packet_list.reset();
	    return true;
	}

	if (!is_blocking) {
	    pthread_mutex_unlock(&this->lock);
	    packet_list.reset();
	    return false;
	}

	pthread_cond_wait(&this->cond_packets_removed, &this->lock);
    }
}


bool udp_packet_ringbuf::get_packet_list(udp_packet_list &packet_list)
{
    packet_list.reset();
    pthread_mutex_lock(&this->lock);

    for (;;) {
	if (ringbuf_size > 0) {
	    int i = ringbuf_pos % ringbuf_capacity;
	    std::swap(this->ringbuf[i], packet_list);
	    this->ringbuf_pos++;
	    this->ringbuf_size--;

	    pthread_cond_broadcast(&this->cond_packets_removed);
	    pthread_mutex_unlock(&this->lock);
	    return true;
	}

	if (stream_ended) {
	    pthread_mutex_unlock(&this->lock);
	    return false;
	}

	pthread_cond_wait(&this->cond_packets_added, &this->lock);
    }
}


void udp_packet_ringbuf::end_stream()
{
    pthread_mutex_lock(&this->lock);
    this->stream_ended = true;
    pthread_cond_broadcast(&this->cond_packets_added);
    pthread_cond_broadcast(&this->cond_packets_removed);
    pthread_mutex_unlock(&this->lock);
}


bool udp_packet_ringbuf::is_alive()
{
    pthread_mutex_lock(&this->lock);
    bool ret = !this->stream_ended;
    pthread_mutex_unlock(&this->lock);
    return ret;
}


}  // namespace ch_frb_io
