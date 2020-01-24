#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


static_assert(constants::max_input_udp_packet_size >= constants::max_output_udp_packet_size,
	      "udp_packet_list.cpp assumes constants::max_input_udp_packet_size >= constants::max_output_udp_packet_size");


udp_packet_list::udp_packet_list(int max_npackets_, int max_nbytes_) :
    max_npackets(max_npackets_), max_nbytes(max_nbytes_)
{
    if (max_npackets <= 0)
	throw runtime_error("udp_packet_list constructor: expected max_npackets > 0");
    if (max_nbytes <= 0)
	throw runtime_error("udp_packet_list constructor: expected max_nbytes > 0");

    // Note: this->curr_npackets and this->curr_nbytes are initialized to zero automatically.
    this->buf = unique_ptr<uint8_t[]> (new uint8_t[max_nbytes + constants::max_input_udp_packet_size]);
    this->off_buf = unique_ptr<int[]> (new int[max_npackets + 1]);
    this->sender = vector<sockaddr_in>(max_npackets);
    this->data_start = buf.get();
    this->data_end = data_start;
    this->packet_offsets = off_buf.get();
    this->packet_offsets[0] = 0;
}


void udp_packet_list::add_packet(int packet_nbytes, const sockaddr_in &sender)
{
    if (_unlikely((packet_nbytes <= 0) || (packet_nbytes > constants::max_input_udp_packet_size)))
	throw runtime_error("udp_packet_list::add_packet(): bad value of 'packet_nbytes'");
    if (_unlikely(is_full))
	throw runtime_error("udp_packet_list::add_packet() called on full packet_list");

    this->sender[this->curr_npackets] = sender;
    this->curr_npackets++;
    this->curr_nbytes += packet_nbytes;
    this->is_full = (curr_npackets >= max_npackets) || (curr_nbytes >= max_nbytes);
    this->data_end = data_start + curr_nbytes;
    this->packet_offsets[curr_npackets] = curr_nbytes;
}


void udp_packet_list::reset()
{
    this->curr_npackets = 0;
    this->curr_nbytes = 0;
    this->is_full = false;
    this->data_end = data_start;
    this->packet_offsets[0] = 0;
}


}  // namespace ch_frb_io
