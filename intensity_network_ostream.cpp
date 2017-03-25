#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netdb.h>

#include <random>
#include <iostream>
#include "ch_frb_io_internals.hpp"

// linux-osx portability nuisance issue
#ifndef AI_V4MAPPED_CFG
#define AI_V4MAPPED_CFG AI_V4MAPPED
#endif


using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif



// -------------------------------------------------------------------------------------------------
//
// class intensity_network_ostream


// static member function
shared_ptr<intensity_network_ostream> intensity_network_ostream::make(const initializer &ini_params_)
{
    intensity_network_ostream *retp = new intensity_network_ostream(ini_params_);
    shared_ptr<intensity_network_ostream> ret(retp);

    ret->_open_socket();

    // Spawn network thread.  Note that we pass a bare pointer to an object ('ret') on our stack
    // and this pointer will be invalid after make() returns.  Therefore, the network thread only
    // dereferences the pointer before setting the network_thread_started flag, and make() waits for this
    // flag to be set before it returns.

    int err = pthread_create(&ret->network_thread, NULL, intensity_network_ostream::network_pthread_main, (void *) &ret);
    if (err < 0)
	throw runtime_error(string("ch_frb_io: pthread_create() failed in intensity_network_ostream constructor: ") + strerror(errno));

    pthread_mutex_lock(&ret->state_lock);
    while (!ret->network_thread_started)
	pthread_cond_wait(&ret->cond_state_changed, &ret->state_lock);
    pthread_mutex_unlock(&ret->state_lock);

    return ret;
}

    
intensity_network_ostream::intensity_network_ostream(const initializer &ini_params_) :
    ini_params(ini_params_),
    nbeams(ini_params.beam_ids.size()),
    nfreq_coarse_per_packet(ini_params.nfreq_coarse_per_packet),
    nfreq_coarse_per_chunk(ini_params.coarse_freq_ids.size()),
    nupfreq(ini_params.nupfreq),
    nt_per_packet(ini_params.nt_per_packet),
    nt_per_chunk(ini_params.nt_per_chunk),
    nbytes_per_packet(intensity_packet::packet_size(nbeams, nfreq_coarse_per_packet, nupfreq, nt_per_packet)),
    npackets_per_chunk((nfreq_coarse_per_chunk / nfreq_coarse_per_packet) * (nt_per_chunk / nt_per_packet)),
    nbytes_per_chunk(nbytes_per_packet * npackets_per_chunk),
    elts_per_chunk(nbeams * nfreq_coarse_per_chunk * nupfreq * nt_per_chunk),
    fpga_counts_per_sample(ini_params.fpga_counts_per_sample),
    fpga_counts_per_packet(fpga_counts_per_sample * nt_per_packet),
    fpga_counts_per_chunk(fpga_counts_per_sample * nt_per_chunk),
    target_gbps(ini_params.target_gbps)
{
    // Tons of argument checking.

    if ((ini_params.beam_ids.size() == 0) || (ini_params.beam_ids.size() >= 65536))
	throw runtime_error("chime intensity_network_ostream constructor: beam_ids vector is empty or too large");

    if ((ini_params.coarse_freq_ids.size() == 0) || (ini_params.coarse_freq_ids.size() >= 65536))
	throw runtime_error("chime intensity_network_ostream constructor: coarse_freq_ids vector is uninitialized or too large");
    if (nfreq_coarse_per_packet <= 0)
	throw runtime_error("chime intensity_network_ostream constructor: nfreq_coarse_per_packet was negative or uninitialized");
    if (nfreq_coarse_per_chunk % nfreq_coarse_per_packet != 0)
	throw runtime_error("chime intensity_network_ostream constructor: expected nfreq_coarse_per_chunk to be a multiple of nfreq_coarse_per_packet");

    if (nt_per_chunk <= 0)
	throw runtime_error("chime intensity_network_ostream constructor: nt_per_chunk was negative or uninitialized");
    if ((nt_per_packet <= 0) || (nt_per_packet > constants::max_allowed_nt_per_packet))
	throw runtime_error("chime intensity_network_ostream constructor: bad value of nt_per_packet (or uninitialized)");
    if (!is_power_of_two(nt_per_packet))
	throw runtime_error("chime intensity_network_ostream constructor: expected nt_per_packet to be a power of two");
    if (nt_per_chunk % nt_per_packet != 0)
	throw runtime_error("chime intensity_network_ostream constructor: expected nt_per_chunk to be a multiple of nt_per_packet");

    if ((nupfreq <= 0) || (nupfreq > constants::max_allowed_nupfreq))
	throw runtime_error("chime intensity_network_ostream constructor: bad value of nupfreq (or uninitalized)");
    if ((fpga_counts_per_sample <= 0) || (fpga_counts_per_sample > constants::max_allowed_fpga_counts_per_sample))
	throw runtime_error("chime intensity_network_ostream constructor: bad value of fpga_counts_per_sample (or uninitialized)");

    if (ini_params.wt_cutoff <= 0.0)
	throw runtime_error("chime intensity_network_ostream constructor: expected wt_cutoff to be > 0");
    if (ini_params.target_gbps < 0.0)
	throw runtime_error("ch_frb_io::intensity_network_ostream::initializer::target_gbps is negative");
    if (ini_params.target_gbps > constants::max_allowed_output_gbps)
	throw runtime_error("ch_frb_io::intensity_network_ostream::initializer::target_gbps is > max_allowed_output_gbps, presumably unintentional");
    if (!ini_params.throttle && (ini_params.target_gbps != 0.0))
	throw runtime_error("ch_frb_io::intensity_network_ostream::initializer::throttle is false, but target_gbps is nonzero, suspect this is a misconfiguration");

    if (nbytes_per_packet > constants::max_output_udp_packet_size)
	throw runtime_error("chime intensity_network_ostream constructor: packet size is too large, you need to decrease nfreq_per_packet or nt_per_packet");

    for (unsigned int i = 0; i < ini_params.beam_ids.size(); i++) {
	if ((ini_params.beam_ids[i] < 0) || (ini_params.beam_ids[i] > constants::max_allowed_beam_id))
	    throw runtime_error("intensity_network_ostream constructor: bad beam_id");
	for (unsigned int j = 0; j < i; j++)
	    if (ini_params.beam_ids[i] == ini_params.beam_ids[j])
		throw runtime_error("intensity_network_ostream constructor: duplicate beam_id");
    }

    for (unsigned int i = 0; i < ini_params.coarse_freq_ids.size(); i++) {
	if ((ini_params.coarse_freq_ids[i] < 0) || (ini_params.coarse_freq_ids[i] >= constants::nfreq_coarse_tot))
	    throw runtime_error("intensity_network_ostream constructor: bad coarse_freq_id");
	for (unsigned int j = 0; j < i; j++)
	    if (ini_params.coarse_freq_ids[i] == ini_params.coarse_freq_ids[j])
		throw runtime_error("intensity_network_ostream constructor: duplicate coarse_freq_id");
    }

    if (ini_params.throttle && (target_gbps == 0.0)) {
	// Infer target_gbps from fpga_counts_per_sample.
	double num = (8.0e-9 * nbytes_per_packet) * (nfreq_coarse_per_chunk / nfreq_coarse_per_packet);   // gigabits
	double den = (fpga_counts_per_packet * constants::dt_fpga);    // seconds
	this->target_gbps = num / den;
    }

    if (target_gbps > constants::max_allowed_output_gbps) {
	throw runtime_error("ch_frb_io::intensity_network_ostream:: inferred target_gbps(=" + to_string(target_gbps)
			    + ") is > max_allowed_output_gbps(=" + to_string(constants::max_allowed_output_gbps) + ")");
    }

    // Parse dstname (expect string of the form HOSTNAME[:PORT])

    this->hostname = ini_params.dstname;
    this->udp_port = constants::default_udp_port;

    size_t i = ini_params.dstname.find(":");

    if (i != std::string::npos) {
	this->hostname = ini_params.dstname.substr(0,i);
	if (!lexical_cast(ini_params.dstname.substr(i+1), this->udp_port))
	    throw runtime_error("ch_frb_io: couldn't convert string '" + ini_params.dstname.substr(i+1) + "' to 16-bit udp port number");
    }

    // Remaining initializations (except socket, which is initialized in intensity_network_ostream::_open_socket())

    this->beam_ids_16bit.resize(nbeams, 0);
    for (int i = 0; i < nbeams; i++)
	beam_ids_16bit[i] = uint16_t(ini_params.beam_ids[i]);

    this->coarse_freq_ids_16bit.resize(nfreq_coarse_per_chunk, 0);
    for (int i = 0; i < nfreq_coarse_per_chunk; i++)
	coarse_freq_ids_16bit[i] = uint16_t(ini_params.coarse_freq_ids[i]);
    
    pthread_mutex_init(&this->state_lock, NULL);
    pthread_mutex_init(&this->statistics_lock, NULL);
    pthread_cond_init(&this->cond_state_changed, NULL);

    int capacity = constants::output_ringbuf_capacity;
    this->ringbuf = make_unique<udp_packet_ringbuf> (capacity, npackets_per_chunk, nbytes_per_chunk);
    this->tmp_packet_list = make_unique<udp_packet_list> (npackets_per_chunk, nbytes_per_chunk);
}


intensity_network_ostream::~intensity_network_ostream()
{
    if (sockfd >= 0) {
	close(sockfd);
	sockfd = -1;
    }

    pthread_cond_destroy(&cond_state_changed);
    pthread_mutex_destroy(&state_lock);
    pthread_mutex_destroy(&statistics_lock);
}


// Socket initialization factored to its own routine, rather than putting it in the constructor,
// so that the socket will always be closed if an exception is thrown somewhere.
void intensity_network_ostream::_open_socket()
{
    if (sockfd >= 0)
	throw runtime_error("ch_frb_io: double call to intensity_network_ostream::_open_socket()");

    struct sockaddr_in saddr;
    memset(&saddr, 0, sizeof(saddr));
    saddr.sin_family = AF_INET;
    saddr.sin_port = htons(this->udp_port);

    struct addrinfo dns_hint;
    memset(&dns_hint, 0, sizeof(dns_hint));
    dns_hint.ai_flags = AI_V4MAPPED_CFG | AI_ADDRCONFIG;
    dns_hint.ai_family = AF_INET;   // IPv4
    dns_hint.ai_socktype = SOCK_DGRAM;
    dns_hint.ai_protocol = IPPROTO_UDP;

    struct addrinfo *dns_result = nullptr;

    int err = getaddrinfo(this->hostname.c_str(), NULL, &dns_hint, &dns_result);
    if (err)
	throw runtime_error("DNS lookup failed for '" + hostname + "': " + string(gai_strerror(err)));
    if (!dns_result)
	throw runtime_error("ch_frb_io: internal error: getaddrinfo() returned success, but result pointer is non-NULL");

    // just use first DNS entry returned (rather than traversing addrinfo::ai_next pointers)
    struct sockaddr_in *p = (struct sockaddr_in *) dns_result->ai_addr;
    saddr.sin_addr = p->sin_addr;

    freeaddrinfo(dns_result);

    this->sockfd = socket(AF_INET, SOCK_DGRAM, 0);
    if (sockfd < 0)
	throw runtime_error(string("ch_frb_io: couldn't create udp socket: ") + strerror(errno));

    int socket_bufsize = constants::send_socket_bufsize;
    err = setsockopt(sockfd, SOL_SOCKET, SO_SNDBUF, (void *) &socket_bufsize, sizeof(socket_bufsize));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: setsockopt(SO_SNDBUF) failed: ") + strerror(errno));
    
    if (ini_params.bind_port) {
        struct sockaddr_in server_address;
        memset(&server_address, 0, sizeof(server_address));
        server_address.sin_family = AF_INET;
        inet_pton(AF_INET, "0.0.0.0", &server_address.sin_addr);
        server_address.sin_port = htons(ini_params.bind_port);
        int err = ::bind(sockfd, (struct sockaddr *) &server_address, sizeof(server_address));
        if (err < 0)
            throw runtime_error(string("ch_frb_io: bind() failed: ") + strerror(errno));
    }

    if (connect(sockfd, reinterpret_cast<struct sockaddr *> (&saddr), sizeof(saddr)) < 0)
	throw runtime_error("ch_frb_io: couldn't connect udp socket to dstname '" + ini_params.dstname + "': " + strerror(errno));
}


// The 'intensity' and 'weights' arrays have shapes (nbeams, nfreq_coarse_per_chunk, nupfreq, nt_per_chunk)
void intensity_network_ostream::_encode_chunk(const float *intensity, const float *weights, int stride, uint64_t fpga_count, const unique_ptr<udp_packet_list> &out)
{
    // The number of packets per chunk is (nf_outer * nt_outer)
    int beam_stride = nfreq_coarse_per_chunk * nupfreq * stride;
    int nf_outer = nfreq_coarse_per_chunk / nfreq_coarse_per_packet;
    int nt_outer = nt_per_chunk / nt_per_packet;

    if (fpga_count % (fpga_counts_per_sample * nt_per_packet) != 0)
	throw runtime_error("intensity_network_ostream::_encode_chunk(): fpga count must be divisible by (fpga_counts_per_sample * nt_per_packet)");
    if (out->curr_npackets > 0)
	throw runtime_error("intensity_network_ostream::_encode_chunk(): internal error: packet_list nonempty");
    if (out->max_npackets < nt_outer * nf_outer)
	throw runtime_error("intensity_network_ostream::_encode_chunk(): internal error: packet_list is underallocated");
    if (out->max_nbytes < nt_outer * nf_outer * nbytes_per_packet)
	throw runtime_error("intensity_network_ostream::_encode_chunk(): internal error: packet_list is underallocated");

    intensity_packet packet;
    
    // Some intensity_packet fields are packet-independent; these are initialized here.
    packet.protocol_version = 1;
    packet.data_nbytes = nbeams * nfreq_coarse_per_packet * nupfreq * nt_per_packet;
    packet.fpga_counts_per_sample = fpga_counts_per_sample;
    packet.nbeams = nbeams;
    packet.nfreq_coarse = nfreq_coarse_per_packet;
    packet.nupfreq = nupfreq;
    packet.ntsamp = nt_per_packet;
    packet.beam_ids = &beam_ids_16bit[0];

    // Loop over packets in chunk.
    for (int it_outer = 0; it_outer < nt_outer; it_outer++) {
	for (int if_outer = 0; if_outer < nf_outer; if_outer++) {
	    int data_offset = (if_outer * nfreq_coarse_per_packet * nupfreq * stride) + (it_outer * nt_per_packet);

	    // Some intensity_packet fields are packet-dependent; these are initialized here.
	    packet.coarse_freq_ids = &coarse_freq_ids_16bit[if_outer * nfreq_coarse_per_packet];
	    packet.fpga_count = fpga_count + it_outer * nt_per_packet * fpga_counts_per_sample;

	    int nbytes_encoded = packet.encode(out->data_end, intensity + data_offset, weights + data_offset, beam_stride, stride, ini_params.wt_cutoff);

	    // A probably-paranoid sanity check
	    if (_unlikely(nbytes_encoded != nbytes_per_packet))
		throw runtime_error("ch_frb_io: internal error in network_ostream: nbytes_encoded != nbytes_per_packet");

	    out->add_packet(nbytes_per_packet);
	}
    }
}


void intensity_network_ostream::send_chunk(const float *intensity, const float *weights, int stride, uint64_t fpga_count)
{
    this->_encode_chunk(intensity, weights, stride, fpga_count, this->tmp_packet_list);

    if (ringbuf->put_packet_list(tmp_packet_list, ini_params.is_blocking))
	return;

    // If we get here, then packet list was dropped!
    if (ini_params.emit_warning_on_buffer_drop)
	cerr << "ch_frb_io: network write thread crashed or is running slow, dropping packets\n";
    if (ini_params.throw_exception_on_buffer_drop)
	throw runtime_error("ch_frb_io: packets were dropped and output stream was constructed with 'throw_exception_on_buffer_drop' flag");

}


void intensity_network_ostream::end_stream(bool join_network_thread)
{
    ringbuf->end_stream();

    if (!join_network_thread)
	return;

    pthread_mutex_lock(&this->state_lock);

    if (network_thread_joined) {
	pthread_mutex_unlock(&this->state_lock);
	throw runtime_error("ch_frb_io: attempt to join ostream output thread twice");
    }
    
    pthread_mutex_unlock(&this->state_lock);

    if (pthread_join(network_thread, NULL))
	throw runtime_error("ch_frb_io: couldn't join network thread [output]");
}


// -------------------------------------------------------------------------------------------------
//
// Network write thread


// static member function
void *intensity_network_ostream::network_pthread_main(void *opaque_arg)
{
    if (!opaque_arg)
	throw runtime_error("ch_frb_io: internal error: NULL opaque pointer passed to network_pthread_main()");

    // Note that the arg/opaque_arg pointer is only dereferenced here, for reasons explained in a comment in make() above.
    shared_ptr<intensity_network_ostream> *arg = (shared_ptr<intensity_network_ostream> *) opaque_arg;
    shared_ptr<intensity_network_ostream> stream = *arg;

    if (!stream)
	throw runtime_error("ch_frb_io: internal error: empty shared_ptr passed to network_pthread_main()");

    try {
	stream->_network_thread_body();
    } catch (...) {
	stream->end_stream(false);   // "false" means "don't join threads" (would deadlock otherwise!)
	throw;
    }

    stream->end_stream(false);   // "false" has same meaning as above
    return NULL;
}


void intensity_network_ostream::_network_thread_body()
{
    pthread_mutex_lock(&state_lock);
    network_thread_started = true;
    pthread_cond_broadcast(&cond_state_changed);
    pthread_mutex_unlock(&state_lock);

    auto packet_list = make_unique<udp_packet_list> (npackets_per_chunk, nbytes_per_chunk);
    
    // To be initialized when first packet is sent.
    struct timeval tv_ini;

    // Thread-local timestamps and counters.  
    //
    // Note that 'class intensity_network_ostream' defines three thread-shared members, which
    // will be kept in sync with their thread-local counterparts:
    //
    //    curr_timestamp <-> tstamp
    //    npackets_sent  <-> npackets_tot
    //    nbytes_sent    <-> nbytes_tot
    //
    // At the bottom of the packet loop below, we acquire the lock and update the thread-shared
    // versions.  In the body of the loop, we use the thread-locals so that we don't need
    // to acquire the lock repeatedly.  This scheme is safe because these variables are only 
    // modified by the network thread (i.e. other threads only access them read-only).

    int64_t prev_packet_nbytes = 0;
    int64_t prev_tstamp = 0;
    int64_t nbytes_tot = 0;     // not to be confused with thread-shared variable this->nbytes_sent
    int64_t npackets_tot = 0;   // not to be confused with this->npackets_sent
    int64_t tstamp = 0;         // not to be confused with this->curr_timestamp
    
    // Loop over packet_lists
    for (;;) {
	if (!ringbuf->get_packet_list(packet_list))
	    break;   // end of stream reached (probably normal termination)

	// Loop over packets
	for (int ipacket = 0; ipacket < packet_list->curr_npackets; ipacket++) {
	    const uint8_t *packet = packet_list->get_packet_data(ipacket);
	    const int packet_nbytes = packet_list->get_packet_nbytes(ipacket);

	    if (npackets_tot == 0)
		tv_ini = xgettimeofday();

	    tstamp = usec_between(tv_ini, xgettimeofday());

	    // Throttling logic: compare actual bandwidth to 'target_gbps' and sleep if necessary.
	    //
	    // t1 = "global" target timestamp based on cumulative transmission so far.  In normal
	    //      operation this will determine the time when the new packet is sent.
	    //
	    // t2 = "local" target timestamp based on transmission time of previous packet.  We
	    //      compute this below with a fudge factor of 0.8.  The idea is that in a situation
	    //      where the cumulative transfer has fallen behind, we can temporarily transmit
	    //      20% faster while we catch up.

	    if ((target_gbps > 0.0) && (npackets_tot > 0)) {
		int64_t t1 = 8.0e-3 * nbytes_tot / target_gbps;
		int64_t t2 = prev_tstamp + int64_t(0.8 * 8.0e-3 * prev_packet_nbytes / target_gbps);
		int64_t target_timestamp = min(t1,t2);

		if (tstamp < target_timestamp) {
		    xusleep(target_timestamp - tstamp);
		    tstamp = target_timestamp;
		}
	    }

            ssize_t n = this->_send(this->sockfd, packet, packet_nbytes, 0);
	    if (n < 0)
		throw runtime_error(string("chime intensity_network_ostream: udp packet send() failed: ") + strerror(errno));
	    if (n != packet_nbytes)
		throw runtime_error(string("chime intensity_network_ostream: udp packet send() sent ") + to_string(n) + "/" + to_string(packet_nbytes) + " bytes?!");

	    prev_packet_nbytes = packet_nbytes;
	    prev_tstamp = tstamp;
	    nbytes_tot += packet_nbytes;
	    npackets_tot++;

	    // Keep thread-shared variables in sync with thread-locals, as described above.
            pthread_mutex_lock(&statistics_lock);
	    this->curr_timestamp = tstamp;
	    this->nbytes_sent = nbytes_tot;
	    this->npackets_sent = npackets_tot;
            pthread_mutex_unlock(&statistics_lock);
	}
    }

    this->_announce_end_of_stream();
}

ssize_t intensity_network_ostream::_send(int socket, const uint8_t* packet, int nbytes, int flags) {
    return send(socket, packet, nbytes, flags);
}

void intensity_network_ostream::get_statistics(int64_t& curr_timestamp,
                                               int64_t& npackets_sent,
                                               int64_t& nbytes_sent) {
    pthread_mutex_lock(&statistics_lock);
    curr_timestamp = this->curr_timestamp;
    npackets_sent  = this->npackets_sent;
    nbytes_sent    = this->nbytes_sent;
    pthread_mutex_unlock(&statistics_lock);
}


void intensity_network_ostream::_announce_end_of_stream()
{
    // Send end-of-stream packets.  (This isn't part of the packet protocol, but the network _input_
    // stream contains an option to shut down gracefully if a special packet with nbeams=nupfreq=nt=0
    // is received.)
    //
    // Since UDP doesn't guarantee delivery, we have no way to ensure that the end-of-stream packet 
    // reaches the other side, but we'll make a best effort by sending 10 packets separated by 0.1 sec.

    cerr << "ch_frb_io: network output thread sending end-of-stream packets\n";

    for (int ipacket = 0; ipacket < 10; ipacket++) {
	vector<uint8_t> packet(24, uint8_t(0));
	*((uint32_t *) &packet[0]) = uint32_t(1);  // protocol number

	ssize_t n = send(this->sockfd, &packet[0], packet.size(), 0);

	if (n == (ssize_t)packet.size()) {
	    usleep(100000);  // 10^5 microseconds
	    continue;
	}

	// Emit warning if we fail on the first packet.  No warning emitted if subsequent packets
	// fail, since a likely explanation is that we're running over the loopback interface and
	// the receiving socket has been closed.

	if (ipacket == 0)
	    cerr << "warning: end-of-stream packets failed to send\n";

	break;
    }

    // End-of-stream summary info

    stringstream ss;
    ss << "ch_frb_io: network output thread exiting, npackets_sent=" << npackets_sent;
    if (npackets_sent >= 2)
	ss << ", gbps=" << (8.0e-3 * nbytes_sent / double(curr_timestamp));
    if (target_gbps > 0.0)
	ss << ", target_gbps=" << target_gbps;
    ss << "\n";

    cerr << ss.str();
}


}  // namespace ch_frb_io
