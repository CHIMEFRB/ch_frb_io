#include <fcntl.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>

#include <functional>
#include <iostream>

#include "ch_frb_io_internals.hpp"
#include "chlog.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


// -------------------------------------------------------------------------------------------------
//
// class intensity_network_stream

// Static member function (de facto constructor)
shared_ptr<intensity_network_stream> intensity_network_stream::make(const initializer &x)
{
    // not using make_shared because constructor is protected
    intensity_network_stream *retp = new intensity_network_stream(x);
    shared_ptr<intensity_network_stream> ret(retp);

    ret->_open_socket();

    // Spawn assembler thread.
    ret->assembler_thread = std::thread(std::bind(&intensity_network_stream::assembler_thread_main, ret));
    // Spawn network thread
    ret->network_thread = std::thread(std::bind(&intensity_network_stream::network_thread_main, ret));

    return ret;
}


intensity_network_stream::intensity_network_stream(const initializer &ini_params_) :
    ini_params(ini_params_),
    nassemblers(ini_params_.beam_ids.size()),
    network_thread_waiting_usec(0),
    network_thread_working_usec(0),
    assembler_thread_waiting_usec(0),
    assembler_thread_working_usec(0)
{
    // Argument checking

    // FIXME I recently "promoted" a lot of compile-time constants to members of
    // intensity_network_stream::initializer.  Now that this has been done, it
    // would make sense to add a lot more checks here!

    if (nassemblers == 0)
	throw runtime_error("ch_frb_io: length-zero beam_id vector passed to intensity_network_stream constructor");

    for (int i = 0; i < nassemblers; i++) {
	if ((ini_params.beam_ids[i] < 0) || (ini_params.beam_ids[i] > constants::max_allowed_beam_id))
	    throw runtime_error("ch_frb_io: bad beam_id passed to intensity_network_stream constructor");
	for (int j = 0; j < i; j++)
	    if (ini_params.beam_ids[i] == ini_params.beam_ids[j])
		throw runtime_error("ch_frb_io: duplicate beam_ids passed to intensity_network_stream constructor");
    }

    if ((ini_params.udp_port <= 0) || (ini_params.udp_port >= 65536))
	throw runtime_error("ch_frb_io: intensity_network_stream constructor: bad udp port " + to_string(ini_params.udp_port));

    if (ini_params.force_fast_kernels && ini_params.force_reference_kernels)
	throw runtime_error("ch_frb_io: both flags force_fast_kernels, force_reference_kernels were set");

#ifndef __AVX2__
    if (ini_params.force_fast_kernels)
	throw runtime_error("ch_frb_io: the 'force_fast_kernels' flag was set, but this machine does not have the AVX2 instruction set");
#endif

    if (ini_params.assembled_ringbuf_capacity <= 0)
	throw runtime_error("ch_frb_io: intensity_network_stream::initializer::assembled_ringbuf_capacity must be > 0");

    for (int n: ini_params.telescoping_ringbuf_capacity) {
	if (n < 2)
	    throw runtime_error("ch_frb_io: all elements of intensity_network_stream::initializer::telescoping_ringbuf_capacity must be >= 2");
    }

    // All initializations except the socket (which is initialized in _open_socket()),
    // and the assemblers (which are initalized when the first packet is received).

    this->unassembled_ringbuf = make_unique<udp_packet_ringbuf> (ini_params.unassembled_ringbuf_capacity, 
								 ini_params.max_unassembled_packets_per_list, 
								 ini_params.max_unassembled_nbytes_per_list);

    this->incoming_packet_list = make_unique<udp_packet_list> (ini_params.max_unassembled_packets_per_list,
							       ini_params.max_unassembled_nbytes_per_list);

    this->cumulative_event_counts = vector<int64_t> (event_type::num_types, 0);
    this->network_thread_event_subcounts = vector<int64_t> (event_type::num_types, 0);
    this->assembler_thread_event_subcounts = vector<int64_t> (event_type::num_types, 0);

    pthread_mutex_init(&state_lock, NULL);
    pthread_mutex_init(&event_lock, NULL);
    pthread_cond_init(&cond_state_changed, NULL);
}


intensity_network_stream::~intensity_network_stream()
{
    pthread_cond_destroy(&cond_state_changed);
    pthread_mutex_destroy(&state_lock);
    pthread_mutex_destroy(&event_lock);

    if (sockfd >= 0) {
	close(sockfd);
	sockfd = -1;
    }
}


// Socket initialization factored to its own routine, rather than putting it in the constructor,
// so that the socket will always be closed if an exception is thrown somewhere.
void intensity_network_stream::_open_socket()
{
    // FIXME assumes timeout < 1 sec
    const struct timeval tv_timeout = { 0, ini_params.socket_timeout_usec };

    this->sockfd = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    if (sockfd < 0)
	throw runtime_error(string("ch_frb_io: socket() failed: ") + strerror(errno));

    // In the CHIME L1 server, it was convenient to set the close-on-exec flag
    // on the socket file descriptor, to avoid corner cases such as a "zombie"
    // L1b process preventing the (ipaddr, port) pair being reused.

    int flags = fcntl(sockfd, F_GETFD);
    flags |= FD_CLOEXEC;

    if (fcntl(sockfd, F_SETFD, flags) < 0)
	throw runtime_error(string("ch_frb_io: couldn't set close-on-exec flag on socket file descriptor") + strerror(errno));

    // bufsize
    int err = setsockopt(sockfd, SOL_SOCKET, SO_RCVBUF, (void *) &ini_params.socket_bufsize, sizeof(ini_params.socket_bufsize));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: setsockopt(SO_RCVBUF) failed: ") + strerror(errno));

    // timeout
    err = setsockopt(sockfd, SOL_SOCKET, SO_RCVTIMEO, &tv_timeout, sizeof(tv_timeout));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: setsockopt(SO_RCVTIMEO) failed: ") + strerror(errno));
}


void intensity_network_stream::_add_event_counts(vector<int64_t> &event_subcounts)
{
    if (cumulative_event_counts.size() != event_subcounts.size())
	throw runtime_error("ch_frb_io: internal error: vector length mismatch in intensity_network_stream::_add_event_counts()");

    pthread_mutex_lock(&this->event_lock);
    for (unsigned int i = 0; i < cumulative_event_counts.size(); i++)
	this->cumulative_event_counts[i] += event_subcounts[i];
    pthread_mutex_unlock(&this->event_lock);

    memset(&event_subcounts[0], 0, event_subcounts.size() * sizeof(event_subcounts[0]));
}


void intensity_network_stream::start_stream()
{
    pthread_mutex_lock(&this->state_lock);

    if (stream_end_requested || join_called) {
	pthread_mutex_unlock(&this->state_lock);
	throw runtime_error("ch_frb_io: intensity_network_stream::start_stream() called on completed or cancelled stream");
    }

    // If stream has already been started, this is not treated as an error.
    this->stream_started = true;
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->state_lock);
}


// Just sets the stream_end_requested flag and returns.  The shutdown logic then proceeeds as follows.
// The network thread will see that the stream_end_requested flag has been set, flush packets to the 
// assembler thread, call unassembled_ringbuf->end_stream(), and exit.  The assembler thread will see 
// that the ringbuf has ended, flush assembled_chunks to the processing threads, and exit.
//
// Note that end_stream() can be called multiple times (this usually happens as part of the shutdown process).

void intensity_network_stream::end_stream()
{
    pthread_mutex_lock(&this->state_lock);
    this->stream_started = true;
    this->first_packet_received = true;
    this->assemblers_initialized = true;
    this->stream_end_requested = true;    
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->state_lock);
}


void intensity_network_stream::join_threads()
{
    pthread_mutex_lock(&this->state_lock);
    
    if (!stream_started) {
	pthread_mutex_unlock(&this->state_lock);
	throw runtime_error("ch_frb_io: intensity_network_stream::join_threads() was called with no prior call to start_stream()");
    }

    if (join_called) {
	while (!threads_joined)
	    pthread_cond_wait(&this->cond_state_changed, &this->state_lock);
	pthread_mutex_unlock(&this->state_lock);
	return;
    }

    this->join_called = true;
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->state_lock);

    network_thread.join();
    assembler_thread.join();

    pthread_mutex_lock(&this->state_lock);
    this->threads_joined = true;
    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->state_lock);    
}


void intensity_network_stream::stream_to_files(const std::string& filename_pattern) {
    stream_filename = filename_pattern;
    pthread_mutex_lock(&this->state_lock);
    bool inited = this->assemblers_initialized;
    pthread_mutex_unlock(&this->state_lock);
    if (inited)
        for (auto it = assemblers.begin(); it != assemblers.end(); it++)
            (*it)->stream_filename_pattern = filename_pattern;
}

void intensity_network_stream::add_assembled_chunk_callback(std::function<void(std::shared_ptr<assembled_chunk>)> callback) {
    pthread_mutex_lock(&this->state_lock);
    bool inited = this->assemblers_initialized;
    pthread_mutex_unlock(&this->state_lock);
    if (inited)
        for (auto it = assemblers.begin(); it != assemblers.end(); it++)
            (*it)->chunk_callbacks.push_back(callback);
    else
        chunk_callbacks.push_back(callback);
}

/*
 std::function has no sensible == operator, so how do we do this?!

void intensity_network_stream::remove_assembled_chunk_callback(std::function<void()> callback) {
    pthread_mutex_lock(&this->state_lock);
    bool inited = this->assemblers_initialized;
    pthread_mutex_unlock(&this->state_lock);
    if (inited)
        for (auto it = assemblers.begin(); it != assemblers.end(); it++) {
            std::vector<std::function<void()> > &vec = (*it)->chunk_callbacks;
            for (auto it2 = vec.begin(); it2 != vec.end(); it2++) {
                if ((*it2) == callback) {
                    vec.erase(it2);
                    break;
                }
            }
        }
    else {
        std::vector<std::function<void()> > &vec = chunk_callbacks;
        for (auto it2 = vec.begin(); it2 != vec.end(); it2++) {
            if ((*it2) == callback) {
                vec.erase(it2);
                break;
            }
        }
    }
}
 */

void intensity_network_stream::print_state() {
    cout << "Intensity network stream state:" << endl;
    for (auto it = assemblers.begin(); it != assemblers.end(); it++) {
        cout << "--Assembler:" << endl;
        (*it)->print_state();
    }
}

// Waits for the "assemblers_initialized" flag.
// If "prelocked = True", assume the state_lock mutex is already held,
// otherwise, it will be locked & unlocked on entry/exit.
void intensity_network_stream::_wait_for_assemblers_initialized(bool prelocked) {
    // wait for assemblers to be initialized...
    if (!prelocked)
        pthread_mutex_lock(&this->state_lock);

    while (!this->assemblers_initialized)
	pthread_cond_wait(&this->cond_state_changed, &this->state_lock);

    if (!prelocked)
        pthread_mutex_unlock(&this->state_lock);
}


shared_ptr<assembled_chunk> intensity_network_stream::get_assembled_chunk(int assembler_index, bool wait)
{
    if ((assembler_index < 0) || (assembler_index >= nassemblers))
	throw runtime_error("ch_frb_io: bad assembler_ix passed to intensity_network_stream::get_assembled_chunk()");

    // Wait for assemblers_initialized flag to be set.
    // After this, it's OK to access 'assemblers' without holding the lock.
    _wait_for_assemblers_initialized();

    // This can happen if end_stream() was called before the assemblers were initialized.
    if (assemblers.size() == 0)
	return shared_ptr<assembled_chunk> ();

    return assemblers[assembler_index]->get_assembled_chunk(wait);
}


intensity_network_stream::initializer intensity_network_stream::get_initializer()
{
    return this->ini_params;
}


vector<int64_t> intensity_network_stream::get_event_counts()
{
    vector<int64_t> ret(event_type::num_types, 0);

    pthread_mutex_lock(&this->event_lock);
    memcpy(&ret[0], &this->cumulative_event_counts[0], ret.size() * sizeof(ret[0]));
    pthread_mutex_unlock(&this->event_lock);    

    return ret;
}

unordered_map<string, uint64_t> intensity_network_stream::get_perhost_packets()
{
    // Quickly grab a copy of perhost_packets
    pthread_mutex_lock(&this->event_lock);
    unordered_map<uint64_t, uint64_t> raw(perhost_packets);
    pthread_mutex_unlock(&this->event_lock);

    // Convert to strings
    unordered_map<string, uint64_t> rtn;
    for (auto it = raw.begin(); it != raw.end(); it++) {
        // IPv4 address in high 32 bits, port in low 16 bits
        uint32_t ip = (it->first >> 32) & 0xffffffff;
        uint32_t port = (it->first & 0xffff);
        string sender = to_string(ip >> 24) + "." + to_string((ip >> 16) & 0xff)
            + "." + to_string((ip >> 8) & 0xff) + "." + to_string(ip & 0xff)
            + ":" + to_string(port);
        rtn[sender] = it->second;
    }
    return rtn;
}

vector<unordered_map<string, uint64_t> >
intensity_network_stream::get_statistics() {
    vector<unordered_map<string, uint64_t> > R;
    bool first_packet;
    bool assemblers_init;
    unordered_map<string, uint64_t> m;

    pthread_mutex_lock(&this->state_lock);
    first_packet = this->first_packet_received;
    assemblers_init = this->assemblers_initialized;
    pthread_mutex_unlock(&this->state_lock);

    // Collect statistics for this stream as a whole:
    m["first_packet_received"]  =  first_packet;
    m["nupfreq"]                = (first_packet ? this->fp_nupfreq : 0);
    m["nt_per_packet"]          = (first_packet ? this->fp_nt_per_packet : 0);
    m["fpga_counts_per_sample"] = (first_packet ? this->fp_fpga_counts_per_sample : 0);
    m["fpga_count"]             = (first_packet ? this->fp_fpga_count : 0);
    m["network_thread_waiting_usec"] = network_thread_waiting_usec;
    m["network_thread_working_usec"] = network_thread_working_usec;
    m["assembler_thread_waiting_usec"] = assembler_thread_waiting_usec;
    m["assembler_thread_working_usec"] = assembler_thread_working_usec;

    vector<int64_t> counts = get_event_counts();
    m["count_bytes_received"     ] = counts[event_type::byte_received];
    m["count_packets_received"   ] = counts[event_type::packet_received];
    m["count_packets_good"       ] = counts[event_type::packet_good];
    m["count_packets_bad"        ] = counts[event_type::packet_bad];
    m["count_packets_dropped"    ] = counts[event_type::packet_dropped];
    m["count_packets_endofstream"] = counts[event_type::packet_end_of_stream];
    m["count_beam_id_mismatch"   ] = counts[event_type::beam_id_mismatch];
    m["count_stream_mismatch"    ] = counts[event_type::first_packet_mismatch];
    m["count_assembler_hits"     ] = counts[event_type::assembler_hit];
    m["count_assembler_misses"   ] = counts[event_type::assembler_miss];
    m["count_assembler_drops"    ] = counts[event_type::assembled_chunk_dropped];
    m["count_assembler_queued"   ] = counts[event_type::assembled_chunk_queued];

    int currsize, maxsize;
    unassembled_ringbuf->get_size(&currsize, &maxsize);
    m["udp_ringbuf_size"] = currsize;
    m["udp_ringbuf_maxsize"] = maxsize;

    m["count_bytes_queued"] = socket_queued_bytes;

    int nbeams = this->ini_params.beam_ids.size();
    m["nbeams"] = nbeams;
    R.push_back(m);

    // Report per-host packet counts
    R.push_back(this->get_perhost_packets());

    // Collect statistics per beam:
    for (int b=0; b<nbeams; b++) {
        m.clear();
        m["beam_id"] = this->ini_params.beam_ids[b];
        if (assemblers_init && (assemblers.size() > 0)) {
            // Grab the ring buffer to find the min & max chunk numbers and size.
            uint64_t fpgacounts_next=0, n_ready=0, capacity=0, nelements=0, fpgacounts_min=0, fpgacounts_max=0;
            this->assemblers[b]->get_ringbuf_size(&fpgacounts_next, &n_ready, &capacity, &nelements, &fpgacounts_min, &fpgacounts_max);
            m["ringbuf_fpga_next"] = fpgacounts_next;
            m["ringbuf_n_ready"] = n_ready;
            m["ringbuf_capacity"] = capacity;
            m["ringbuf_ntotal"] = nelements;
            if (nelements == 0) {
                m["ringbuf_fpga_min"] = 0;
                m["ringbuf_fpga_max"] = 0;
            } else {
                m["ringbuf_fpga_min"] = fpgacounts_min;
                m["ringbuf_fpga_max"] = fpgacounts_max;
            }
        }
        R.push_back(m);
    }

    return R;
}

vector< vector< pair<shared_ptr<assembled_chunk>, uint64_t> > >
intensity_network_stream::get_ringbuf_snapshots(vector<uint64_t> &beams,
                                                uint64_t min_fpga_counts,
                                                uint64_t max_fpga_counts)
{
    vector< vector< pair<shared_ptr<assembled_chunk>, uint64_t> > > R;

    pthread_mutex_lock(&this->state_lock);
    bool assemblers_init = this->assemblers_initialized;
    pthread_mutex_unlock(&this->state_lock);

    if (!assemblers_init)
        return R;

    int nbeams = this->ini_params.beam_ids.size();

    for (size_t ib=0; ib<beams.size(); ib++) {
        uint64_t beam = beams[ib];
        vector<pair<shared_ptr<assembled_chunk>, uint64_t > > chunks;
        // Which of my assemblers (if any) is handling the requested beam?
        for (int i=0; i<nbeams; i++) {
	    if (this->ini_params.beam_ids[i] != (int)beam)
                continue;
            chunks = this->assemblers[i]->get_ringbuf_snapshot(min_fpga_counts,
                                                               max_fpga_counts);
        }
        R.push_back(chunks);
    }
    return R;
}

bool intensity_network_stream::get_first_packet_params(int &nupfreq, int &nt_per_packet, uint64_t &fpga_counts_per_sample, uint64_t &fpga_count)
{
    pthread_mutex_lock(&this->state_lock);

    while (!this->first_packet_received) {
	pthread_cond_wait(&this->cond_state_changed, &this->state_lock);
    
        if (this->stream_end_requested) {
            // end_stream() was called before first packet was received
            pthread_mutex_unlock(&this->state_lock);
            nupfreq = 0;
            nt_per_packet = 0;
            fpga_counts_per_sample = 0;
            fpga_count = 0;
            return false;
        }
    }
    
    pthread_mutex_unlock(&this->state_lock);

    // After the first_packet_received flag is set, it's OK to access the fp_* fields without holding the lock.
    nupfreq = fp_nupfreq;
    nt_per_packet = fp_nt_per_packet;
    fpga_counts_per_sample = fp_fpga_counts_per_sample;
    fpga_count = fp_fpga_count;
    return true;
}


// -------------------------------------------------------------------------------------------------
//
// Network thread


void intensity_network_stream::network_thread_main() {

    pin_thread_to_cores(ini_params.network_thread_cores);

    // We use try..catch to ensure that _network_thread_exit() always gets called, even if an exception is thrown.
    // We also print the exception so that it doesn't get "swallowed".

    try {
	_network_thread_body();
    } catch (exception &e) {
	cout << e.what() << endl;
	_network_thread_exit();
	throw;
    }
    _network_thread_exit();
}

void intensity_network_stream::_network_thread_body()
{
    pthread_mutex_lock(&this->state_lock);

    // Wait for "stream_started"
    for (;;) {
	if (this->stream_end_requested) {
	    // This case can arise if end_stream() is called early
	    pthread_mutex_unlock(&this->state_lock);
	    return;
	}
	if (this->stream_started) {
	    pthread_mutex_unlock(&this->state_lock);
	    break;
	}
	pthread_cond_wait(&this->cond_state_changed, &this->state_lock);
    }

    // Start listening on socket 

    string listening_msg = "ch_frb_io: listening for packets (ip_addr=" + ini_params.ipaddr + ", udp_port=" + to_string(ini_params.udp_port) + ")";
    string receiving_msg = "ch_frb_io: receiving packets! (ip_addr=" + ini_params.ipaddr + ", udp_port=" + to_string(ini_params.udp_port) + ")\n";

    struct sockaddr_in server_address;
    memset(&server_address, 0, sizeof(server_address));
    server_address.sin_family = AF_INET;
    server_address.sin_port = htons(ini_params.udp_port);

    int err = inet_pton(AF_INET, ini_params.ipaddr.c_str(), &server_address.sin_addr);
    if (err <= 0)
	throw runtime_error(ini_params.ipaddr + ": inet_pton() failed (note that no DNS lookup is done, the argument must be a numerical IP address)");

    err = ::bind(sockfd, (struct sockaddr *) &server_address, sizeof(server_address));
    if (err < 0)
	throw runtime_error(string("ch_frb_io: bind() failed (" + ini_params.ipaddr + ":" + to_string(ini_params.udp_port) + "): " + strerror(errno)));

    cout << listening_msg << endl;

    // Main packet loop

    int64_t *event_subcounts = &this->network_thread_event_subcounts[0];
    struct timeval tv_ini = xgettimeofday();
    uint64_t incoming_packet_list_timestamp = 0;
    uint64_t cancellation_check_timestamp = 0;
    bool is_first_packet = true;

    // All timestamps are in microseconds relative to tv_ini.
    uint64_t curr_timestamp = 0;

    for (;;) {
        uint64_t timestamp;

	// Periodically check whether stream has been cancelled by end_stream().
	if (curr_timestamp > cancellation_check_timestamp + ini_params.stream_cancellation_latency_usec) {
	    pthread_mutex_lock(&this->state_lock);

	    if (this->stream_end_requested) {
		pthread_mutex_unlock(&this->state_lock);    
                _network_flush_packets();
		return;
	    }

	    pthread_mutex_unlock(&this->state_lock);

	    // We call _add_event_counts() in a few different places in this routine, to ensure that
	    // the network thread's event counts are always regularly accumulated.
	    this->_add_event_counts(network_thread_event_subcounts);

	    cancellation_check_timestamp = curr_timestamp;
	}

	// Periodically flush packets to assembler thread (only happens if packet rate is low; normal case is that the packet_list fills first)
	if (curr_timestamp > incoming_packet_list_timestamp + ini_params.unassembled_ringbuf_timeout_usec) {
            _network_flush_packets();
	    incoming_packet_list_timestamp = curr_timestamp;
	}

	timestamp = usec_between(tv_ini, xgettimeofday());
        network_thread_working_usec += (timestamp - curr_timestamp);

	// Read new packet from socket (note that socket has a timeout, so this call can time out)
	uint8_t *packet_data = incoming_packet_list->data_end;
        // Record the sender IP & port here
        sockaddr_in sender_addr;
        int slen = sizeof(sender_addr);
        int packet_nbytes = ::recvfrom(sockfd, packet_data, ini_params.max_packet_size + 1, 0,
                                       (struct sockaddr *)&sender_addr, (socklen_t *)&slen);

	curr_timestamp = usec_between(tv_ini, xgettimeofday());
        network_thread_waiting_usec += (curr_timestamp - timestamp);

	// Check for error or timeout in read()
	if (packet_nbytes < 0) {
	    if ((errno == EAGAIN) || (errno == ETIMEDOUT))
		continue;  // normal timeout
	    if (errno == EINTR)
                continue; // this can happen when running in gdb
            throw runtime_error(string("ch_frb_io network thread: read() failed: ") + strerror(errno));
	}

        {
            int nqueued = 0;
            if (ioctl(sockfd, FIONREAD, &nqueued) == -1) {
                cout << "Failed to call ioctl(FIONREAD)" << endl;
            }
            socket_queued_bytes = nqueued;
            //cout << "recv: now " << nqueued << " bytes queued in UDP socket" << endl;
        }

        // Increment the number of packets we've received from this sender:
        uint64_t ip = ntohl(sender_addr.sin_addr.s_addr);
        uint64_t port = ntohs(sender_addr.sin_port);
        // IPv4 address in high 32 bits, port in low 16 bits
        uint64_t sender = (ip << 32) | port;
        network_thread_perhost_packets[sender]++;

	event_subcounts[event_type::byte_received] += packet_nbytes;
	event_subcounts[event_type::packet_received]++;

	// If we receive a special "short" packet (length 24), it indicates end-of-stream.
	if (_unlikely(packet_nbytes == 24)) {
	    event_subcounts[event_type::packet_end_of_stream]++;
	    if (ini_params.accept_end_of_stream_packets)
		return;   // triggers shutdown of entire stream
	    continue;
	}

	// Special logic for processing first packet, which triggers some initializations.
	if (is_first_packet) {
	    intensity_packet packet;

	    if (!packet.decode(packet_data, packet_nbytes)) {
		event_subcounts[event_type::packet_bad]++;
		continue;
	    }

	    this->fp_nupfreq = packet.nupfreq;
	    this->fp_nt_per_packet = packet.ntsamp;
	    this->fp_fpga_counts_per_sample = packet.fpga_counts_per_sample;
	    this->fp_fpga_count = packet.fpga_count;
	    is_first_packet = false;

	    // This will unblock the assembler thread.  (After being spawned, the assembler thread
	    // waits for the first_packet_received flag before doing anything else.)
	    pthread_mutex_lock(&this->state_lock);
	    this->first_packet_received = true;
	    pthread_cond_broadcast(&this->cond_state_changed);
	    pthread_mutex_unlock(&this->state_lock);

	    cout << receiving_msg;
	}

	// The incoming_packet_list is timestamped with the arrival time of its first packet.
	if (incoming_packet_list->curr_npackets == 0)
	    incoming_packet_list_timestamp = curr_timestamp;

	incoming_packet_list->add_packet(packet_nbytes);

	if (incoming_packet_list->is_full)
            _network_flush_packets();
    }
}

// This gets called from the network thread to flush packets to the assembler
// threads.
void intensity_network_stream::_network_flush_packets() {
    this->_put_unassembled_packets();
    this->_add_event_counts(network_thread_event_subcounts);

    // Update the "perhost_packets" counter from "network_thread_perhost_packets"
    pthread_mutex_lock(&this->event_lock);
    for (auto it = network_thread_perhost_packets.begin();
         it != network_thread_perhost_packets.end(); it++) {
        perhost_packets[it->first] = it->second;
    }
    pthread_mutex_unlock(&this->event_lock);
}

// This gets called when the network thread exits (on all exit paths).
void intensity_network_stream::_network_thread_exit()
{
    // This just sets the stream_end_requested flag, if it hasn't been set already.
    this->end_stream();
    
    // Flush any pending packets to assembler thread.
    this->_put_unassembled_packets();
    
    // Make sure all event counts are accumulated.
    this->_add_event_counts(network_thread_event_subcounts);

    // Set end-of-stream flag in the unassembled_ringbuf, so that the assembler knows there are no more packets.
    unassembled_ringbuf->end_stream();

    // Make sure socket is closed.
    if (sockfd >= 0) {
	close(sockfd);
	sockfd = -1;
    }
}


void intensity_network_stream::_put_unassembled_packets()
{
    int npackets = incoming_packet_list->curr_npackets;

    if (!npackets)
	return;

    bool success = unassembled_ringbuf->put_packet_list(incoming_packet_list, false);

    if (!success) {
	network_thread_event_subcounts[event_type::packet_dropped] += npackets;

	if (ini_params.emit_warning_on_buffer_drop)
	    cout << "ch_frb_io: assembler thread crashed or is running slow, dropping packets" << endl;
	if (ini_params.throw_exception_on_buffer_drop)
	    throw runtime_error("ch_frb_io: unassembled packets were dropped and stream was constructed with 'throw_exception_on_buffer_drop' flag");
    }
}


// -------------------------------------------------------------------------------------------------
//
// assembler thread


void intensity_network_stream::assembler_thread_main() {

    pin_thread_to_cores(ini_params.assembler_thread_cores);

    // We use try..catch to ensure that _assembler_thread_exit() always gets called, even if an exception is thrown.
    // We also print the exception so that it doesn't get "swallowed".

    try {
	_assembler_thread_body();
    } catch (exception &e) {
	cout << e.what() << endl;
	_assembler_thread_exit();
	throw;
    }
    _assembler_thread_exit();
}

void intensity_network_stream::_assembler_thread_body()
{
    pthread_mutex_lock(&this->state_lock);
    // wait for first_packet_received flag to be set
    for (;;) {
   	if (this->stream_end_requested) {
	    // This case can arise if end_stream() is called early
	    pthread_mutex_unlock(&this->state_lock);
	    return;
	}
	if (this->first_packet_received) {
	    pthread_mutex_unlock(&this->state_lock);
	    break;
	}
	pthread_cond_wait(&this->cond_state_changed, &this->state_lock);
    }

    uint16_t nupfreq = this->fp_nupfreq;
    uint16_t nt_per_packet = this->fp_nt_per_packet;
    uint16_t fpga_counts_per_sample = fp_fpga_counts_per_sample;

    // Allocate assemblers, and set the assemblers_initialized flag.  (Note: it would be simpler to do this
    // in the network thread, since we could coalesce the 'first_packet_received' and 'assemblers_initialized'
    // states into a single state, but this led to transient packet loss which was problematic for unit tests.)

    this->assemblers.resize(nassemblers);
    for (int ix = 0; ix < nassemblers; ix++) {
	assemblers[ix] = make_unique<assembled_chunk_ringbuf>(ini_params, ini_params.beam_ids[ix], nupfreq, nt_per_packet, fpga_counts_per_sample, this->fp_fpga_count);
    }

    pthread_mutex_lock(&this->state_lock);
    this->assemblers_initialized = true;

    // did a call to stream_to_file() come in while we were waiting
    // for the assemblers to be initialized?  If so, dispatch that
    // now.
    if (stream_filename.length()) {
        for (auto it = assemblers.begin(); it != assemblers.end(); it++)
            (*it)->stream_filename_pattern = stream_filename;
    }

    // did calls to add_assembled_chunk_callback come in while we were
    // waiting for the assembers?  If so, dispatch now.
    if (chunk_callbacks.size()) {
        for (auto it = assemblers.begin(); it != assemblers.end(); it++)
            (*it)->chunk_callbacks.insert((*it)->chunk_callbacks.end(),
                                          chunk_callbacks.begin(),
                                          chunk_callbacks.end());
        chunk_callbacks.clear();
    }

    pthread_cond_broadcast(&this->cond_state_changed);
    pthread_mutex_unlock(&this->state_lock);

    auto packet_list = make_unique<udp_packet_list> (ini_params.max_unassembled_packets_per_list, ini_params.max_unassembled_nbytes_per_list);

    int64_t *event_subcounts = &this->assembler_thread_event_subcounts[0];

    struct timeval tva, tvb;
    tva = xgettimeofday();

    // Main packet loop

    while (1) {
        tvb = xgettimeofday();
        assembler_thread_working_usec += usec_between(tva, tvb);

        if (!unassembled_ringbuf->get_packet_list(packet_list))
            break;

        tva = xgettimeofday();
        assembler_thread_waiting_usec += usec_between(tvb, tva);

	for (int ipacket = 0; ipacket < packet_list->curr_npackets; ipacket++) {
            uint8_t *packet_data = packet_list->get_packet_data(ipacket);
            int packet_nbytes = packet_list->get_packet_nbytes(ipacket);
	    intensity_packet packet;

	    if (!packet.decode(packet_data, packet_nbytes)) {
		event_subcounts[event_type::packet_bad]++;
		continue;
	    }

	    bool mismatch = (packet.nupfreq != nupfreq) || (packet.ntsamp != nt_per_packet) || (packet.fpga_counts_per_sample != fpga_counts_per_sample);

	    if (_unlikely(mismatch)) {
		event_subcounts[event_type::first_packet_mismatch]++;
		if (ini_params.throw_exception_on_first_packet_mismatch)
		    throw runtime_error("ch_frb_io: first-packet mismatch occurred and stream was constructed with 'throw_exception_on_first_packet_mismatch' flag");
		continue;
	    }

	    // All checks passed.  Packet is declared "good" here.  
	    //
	    // The following checks have been performed, either in this routine or in intensity_packet::read().
	    //   - dimensions (nbeams, nfreq_coarse, nupfreq, ntsamp) are positive,
	    //     and not large enough to lead to integer overflows
	    //   - packet and data byte counts are correct
	    //   - coarse_freq_ids are valid (didn't check for duplicates but that's ok)
	    //   - ntsamp is a power of two
	    //   - fpga_counts_per_sample is > 0
	    //   - fpga_count is a multiple of (fpga_counts_per_sample * ntsamp)
	    //
	    // These checks are assumed by assembled_chunk::add_packet(), and mostly aren't rechecked, 
	    // so it's important that they're done here!

	    event_subcounts[event_type::packet_good]++;
	    
	    int nbeams = packet.nbeams;
	    int nfreq_coarse = packet.nfreq_coarse;
	    int new_data_nbytes = nfreq_coarse * packet.nupfreq * packet.ntsamp;
	    const int *assembler_beam_ids = &ini_params.beam_ids[0];  // bare pointer for speed
    
	    // Danger zone: we modify the packet by leaving its pointers in place, but shortening its
	    // length fields.  The new packet corresponds to a subset of the original packet containing
	    // only beam index zero.  This scheme avoids the overhead of copying the packet.
	    
	    packet.data_nbytes = new_data_nbytes;
	    packet.nbeams = 1;
	    
	    for (int ibeam = 0; ibeam < nbeams; ibeam++) {
		// Loop invariant: at the top of this loop, 'packet' corresponds to a subset of the
		// original packet containing only beam index 'ibeam'.
		
		// Loop over assembler ids, to find a match for the packet_id.
		int packet_id = packet.beam_ids[0];
		int assembler_ix = 0;

		for (;;) {
		    if (assembler_ix >= nassemblers) {
			// No match found
			event_subcounts[event_type::beam_id_mismatch]++;
			if (ini_params.throw_exception_on_beam_id_mismatch)
			    throw runtime_error("ch_frb_io: beam_id mismatch occurred and stream was constructed with 'throw_exception_on_beam_id_mismatch' flag");
			break;
		    }

		    if (assembler_beam_ids[assembler_ix] != packet_id) {
			assembler_ix++;
			continue;
		    }

		    // Match found
		    assemblers[assembler_ix]->put_unassembled_packet(packet, event_subcounts);
		    break;
		}
		
		// Danger zone: we do some pointer arithmetic, to modify the packet so that it now
		// corresponds to a new subset of the original packet, corresponding to beam index (ibeam+1).
		
		packet.beam_ids += 1;
		packet.scales += nfreq_coarse;
		packet.offsets += nfreq_coarse;
		packet.data += new_data_nbytes;
	    }
	}

	// We accumulate event counts once per udp_packet_list.
	this->_add_event_counts(assembler_thread_event_subcounts);
    }
}

bool intensity_network_stream::inject_assembled_chunk(assembled_chunk* chunk) {

    pthread_mutex_lock(&this->state_lock);
    if (!first_packet_received) {
        cout << "Injecting assembled_chunk... initializing assembler threads." << endl;
        // Init...
        this->fp_nupfreq = chunk->nupfreq;
        this->fp_nt_per_packet = chunk->nt_per_packet;
        this->fp_fpga_counts_per_sample = chunk->fpga_counts_per_sample;
        this->fp_fpga_count = chunk->isample * chunk->fpga_counts_per_sample;
        this->first_packet_received = true;
        pthread_cond_broadcast(&this->cond_state_changed);

        // we already hold the state_lock
        _wait_for_assemblers_initialized(true);
    }
    pthread_mutex_unlock(&this->state_lock);

    // Find the right assembler and inject the chunk there.
    const int *assembler_beam_ids = &ini_params.beam_ids[0];
    for (int i = 0; i < nassemblers; i++) {
        if (assembler_beam_ids[i] == chunk->beam_id) {
            return assemblers[i]->inject_assembled_chunk(chunk);
        }
    }
    cout << "inject_assembled_chunk: no match for beam " << chunk->beam_id << endl;
    return false;
}

// Called whenever the assembler thread exits (on all exit paths)
void intensity_network_stream::_assembler_thread_exit()
{
    this->end_stream();
    unassembled_ringbuf->end_stream();

    // Use assemblers.size() instead of 'nassemblers', to correctly handle the case where end_stream()
    // gets called early, and the assemblers never get allocated.  (We also test for empty pointers, to
    // handle a corner case where an exception is thrown in the middle of the ring buffer allocation.)

    for (unsigned int i = 0; i < assemblers.size(); i++) {
	if (assemblers[i])
	    assemblers[i]->end_stream(&assembler_thread_event_subcounts[0]);
    }

    // Make sure all event counts are accumulated.
    this->_add_event_counts(assembler_thread_event_subcounts);

#if 0
    // Decided to remove this summary info!
    vector<int64_t> counts = this->get_event_counts();

    stringstream ss;
    ss << "ch_frb_io: assembler thread exiting\n"
       << "    bytes received (GB): " << (1.0e-9 * counts[event_type::byte_received]) << "\n"
       << "    packets received: " << counts[event_type::packet_received] << "\n"
       << "    good packets: " << counts[event_type::packet_good] << "\n"
       << "    bad packets: " << counts[event_type::packet_bad] << "\n"
       << "    dropped packets: " << counts[event_type::packet_dropped] << "\n"
       << "    end-of-stream packets: " << counts[event_type::packet_end_of_stream] << "\n"
       << "    beam id mismatches: " << counts[event_type::beam_id_mismatch] << "\n"
       << "    first-packet mismatches: " << counts[event_type::first_packet_mismatch] << "\n"
       << "    assembler hits: " << counts[event_type::assembler_hit] << "\n"
       << "    assembler misses: " << counts[event_type::assembler_miss] << "\n"
       << "    assembled chunks dropped: " << counts[event_type::assembled_chunk_dropped] << "\n"
       << "    assembled chunks queued: " << counts[event_type::assembled_chunk_queued];

    cout << ss.str().c_str() << endl;
#endif
}


}  // namespace ch_frb_io
