#include <fcntl.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>

#if defined(__linux__)
// for forking stats (SIOCOUTQ)
#include <linux/sockios.h>
#endif

#include <functional>
#include <algorithm>
#include <iostream>

#include <curl/curl.h>
#include <json/json.h>

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

// This is a lightweight scoped lock
typedef std::lock_guard<std::mutex> guard_t;
// This is also a scoped lock that supports use of a condition variable.
typedef std::unique_lock<std::mutex> ulock_t;


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
    packet_max_fpga_seen(0),
    network_thread_waiting_usec(0),
    network_thread_working_usec(0),
    assembler_thread_waiting_usec(0),
    assembler_thread_working_usec(0),
    frame0_nano(0),
    stream_priority(0),
    stream_chunks_written(0),
    stream_bytes_written(0),
    forking_paused(false)
{
    // Argument checking

    // FIXME I recently "promoted" a lot of compile-time constants to members of
    // intensity_network_stream::initializer.  Now that this has been done, it
    // would make sense to add a lot more checks here!

    if (ini_params.nbeams == 0)
	throw runtime_error("ch_frb_io: length-zero beam_id vector passed to intensity_network_stream constructor");

    if ((ini_params.nupfreq <= 0) || (ini_params.nupfreq > constants::max_allowed_nupfreq))
	throw runtime_error("ch_frb_io: bad value of 'nupfreq'");
    
    if (ini_params.nrfifreq < 0)
	throw runtime_error("ch_frb_io: bad value of 'nrfifreq'");
    
    if ((ini_params.nrfifreq > 0) && ((constants::nfreq_coarse_tot * ini_params.nupfreq) % ini_params.nrfifreq))
	throw runtime_error("ch_frb_io: bad value of 'nrfifreq'");

    if ((ini_params.nt_per_packet <= 0) || (ini_params.nt_per_packet > constants::max_allowed_nt_per_packet))
	throw runtime_error("ch_frb_io: bad value of 'nt_per_packet'");

    if ((ini_params.fpga_counts_per_sample <= 0) || (ini_params.fpga_counts_per_sample > constants::max_allowed_fpga_counts_per_sample))
	throw runtime_error("ch_frb_io: bad value of 'fpga_counts_per_sample'");	

    if ((ini_params.stream_id < 0) || (ini_params.stream_id > 9))
	throw runtime_error("ch_frb_io: bad value of 'stream_id'");

    if ((ini_params.nt_align < 0) || (ini_params.nt_align % constants::nt_per_assembled_chunk))
	throw runtime_error("ch_frb_io: 'nt_align' must be a multiple of nt_per_assembled_chunk(=" + to_string(constants::nt_per_assembled_chunk) + ")");
	 
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

    // Note: the socket is initialized in _open_socket().

    this->unassembled_ringbuf = make_unique<udp_packet_ringbuf> (ini_params.unassembled_ringbuf_capacity, 
								 ini_params.max_unassembled_packets_per_list, 
								 ini_params.max_unassembled_nbytes_per_list);

    this->incoming_packet_list = make_unique<udp_packet_list> (ini_params.max_unassembled_packets_per_list,
							       ini_params.max_unassembled_nbytes_per_list);

    this->cumulative_event_counts = vector<int64_t> (event_type::num_types, 0);
    this->network_thread_event_subcounts = vector<int64_t> (event_type::num_types, 0);
    this->assembler_thread_event_subcounts = vector<int64_t> (event_type::num_types, 0);

    network_thread_perhost_packets = make_shared<packet_counts>();
    perhost_packets = make_shared<packet_counts>();
}


intensity_network_stream::~intensity_network_stream()
{
    if (sockfd >= 0) {
	close(sockfd);
	sockfd = -1;
    }
}

void intensity_network_stream::add_first_packet_listener(first_packet_listener f) {
    first_packet_listeners.push_back(f);
}

vector<int> intensity_network_stream::get_beam_ids() {
    return beam_ids;
}

uint64_t intensity_network_stream::get_first_fpgacount() {
    if (!first_packet_received)
        throw runtime_error("ch_frb_io: get_first_fpgacount called, but first packet has not been received yet.");
    return first_fpgacount;
}

uint64_t intensity_network_stream::get_frame0_nano() {
    return frame0_nano;
}

shared_ptr<assembled_chunk_ringbuf> intensity_network_stream::_assembler_for_beam(int beam_id) {
    auto it = beam_to_assembler.find(beam_id);
    if (it == beam_to_assembler.end()) {
        return shared_ptr<assembled_chunk_ringbuf>();
    }
    return it->second;
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

    {
        guard_t lock(this->event_mutex);
        for (unsigned int i = 0; i < cumulative_event_counts.size(); i++)
            this->cumulative_event_counts[i] += event_subcounts[i];
    }

    memset(&event_subcounts[0], 0, event_subcounts.size() * sizeof(event_subcounts[0]));
}


void intensity_network_stream::start_stream()
{
    ulock_t lock(this->state_mutex);

    if (stream_end_requested || join_called)
	throw runtime_error("ch_frb_io: intensity_network_stream::start_stream() called on completed or cancelled stream");

    // If stream has already been started, this is not treated as an error.
    this->stream_started = true;
    this->cond_state_changed.notify_all();
}

void intensity_network_stream::wait_for_first_packet() {
    if (stream_end_requested || join_called)
	throw runtime_error("ch_frb_io: intensity_network_stream::wait_for_first_packet() called on completed or cancelled stream");
    ulock_t lock(this->state_mutex);
    // Wait for "stream_started"
    for (;;) {
	if (this->stream_end_requested)
	    // This case can arise if end_stream() is called early
	    return;
        if (this->first_packet_received)
            return;
	this->cond_state_changed.wait(lock);
    }
}

// Just sets the stream_end_requested flag and returns.  The shutdown logic then proceeeds as follows.
// The network thread will see that the stream_end_requested flag has been set, flush packets to the 
// assembler thread, call unassembled_ringbuf->end_stream(), and exit.  The assembler thread will see 
// that the ringbuf has ended, flush assembled_chunks to the processing threads, and exit.
//
// Note that end_stream() can be called multiple times (this usually happens as part of the shutdown process).

void intensity_network_stream::end_stream()
{
    ulock_t lock(this->state_mutex);
    this->stream_started = true;
    this->stream_end_requested = true;    
    this->cond_state_changed.notify_all();
}

void intensity_network_stream::flush_end_of_stream() {
    ulock_t lock(this->state_mutex);
    this->flush_end_of_stream_requested = true;
    this->cond_state_changed.notify_all();
}

void intensity_network_stream::reset_stream() {
    chlog("intensity_network_stream::reset_stream()!");
    {
        ulock_t lock(this->state_mutex);
        this->stream_started = false;
        this->first_packet_received = false;
        this->stream_end_requested = false;
        this->stream_restart = true;
        this->flush_end_of_stream_requested = true;
        
        chlog("i-n-s::reset_stream: joining assembler thread");
        this->assembler_thread.join();
        
        this->assemblers.clear();
        this->beam_ids.clear();
        this->beam_to_assembler.clear();
        this->first_fpgacount = 0;
        this->frame0_nano = 0;
        {
            ulock_t slock(this->stream_lock);
            this->stream_beam_ids.clear();
        }
        {
            ulock_t flock(this->forking_mutex);
            this->forking_packets.clear();
        }

        this->unassembled_ringbuf = make_unique<udp_packet_ringbuf> (ini_params.unassembled_ringbuf_capacity, 
                                                                     ini_params.max_unassembled_packets_per_list, 
                                                                     ini_params.max_unassembled_nbytes_per_list);
        
        chlog("i-n-s::reset_stream: restarting assembler thread");
        this->assembler_thread = std::thread(std::bind(&intensity_network_stream::assembler_thread_main, this));
        
        this->cond_state_changed.notify_all();

        // assembler thread quits when the unassembled_ringbuf runs dry
        // (happens when we receive an end-of-stream packet)
        
        /// shut down assembler thread
        /// reset first_fpgacount
        /// clear beam_to_assembmler
        /// reset frame0_nano ?
        /// reset event counts??
        // shut down streaming
        // shut down forking
        /// restart assembler thread
    }
}

void intensity_network_stream::join_threads()
{
    ulock_t lock(this->state_mutex);
    
    if (!stream_started)
	throw runtime_error("ch_frb_io: intensity_network_stream::join_threads() was called with no prior call to start_stream()");

    if (join_called) {
	while (!threads_joined)
            this->cond_state_changed.wait(lock);
	return;
    }

    this->join_called = true;
    this->cond_state_changed.notify_all();
    lock.unlock();

    network_thread.join();
    assembler_thread.join();

    lock.lock();
    this->threads_joined = true;
    this->cond_state_changed.notify_all();
}


void intensity_network_stream::stream_to_files(const string &filename_pattern, const vector<int> &stream_beam_ids, int priority, bool need_rfi)
{
    // Throw exception if 'beam_ids' argument contains a beam_id which is not
    // actually processed by this stream.
    //
    // Note: The 'nbeams', 'assemblers', and 'ini_params.beam_ids' fields are constant
    // after initialization, so we don't need to acquire a lock for these fields.

    unique_lock<mutex> ulock(stream_lock);

    for (int b: stream_beam_ids) {
        shared_ptr<assembled_chunk_ringbuf> assembler = _assembler_for_beam(b);
        if (!assembler)
	    throw runtime_error("intensity_network_stream::stream_to_files(): specified beam_id " + to_string(b) + " is not contained in by beam_ids: " + vstr(this->beam_ids));
    }
    
    this->stream_filename_pattern = filename_pattern;
    this->stream_beam_ids = stream_beam_ids;
    this->stream_priority = priority;
    this->stream_rfi_mask = need_rfi;
    this->stream_chunks_written = 0;
    this->stream_bytes_written = 0;

    // We're going to set the streaming status for *all* our assemblers -- clearing those
    // not in stream_beam_ids.
    for (int ibeam = 0; ibeam < (int)this->beam_ids.size(); ibeam++) {
	if (vcontains(stream_beam_ids, this->beam_ids[ibeam]))
	    assemblers[ibeam]->stream_to_files(filename_pattern, priority, need_rfi);
	else
	    assemblers[ibeam]->stream_to_files("", 0, false);
    }
}

void intensity_network_stream::get_streaming_status(std::string &filename_pattern,
                                                    std::vector<int> &beam_ids,
                                                    int &priority,
                                                    int &chunks_written,
                                                    size_t &bytes_written) {
    unique_lock<mutex> ulock(stream_lock);
    filename_pattern = this->stream_filename_pattern;
    beam_ids.insert(beam_ids.end(),
                    this->stream_beam_ids.begin(), this->stream_beam_ids.end());
    priority = this->stream_priority;

    chunks_written = 0;
    bytes_written = 0;
    for (int ibeam = 0; ibeam < (int)this->beam_ids.size(); ibeam++) {
	if (!vcontains(beam_ids, this->beam_ids[ibeam]))
            continue;
        int achunks = 0;
        size_t abytes = 0;
        assemblers[ibeam]->get_streamed_chunks(achunks, abytes);
        chunks_written += achunks;
        bytes_written += abytes;
    }
}

void intensity_network_stream::print_state() {
    cout << "Intensity network stream state:" << endl;
    for (auto it = assemblers.begin(); it != assemblers.end(); it++) {
        cout << "--Assembler:" << endl;
        (*it)->print_state();
    }
}


shared_ptr<assembled_chunk> intensity_network_stream::get_assembled_chunk(int assembler_index, bool wait)
{
    if ((assembler_index < 0) || (assembler_index >= (int)assemblers.size()))
        throw runtime_error("ch_frb_io: bad assembler_ix " + std::to_string(assembler_index) + " passed to intensity_network_stream::get_assembled_chunk() -- allowable range [0, " + std::to_string(assemblers.size()) + ")");

    auto ret = assemblers[assembler_index]->get_assembled_chunk(wait);

    // Note that we wait for data before crashing.
    if (ini_params.deliberately_crash)
	throw runtime_error("dedispersion thread will deliberately crash now!");

    return ret;
}


vector<int64_t> intensity_network_stream::get_event_counts()
{
    vector<int64_t> ret(event_type::num_types, 0);

    guard_t lock(this->event_mutex);
    memcpy(&ret[0], &this->cumulative_event_counts[0], ret.size() * sizeof(ret[0]));

    return ret;
}

unordered_map<string, uint64_t> intensity_network_stream::get_perhost_packets()
{
    // Quickly grab a copy of perhost_packets
    ulock_t lock(this->event_mutex);
    packet_counts pc(*perhost_packets);
    lock.unlock();
    return pc.to_string();
}

packet_counts::packet_counts() {}

packet_counts::packet_counts(const packet_counts& other) :
    tv(other.tv),
    period(other.period),
    counts(other.counts)
{}

unordered_map<string, uint64_t> packet_counts::to_string() const {
    // Convert to strings
    unordered_map<string, uint64_t> rtn;
    for (auto it = counts.begin(); it != counts.end(); it++) {
        // IPv4 address in high 32 bits, port in low 16 bits
        uint32_t ip = (it->first >> 32) & 0xffffffff;
        uint32_t port = (it->first & 0xffff);
        string sender = std::to_string(ip >> 24) + "." + std::to_string((ip >> 16) & 0xff)
            + "." + std::to_string((ip >> 8) & 0xff) + "." + std::to_string(ip & 0xff)
            + ":" + std::to_string(port);
        rtn[sender] = it->second;
    }
    return rtn;
}

void packet_counts::increment(const struct sockaddr_in& sender_addr,
                              int nbytes) {
    uint64_t ip = ntohl(sender_addr.sin_addr.s_addr);
    uint64_t port = ntohs(sender_addr.sin_port);
    // IPv4 address in high 32 bits, port in low 16 bits
    uint64_t sender = (ip << 32) | port;
    counts[sender]++;
}

void packet_counts::update(const packet_counts& other) {
    for (auto it=other.counts.begin(); it!=other.counts.end(); it++) {
        counts[it->first] = it->second;
    }
}

double packet_counts::start_time() const {
    return (double)tv.tv_sec + 1e-6 * (double)tv.tv_usec;
}

static void _get_history(double start, double end,
                         const map<double, shared_ptr<packet_counts> >& history,
                         vector<shared_ptr<packet_counts> >& counts) {
    if (history.size() == 0)
        return;
    // Most recent?
    if (start == 0) {
        auto it = history.end();
        it--;
        counts.push_back(it->second);
        return;
    }
    struct timeval now = xgettimeofday();
    double fnow = (double)now.tv_sec + 1e-6 * (double)now.tv_usec;
    // Negative values: N seconds ago.
    if (start <= 0.)
        start = fnow + start;
    if (end <= 0.)
        end = fnow + end;

    // find first history entry overlapping requested time period.
    // upper_bound() gives the first history entry > the start
    auto it = history.upper_bound(start);
    // the one before that, if it exists, should actually contain *start*.
    if (it == history.begin())
        return;
    it--;

    // likewise, upper_bound(end) gives us the history entry one past
    // the last one we want.
    //auto it2 = history.upper_bound(end);
    // for (; it != it2; it++)
    //counts.push_back(*it);
    
    for (; it!=history.end(); it++) {
        // Append until we find a packet that doesn't overlap *end*
        double t = it->second->start_time();
        if (t > end)
            break;
        counts.push_back(it->second);
    }
}

shared_ptr<packet_counts>
intensity_network_stream::get_packet_rates(double start, double period) {
    // This returns a single history entry.
    vector<shared_ptr<packet_counts> > counts;
    {
        guard_t lock(this->packet_history_mutex);
        _get_history(start, start+period, packet_history, counts);
    }
    // FIXME -- if *period* is specified, we could sum over the requested period...
    if (counts.size() == 0)
        return shared_ptr<packet_counts>();
    if (counts.size() == 1)
        return counts[0];
    shared_ptr<packet_counts> avg = make_shared<packet_counts>();
    // Just sum the samples, ignoring the possibly imperfect overlap between desired start+period and samples.
    avg->tv = counts[0]->tv;
    avg->period = 0.;
    for (auto it = counts.begin(); it != counts.end(); it++) {
        avg->period += (*it)->period;
        for (auto it2 = (*it)->counts.begin(); it2 != (*it)->counts.end(); it2++) {
            avg->counts[it2->first] += it2->second;
        }
    }
    return avg;
}

vector<shared_ptr<packet_counts> >
intensity_network_stream::get_packet_rate_history(double start, double end, double period) {
    vector<shared_ptr<packet_counts> > counts;
    guard_t lock(this->packet_history_mutex);
    _get_history(start, end, packet_history, counts);
    return counts;
}

void intensity_network_stream::fake_packet_from(const struct sockaddr_in& sender, int nbytes) {
    // network_thread_perhost_packets is normally only touched by the
    // network thread so is not lock-protected, but when updating the
    // history this is the lock used before reading the
    // network_thread_perhost_packets, so it should work...
    guard_t lock(this->event_mutex);
    network_thread_perhost_packets->increment(sender, nbytes);
    network_thread_perhost_packets->tv = xgettimeofday();

    cumulative_event_counts[event_type::packet_received] ++;
    cumulative_event_counts[event_type::packet_good] ++;
    cumulative_event_counts[event_type::byte_received] += nbytes;
}

vector<unordered_map<string, uint64_t> >
intensity_network_stream::get_statistics() {
    vector<unordered_map<string, uint64_t> > R;
    unordered_map<string, uint64_t> m;

    vector<int64_t> counts = get_event_counts();

    // Collect statistics for this stream as a whole:
    m["first_packet_received"]  = (counts[event_type::packet_received] > 0);
    m["nupfreq"]                = ini_params.nupfreq;
    m["nt_per_packet"]          = ini_params.nt_per_packet;
    m["fpga_counts_per_sample"] = ini_params.fpga_counts_per_sample;
    m["frame0_nano"]            = frame0_nano;
    m["fpga_count"]             = 0;    // XXX FIXME XXX
    m["network_thread_waiting_usec"] = network_thread_waiting_usec;
    m["network_thread_working_usec"] = network_thread_working_usec;
    m["assembler_thread_waiting_usec"] = assembler_thread_waiting_usec;
    m["assembler_thread_working_usec"] = assembler_thread_working_usec;

    m["count_bytes_received"     ] = counts[event_type::byte_received];
    m["count_packets_received"   ] = counts[event_type::packet_received];
    m["count_packets_good"       ] = counts[event_type::packet_good];
    m["count_packets_bad"        ] = counts[event_type::packet_bad];
    m["count_packets_dropped"    ] = counts[event_type::packet_dropped];
    m["count_packets_endofstream"] = counts[event_type::packet_end_of_stream];
    m["count_beam_id_mismatch"   ] = counts[event_type::beam_id_mismatch];
    m["count_stream_mismatch"    ] = counts[event_type::stream_mismatch];
    m["count_assembler_hits"     ] = counts[event_type::assembler_hit];
    m["count_assembler_misses"   ] = counts[event_type::assembler_miss];
    m["count_assembler_drops"    ] = counts[event_type::assembled_chunk_dropped];
    m["count_assembler_queued"   ] = counts[event_type::assembled_chunk_queued];

    int currsize, maxsize;
    unassembled_ringbuf->get_size(&currsize, &maxsize);
    m["udp_ringbuf_size"] = currsize;
    m["udp_ringbuf_maxsize"] = maxsize;

    m["count_bytes_queued"] = socket_queued_bytes;

    m["nbeams"] = ini_params.nbeams;

    // output_device status
    int nqueued = 0;
    for (size_t i=0; i<this->ini_params.output_devices.size(); i++) {
        string name = this->ini_params.output_devices[i]->ini_params.device_name;
        int chunks = this->ini_params.output_devices[i]->count_queued_write_requests();
        nqueued += chunks;
        for (size_t j=0; j<name.size(); j++)
            if ((name[j] == '/') || (name[j] == '-'))
                name[j] = '_';
        m["output_chunks_queued_" + name] = chunks;
    }
    m["output_chunks_queued"] = nqueued;
    
    // memory pool stats
    if (this->ini_params.memory_pool) {
        m["memory_pool_slab_nbytes"] = this->ini_params.memory_pool->nbytes_per_slab;
        m["memory_pool_slab_size"] = this->ini_params.memory_pool->nslabs;
        m["memory_pool_slab_avail"] = this->ini_params.memory_pool->count_slabs_available();
    } else {
        m["memory_pool_slab_nbytes"] = 0;
        m["memory_pool_slab_size"] = 0;
        m["memory_pool_slab_avail"] = 0;
    }
    
    // Streaming data to disk status
    {
        std::string streaming_filename_pattern;
        std::vector<int> streaming_beam_ids;
        int streaming_priority;
        int streaming_chunks_written;
        size_t streaming_bytes_written;
        this->get_streaming_status(streaming_filename_pattern, streaming_beam_ids,
                                   streaming_priority, streaming_chunks_written,
                                   streaming_bytes_written);
        m["streaming_n_beams"] = this->stream_beam_ids.size();
        for (unsigned int i=0; i<this->stream_beam_ids.size(); i++)
            m[stringprintf("streaming_beam_%i", i)] = this->stream_beam_ids[i];
        m["streaming_priority"] = this->stream_priority;
        m["streaming_bytes_written"] = this->stream_chunks_written;
        m["streaming_chunks_written"] = this->stream_bytes_written;
    }
    
    R.push_back(m);

    // Report per-host packet counts
    R.push_back(this->get_perhost_packets());

    // Collect statistics per beam:
    for (int b=0; b<this->beam_ids.size(); b++) {
        m.clear();
        m["beam_id"] = this->beam_ids[b];

        int streamed_chunks = 0;
        size_t streamed_bytes = 0;
        this->assemblers[b]->get_streamed_chunks(streamed_chunks, streamed_bytes);
        m["streaming_chunks_written"] = streamed_chunks;
        m["streaming_bytes_written"] = streamed_bytes;
        
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

        int lev;
        for (lev=1;; lev++) {
            this->assemblers[b]->get_ringbuf_size(NULL, NULL, &capacity, &nelements, &fpgacounts_min, &fpgacounts_max, lev);
            if (capacity == 0)
                break;
            m[stringprintf("ringbuf_capacity_level%i", lev)] = capacity;
            m[stringprintf("ringbuf_ntotal_level%i",   lev)] = nelements;
            m[stringprintf("ringbuf_fpga_min_level%i", lev)] = fpgacounts_min;
            m[stringprintf("ringbuf_fpga_max_level%i", lev)] = fpgacounts_max;
        }
        m["ringbuf_nlevels"] = lev;
        
        R.push_back(m);
    }

    return R;
}

std::shared_ptr<assembled_chunk>
intensity_network_stream::find_assembled_chunk(int beam, uint64_t fpga_counts, bool toplevel)
{
    // Which of my assemblers (if any) is handling the requested beam?
    shared_ptr<assembled_chunk_ringbuf> assembler = _assembler_for_beam(beam);
    if (!assembler)
        throw runtime_error("ch_frb_io internal error: beam_id mismatch in intensity_network_stream::find_assembled_chunk()");
    return assembler->find_assembled_chunk(fpga_counts, toplevel);
}

void intensity_network_stream::get_max_fpga_count_seen(vector<uint64_t> &flushed,
                                                       vector<uint64_t> &retrieved) {
    for (auto assembler : this->assemblers) {
        flushed.push_back(assembler->max_fpga_flushed);
        retrieved.push_back(assembler->max_fpga_retrieved);
    }
}

vector< vector< pair<shared_ptr<assembled_chunk>, uint64_t> > >
intensity_network_stream::get_ringbuf_snapshots(const vector<int> &beams,
                                                uint64_t min_fpga_counts,
                                                uint64_t max_fpga_counts)
{
    vector< vector< pair<shared_ptr<assembled_chunk>, uint64_t> > > R;
    if (beams.size()) {
        // Grab a snapshot for each requested beam (empty if we don't
        // have that beam number).
        for (int beam : beams) {
            // Which of my assemblers (if any) is handling the requested beam?
            shared_ptr<assembled_chunk_ringbuf> assembler = _assembler_for_beam(beam);
            if (assembler)
                R.push_back(assembler->get_ringbuf_snapshot(min_fpga_counts, max_fpga_counts));
            else
                // add empty list
                R.push_back(vector<pair<shared_ptr<assembled_chunk>, uint64_t > >());
        }
    } else {
        // Grab a snapshot from each of my assemblers.
        for (auto assembler : this->assemblers)
            R.push_back(assembler->get_ringbuf_snapshot(min_fpga_counts, max_fpga_counts));
    }
    return R;
}


// -------------------------------------------------------------------------------------------------
//
// Network thread


void intensity_network_stream::network_thread_main() 
{
    // We use try..catch to ensure that _network_thread_exit() always gets called, even if an exception is thrown.
    // We also print the exception so that it doesn't get "swallowed".

    chime_log_set_thread_name("Network-" + std::to_string(ini_params.stream_id));

    try {
	_network_thread_body();   // calls pin_thread_to_cores()
    } catch (exception &e) {
	cout << e.what() << endl;
	_network_thread_exit();
	throw;
    }
    _network_thread_exit();
}


void intensity_network_stream::_network_thread_body()
{
    pin_thread_to_cores(ini_params.network_thread_cores);
    ulock_t lock(this->state_mutex);

    // Wait for "stream_started"
    for (;;) {
	if (this->stream_end_requested)
	    // This case can arise if end_stream() is called early
	    return;
	if (this->stream_started) {
            lock.unlock();
	    break;
	}
	this->cond_state_changed.wait(lock);
    }
    // unlocked at this point.

    // Start listening on socket 

    string listening_msg = "ch_frb_io: listening for packets (ip_addr=" + ini_params.ipaddr + ", udp_port=" + to_string(ini_params.udp_port) + ")\n";
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

    cout << listening_msg;

    // Main packet loop

    for (;;) {
        _network_thread_one_stream();

        // wait for something to happen?!
        lock.lock();
        for (;;) {
            chlog("i_n_s network thread: waiting for state change!  join_called: " << join_called << ", stream_restart: " << stream_restart);
            if (join_called) {
                lock.unlock();
                return;
            }
            if (flush_end_of_stream_requested) {
                const int maxsize = 16384;
                uint8_t packet_data[maxsize];
                for (;;) {
                    int packet_nbytes = ::recv(sockfd, packet_data, maxsize, MSG_PEEK | MSG_DONTWAIT);
                    chlog("Flushing end-of-stream packets: peeked at a packet with " << packet_nbytes << " bytes.");
                    if (packet_nbytes == 24) {
                        packet_nbytes = ::recv(sockfd, packet_data, maxsize, 0);
                        chlog("dumped a packet with " << packet_nbytes << " bytes");
                    } else
                        break;
                }
                flush_end_of_stream_requested = false;
            }
            // !started means restarted
            //if (!stream_started) {
            if (stream_restart) {
                chlog("i-n-s network thread: found stream_restart!");
                stream_restart = false;
                lock.unlock();
                break;
            }
            this->cond_state_changed.wait(lock);
        }
    }
}

void intensity_network_stream::_network_thread_one_stream() {
    int64_t *event_subcounts = &this->network_thread_event_subcounts[0];
    struct timeval tv_ini = xgettimeofday();
    uint64_t packet_history_timestamp = 0;
    uint64_t incoming_packet_list_timestamp = 0;
    uint64_t cancellation_check_timestamp = 0;

    // All timestamps are in microseconds relative to tv_ini.
    uint64_t curr_timestamp = 0;
    struct timeval curr_tv = tv_ini;

    std::shared_ptr<packet_counts> last_packet_counts = make_shared<packet_counts>();
    // the previous counts added to the history list
    last_packet_counts->tv = tv_ini;

    //uint64_t rate_logging_timestamp = 0;
    //uint64_t rate_nbytes = 0;
    //uint64_t rate_npackets = 0;

    ulock_t lock(this->state_mutex);
    lock.unlock();

    for (;;) {
        uint64_t timestamp;

	// Periodically check whether stream has been cancelled by end_stream().
	if (curr_timestamp > cancellation_check_timestamp + ini_params.stream_cancellation_latency_usec) {
            lock.lock();
	    if (this->stream_end_requested) {
                lock.unlock();
                _network_flush_packets();
		return;
	    }
            lock.unlock();

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

        // Periodically store our per-sender packet counts
        if (curr_timestamp > packet_history_timestamp + ini_params.packet_count_period_usec) {
            _update_packet_rates(last_packet_counts);
            packet_history_timestamp = curr_timestamp;
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

        curr_tv = xgettimeofday();
        curr_timestamp = usec_between(tv_ini, curr_tv);
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
        }

        // Increment the number of packets we've received from this sender:
        network_thread_perhost_packets->increment(sender_addr, packet_nbytes);
        network_thread_perhost_packets->tv = curr_tv;

	event_subcounts[event_type::byte_received] += packet_nbytes;
	event_subcounts[event_type::packet_received]++;

	// If we receive a special "short" packet (length 24), it indicates end-of-stream.
	if (_unlikely(packet_nbytes == 24)) {
	    event_subcounts[event_type::packet_end_of_stream]++;
	    if (ini_params.accept_end_of_stream_packets) {
                chlog("i_n_s: Received end-of-stream packet.");
                // from network_thread_exit()...
                // This just sets the stream_end_requested flag, if it hasn't been set already.
                this->end_stream();
                // Flush any pending packets to assembler thread.
                this->_put_unassembled_packets();
                // Make sure all event counts are accumulated.
                this->_add_event_counts(network_thread_event_subcounts);
                // Set end-of-stream flag in the unassembled_ringbuf, so that the assembler knows there are no more packets.
                unassembled_ringbuf->end_stream();
                return;
            }
	    continue;
	}

	// The incoming_packet_list is timestamped with the arrival time of its first packet.
	if (incoming_packet_list->curr_npackets == 0)
	    incoming_packet_list_timestamp = curr_timestamp;

	incoming_packet_list->add_packet(packet_nbytes, sender_addr);

        //rate_nbytes += packet_nbytes;
        //rate_npackets++;

	if (incoming_packet_list->is_full) {
            _network_flush_packets();
            /*
            // Log our packet receive rate each flush!
            float dt = (float)(curr_timestamp - rate_logging_timestamp) / 1e6;
            chlogf("Packet receive rate: t %.3f : %.3f packets/sec, %.0f bits/sec",
                   curr_timestamp * 1e-6,
                   (float)rate_npackets / dt, (float)rate_nbytes * 8 / dt);
            rate_logging_timestamp = curr_timestamp;
            rate_nbytes = 0;
            rate_npackets = 0;
             */
        }
    }
}


 
// This gets called from the network thread to flush packets to the assembler threads.
void intensity_network_stream::_network_flush_packets() 
{
    this->_put_unassembled_packets();
    this->_add_event_counts(network_thread_event_subcounts);

    // Update the "perhost_packets" counter from "network_thread_perhost_packets"
    guard_t lock(this->event_mutex);
    perhost_packets->update(*network_thread_perhost_packets);
    perhost_packets->tv = network_thread_perhost_packets->tv;
}

// This gets called from the network thread to update the "perhost_packets" counter from "network_thread_perhost_packets".
void intensity_network_stream::_update_packet_rates(std::shared_ptr<packet_counts> last_packet_counts)
{
    std::shared_ptr<packet_counts> this_packet_counts;
    {
        guard_t lock(this->event_mutex);
        perhost_packets->update(*network_thread_perhost_packets);
        perhost_packets->tv = network_thread_perhost_packets->tv;
        // deep copy
        this_packet_counts = make_shared<packet_counts>(*perhost_packets);
    }
    // Build new packet_counts structure with differences vs last_
    shared_ptr<packet_counts> count_diff = make_shared<packet_counts>();
    count_diff->tv = this_packet_counts->tv;
    count_diff->period = usec_between(last_packet_counts->tv, this_packet_counts->tv) * 1e-6;
    for (auto it : this_packet_counts->counts) {
        auto it2 = last_packet_counts->counts.find(it.first);
        if (it2 == last_packet_counts->counts.end())
            // not found in last_packet_counts; new sender
            count_diff->counts[it.first] = it.second;
        else
            count_diff->counts[it.first] = it.second - it2->second;
    }
    last_packet_counts.reset();
    last_packet_counts.swap(this_packet_counts);

    // Add *count_diff* to the packet history
    {
        guard_t lock(this->packet_history_mutex);
        if (packet_history.size() >= (size_t)ini_params.max_packet_history_size)
            packet_history.erase(--packet_history.end());
        packet_history.insert(make_pair(count_diff->start_time(), count_diff));
    }
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

    // We use try..catch to ensure that _assembler_thread_exit() always gets called, even if an exception is thrown.
    // We also print the exception so that it doesn't get "swallowed".

    chime_log_set_thread_name("Assembler-" + std::to_string(ini_params.stream_id));

    try {
	_assembler_thread_body();  // calls pin_thread_to_cores()
    } catch (exception &e) {
	cout << e.what() << endl;
	_assembler_thread_exit();
	throw;
    }
    _assembler_thread_exit();
}

void intensity_network_stream::start_forking_packets(int beam, int destbeam, const struct sockaddr_in& dest) {
    unique_lock<mutex> ulock(forking_mutex);

    if (forking_socket == 0) {
        forking_socket = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
        if (forking_socket < 0)
            throw runtime_error(string("ch_frb_io: couldn't create udp socket for forking: ") + strerror(errno));
        int flags = fcntl(forking_socket, F_GETFD);
        flags |= FD_CLOEXEC;
        if (fcntl(forking_socket, F_SETFD, flags) < 0)
            throw runtime_error(string("ch_frb_io: couldn't set close-on-exec flag on packet-forking socket file descriptor") + strerror(errno));

        // bufsize
	int bufsize = 1024*1024;
	if (setsockopt(forking_socket, SOL_SOCKET, SO_SNDBUF, (void *) &bufsize, sizeof(bufsize)) < 0)
            throw runtime_error(string("ch_frb_io: setsockopt(SO_SNDBUF) failed for forking socket: ") + strerror(errno));
        // timeout
        const struct timeval tv_timeout = { 0, ini_params.socket_timeout_usec };
        if (setsockopt(forking_socket, SOL_SOCKET, SO_SNDTIMEO, &tv_timeout, sizeof(tv_timeout)) < 0)
            throw runtime_error(string("ch_frb_io: setsockopt(SO_SNDTIMEO) failed for forking socket: ") + strerror(errno));
    }
    intensity_network_stream::packetfork pf;
    // FIXME -- check that we actually hold the requested beam!!
    // FIXME -- check that the dest beam is valid!
    pf.beam = beam;
    pf.destbeam = destbeam;
    pf.dest = dest;
    forking_packets.push_back(pf);

    cout << "Forking packet list:" << endl;
    for (auto it = forking_packets.begin(); it != forking_packets.end(); it++)
        cout << "  beam " << it->beam << " -> " << it->destbeam << " to " <<
            inet_ntoa(it->dest.sin_addr) << " port " << ntohs(it->dest.sin_port) << endl;
}

void intensity_network_stream::stop_forking_packets(int beam, int destbeam, const struct sockaddr_in& dest) {

    // end-of-stream packet
    vector<uint8_t> packet(24, uint8_t(0));
    *((uint32_t *) &packet[0]) = uint32_t(1);  // protocol number

    unique_lock<mutex> ulock(forking_mutex);
    if ((beam == -1) && (destbeam == -1)) {
        // Stop all!
        for (auto& it : forking_packets) {
            if (forking_socket)
                sendto(forking_socket, &packet[0], packet.size(), 0,
                       reinterpret_cast<struct sockaddr*>(&(it.dest)), sizeof(it.dest));
        }
        forking_packets.clear();
    } else
        for (auto it = forking_packets.begin(); it != forking_packets.end(); it++)
            if ((it->beam == beam) && (it->destbeam == destbeam) &&
                (it->dest.sin_port == dest.sin_port) &&
                (it->dest.sin_addr.s_addr == dest.sin_addr.s_addr)) {

                if (forking_socket)
                    sendto(forking_socket, &packet[0], packet.size(), 0,
                           reinterpret_cast<struct sockaddr*>(&(it->dest)), sizeof(it->dest));

                forking_packets.erase(it);
                it--;
            }

    cout << "Forking packet list:" << endl;
    for (auto it = forking_packets.begin(); it != forking_packets.end(); it++)
        cout << "  beam " << it->beam << " -> " << it->destbeam << " to " <<
            inet_ntoa(it->dest.sin_addr) << " port " << ntohs(it->dest.sin_port) << endl;
}

void intensity_network_stream::pause_forking_packets() {
    forking_paused = true;
}

void intensity_network_stream::resume_forking_packets() {
    forking_paused = false;
}


// RAII
class forking_state {
public:
    forking_state(bool forking_active, int socket, mutex& mutx):
        active(forking_active),
        forking_socket(socket),
        flock(mutx, std::defer_lock)
    {
        if (!active)
            return;
        flock.lock();
#if defined(FIONWRITE)
        if (ioctl(forking_socket, FIONWRITE, &fork_sendqueue_start) == -1) {
            chlog("Failed to call ioctl(FIONWRITE): " << strerror(errno));
        }
#elif defined(SIOCOUTQ)
        // linux
        if (ioctl(forking_socket, SIOCOUTQ, &fork_sendqueue_start) == -1) {
            chlog("Failed to call ioctl(SIOCOUTQ): " << strerror(errno));
        }
#elif defined(SO_NWRITE)
        int sz = sizeof(int);
        if (getsockopt(forking_socket, SOL_SOCKET, SO_NWRITE,
                       &fork_sendqueue_start, reinterpret_cast<socklen_t*>(&sz))) {
            chlog("Failed to call getsockopt(SO_NWRITE): " << strerror(errno));
        }
#else
        fork_sendqueue_start = -1;
#endif
#if defined(FIONSPACE)
        if (ioctl(forking_socket, FIONSPACE, &fork_sendspace_start) == -1) {
            chlog("Failed to call ioctl(FIONSPACE): " << strerror(errno));
        }
#else
        fork_sendspace_start = -1;
#endif
    }

    void copy_packet(const intensity_packet& packet) {
        if (active)
            memcpy(&packetcopy, &packet, sizeof(intensity_packet));
    }
    void revert_packet(intensity_packet* packet) {
        if (active)
            memcpy(packet, &packetcopy, sizeof(intensity_packet));
    }

    void finish() {
#if defined(FIONWRITE)
        if (ioctl(forking_socket, FIONWRITE, &fork_sendqueue_end) == -1) {
            chlog("Failed to call ioctl(FIONWRITE): " << strerror(errno));
        }
#elif defined(SIOCOUTQ)
        if (ioctl(forking_socket, SIOCOUTQ, &fork_sendqueue_end) == -1) {
            chlog("Failed to call ioctl(SIOCOUTQ): " << strerror(errno));
        }
#elif defined(SO_NWRITE)
        int sz = sizeof(int);
        if (getsockopt(forking_socket, SOL_SOCKET, SO_NWRITE,
                       &fork_sendqueue_end, reinterpret_cast<socklen_t*>(&sz))) {
            chlog("Failed to call getsockopt(SO_NWRITE): " << strerror(errno));
        }
#else
        fork_sendqueue_end = -1;
#endif
#if defined(FIONSPACE)
        if (ioctl(forking_socket, FIONSPACE, &fork_sendspace_end) == -1) {
            chlog("Failed to call ioctl(FIONSPACE): " << strerror(errno));
        }
#else
        int isz = sizeof(int);
        int bufsize = 0;
        if (getsockopt(forking_socket, SOL_SOCKET, SO_SNDBUF,
                       &bufsize, reinterpret_cast<socklen_t*>(&isz))) {
            chlog("Failed to call getsockopt(SO_SNDBUF): " << strerror(errno));
        }
        fork_sendspace_end = bufsize;
        //fork_sendspace_end = -1;
#endif
        int esz = sizeof(int);
        int forking_error = 0;
        if (getsockopt(forking_socket, SOL_SOCKET, SO_ERROR,
                       &forking_error, reinterpret_cast<socklen_t*>(&esz))) {
            chlog("Failed to call getsockopt(SO_ERROR): " << strerror(errno));
        }
    }

    //chlog("Packet list: " << packet_list->curr_npackets << ", forwarded " << fork_packets_sent << " packets, " << fork_bytes_sent << " bytes in " << worktime/1000 << " ms / " << tperiod/1000 << " ms -> " << ((fork_packets_sent * 1e6) / tperiod) << " packets/sec, " << ((fork_bytes_sent * 1e6 * 8) / tperiod) << " bits/sec.  Send queue: " << fork_sendqueue_start << " | " << fork_sendspace_start << " at start, " << fork_sendqueue_end << " | " << fork_sendspace_end << " at end.  Socket error: " << forking_error);
    
    int fork_sendqueue_start = 0;
    int fork_sendqueue_end = 0;
    int fork_sendspace_start = 0;
    int fork_sendspace_end = 0;
    int fork_packets_sent = 0;
    int fork_bytes_sent = 0;

    intensity_packet packetcopy;

protected:
    bool active;
    int forking_socket;

    unique_lock<mutex> flock;
};


void intensity_network_stream::_assembler_thread_body()
{
    pin_thread_to_cores(ini_params.assembler_thread_cores);

    int nupfreq = this->ini_params.nupfreq;
    int nt_per_packet = this->ini_params.nt_per_packet;
    int fpga_counts_per_sample = this->ini_params.fpga_counts_per_sample;
    int nbeams = 0;

    auto packet_list = make_unique<udp_packet_list> (ini_params.max_unassembled_packets_per_list, ini_params.max_unassembled_nbytes_per_list);

    int64_t *event_subcounts = &this->assembler_thread_event_subcounts[0];

    uint8_t* forked_packet_data1 = reinterpret_cast<uint8_t*>(malloc(16384));
    uint8_t* forked_packet_data2 = reinterpret_cast<uint8_t*>(malloc(16384));
    uint8_t* forked_packet_data3 = reinterpret_cast<uint8_t*>(malloc(16384));
    uint8_t* forked_packet_data[3] = { forked_packet_data1,
                                       forked_packet_data2,
                                       forked_packet_data3 };

    struct timeval tva, tvb;
    tva = xgettimeofday();

    struct timeval tvf = tva;
    // Main packet loop

    while (1) {
        tvb = xgettimeofday();
        assembler_thread_working_usec += usec_between(tva, tvb);

        if (!unassembled_ringbuf->get_packet_list(packet_list)) {
            chlog("assembler thread: unassembled_ringbuf empty; quitting");
            break;
        }
        
        if (!first_packet_received) {
            // We just received our first packet.  Set beam ids!
            int ipacket = 0;
            uint8_t *packet_data = packet_list->get_packet_data(ipacket);
            int packet_nbytes = packet_list->get_packet_nbytes(ipacket);
            intensity_packet packet;
            packet.sender = packet_list->sender[ipacket];
            if (!packet.decode(packet_data, packet_nbytes)) {
                chlog("Bad packet (header, first packet) from " << ip_to_string(packet.sender));
                chlog("packet nbeams " << packet.nbeams << ", nf_c" << packet.nfreq_coarse <<
                      ", nupfreq " << packet.nupfreq << ", ntsamp " << packet.ntsamp << 
                      ", src_nb " << packet_nbytes);
                event_subcounts[event_type::packet_bad]++;
                continue;
            }

            if (this->ini_params.frame0_url.size()) {
                // After we receive our first packet, we will go fetch the frame0_ctime
                // via curl.  This is usually fast, so we'll do it in blocking mode.
                chlog("Retrieving frame0_ctime from " << this->ini_params.frame0_url);
                _fetch_frame0();    // raises runtime_error on failure
            }
            chlog("Received first packet.  Beams:" << packet.nbeams);
            for (int i=0; i<packet.nbeams; i++) {
                int beam = packet.beam_ids[i];
                chlog("  beam: " << beam);
                if ((beam < 0) || (beam > constants::max_allowed_beam_id))
                    throw runtime_error("ch_frb_io: bad beam_id received in first packet");
                auto assembler = make_shared<assembled_chunk_ringbuf>(ini_params, beam, ini_params.stream_id);
                assembler->set_frame0(frame0_nano);
                assemblers.push_back(assembler);
                beam_ids.push_back(beam);
                beam_to_assembler[beam] = assembler;
            }
            nbeams = this->beam_ids.size();
            // Compute first_fpgacount, used by rf_pipelines::chime_network_stream
            {
                // This is from assembled_chunk_ringbuf
                uint64_t packet_t0 = packet.fpga_count / packet.fpga_counts_per_sample;
                uint64_t first_ichunk = packet_t0 / constants::nt_per_assembled_chunk;
                if (ini_params.nt_align > 0) {
                    uint64_t chunk_align = ini_params.nt_align / constants::nt_per_assembled_chunk;
                    first_ichunk = ((first_ichunk + chunk_align - 1) / chunk_align) * chunk_align;
                }
                this->first_fpgacount = first_ichunk * constants::nt_per_assembled_chunk * ini_params.fpga_counts_per_sample;
            }
            
            {
                ulock_t lock(this->state_mutex);
                this->first_packet_received = true;
            }
            // you don't need to (and actually shouldn't) hold the lock when doing a notify.
            this->cond_state_changed.notify_all();
        }

        tva = xgettimeofday();
        assembler_thread_waiting_usec += usec_between(tvb, tva);

        bool forking_active = !forking_paused && (this->forking_packets.size() > 0);
        forking_state fs(forking_active, forking_socket, forking_mutex);

        for (int ipacket = 0; ipacket < packet_list->curr_npackets; ipacket++) {
            uint8_t *packet_data = packet_list->get_packet_data(ipacket);
            int packet_nbytes = packet_list->get_packet_nbytes(ipacket);
            intensity_packet packet;
            packet.sender = packet_list->sender[ipacket];

            if (!packet.decode(packet_data, packet_nbytes)) {
                chlog("Bad packet (header) from " << ip_to_string(packet.sender));
                chlog("packet nbeams " << packet.nbeams << ", nf_c" << packet.nfreq_coarse <<
                      ", nupfreq " << packet.nupfreq << ", ntsamp " << packet.ntsamp << 
                      ", src_nb " << packet_nbytes);
                event_subcounts[event_type::packet_bad]++;
                continue;
            }

            bool mismatch = ((packet.nbeams != ini_params.nbeams) ||
                             (packet.nupfreq != nupfreq) || 
                             (packet.ntsamp != nt_per_packet) || 
                             (packet.fpga_counts_per_sample != fpga_counts_per_sample));

            if (_unlikely(mismatch)) {
                if (ini_params.throw_exception_on_packet_mismatch) {
                    stringstream ss;
                    ss << "ch_frb_io: fatal: packet (nbeams, nupfreq, nt_per_packet, fpga_counts_per_sample) = ("
                       << packet.nbeams << "," << packet.nupfreq << "," << packet.ntsamp << "," << packet.fpga_counts_per_sample 
                       << "), expected ("
                       << ini_params.nbeams << "," << nupfreq << "," << nt_per_packet << "," << fpga_counts_per_sample << ")";
                    throw runtime_error(ss.str());
                }
                chlog("Packet stream mismatch from " << ip_to_string(packet.sender));
                event_subcounts[event_type::stream_mismatch]++;
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

            this->packet_max_fpga_seen = std::max(this->packet_max_fpga_seen.load(),
                                                  packet.fpga_count + (uint64_t)packet.ntsamp * (uint64_t)packet.fpga_counts_per_sample);

            int nfreq_coarse = packet.nfreq_coarse;
            int new_data_nbytes = nfreq_coarse * packet.nupfreq * packet.ntsamp;

            // Danger zone: we modify the packet by leaving its pointers in place, but shortening its
            // length fields.  The new packet corresponds to a subset of the original packet containing
            // only beam index zero.  This scheme avoids the overhead of copying the packet.

            // If we're re-sending packets, make a copy, because pointers get
            // modified in this loop!
            if (forking_active)
                fs.copy_packet(packet);

            packet.data_nbytes = new_data_nbytes;
            packet.nbeams = 1;

            for (int ibeam = 0; ibeam < nbeams; ibeam++) {
                // Loop invariant: at the top of this loop, 'packet' corresponds to a subset of the
                // original packet containing only beam index 'ibeam'.

                // Packet_id is beam number
                int packet_id = packet.beam_ids[0];

                shared_ptr<assembled_chunk_ringbuf> assembler = _assembler_for_beam(packet_id);
                if (_unlikely(!assembler)) {
                    // No match found
                    chlog("Beam id mismatch from " << ip_to_string(packet.sender));
                    event_subcounts[event_type::beam_id_mismatch]++;
                    if (ini_params.throw_exception_on_beam_id_mismatch)
                        throw runtime_error("ch_frb_io: beam_id mismatch occurred and stream was constructed with 'throw_exception_on_beam_id_mismatch' flag.  packet's beam_id: " + std::to_string(packet_id));
                } else {
                    // Match found
                    assembler->put_unassembled_packet(packet, event_subcounts);
                }

                // Danger zone: we do some pointer arithmetic, to modify the packet so that it now
                // corresponds to a new subset of the original packet, corresponding to beam index (ibeam+1).
                packet.beam_ids += 1;
                packet.scales += nfreq_coarse;
                packet.offsets += nfreq_coarse;
                packet.data += new_data_nbytes;
            }

            if (forking_active) {
                // Revert the packet header to its original state!
                fs.revert_packet(&packet);

                int nc = packet.nfreq_coarse;
                // Create headers for 1, 2, or 3-beam sub-packets.
                const int NS = 3;
                intensity_packet subpackets[NS];
                int nsubs[NS];
                for (int s=0; s<NS; s++) {
                    memcpy(&(subpackets[s]), &fs.packetcopy, sizeof(intensity_packet));
                    subpackets[s].nbeams = s+1;
                    subpackets[s].data_nbytes = subpackets[s].nbeams * nc * subpackets[s].nupfreq * subpackets[s].ntsamp;
                    nsubs[s] = subpackets[s].set_pointers(forked_packet_data[s]);
                }
                for (auto it : forking_packets) {
                    if (it.beam == 0) {
                        // send all (four) beams!
                        for (int i=0; i<packet.nbeams; i++)
                            packet.beam_ids[i] += it.destbeam;
                        int nsent = sendto(forking_socket, packet_data, packet_nbytes, 0,
                                           reinterpret_cast<struct sockaddr*>(&(it.dest)), sizeof(it.dest));
                        fs.fork_packets_sent++;
                        fs.fork_bytes_sent += nsent;
                        for (int i=0; i<packet.nbeams; i++)
                            packet.beam_ids[i] -= it.destbeam;
                        if (nsent == -1)
                            chlog("Failed to send forked packet data: " << strerror(errno));
                        if (nsent < packet_nbytes)
                            chlog("Sent " << nsent << " < " << packet_nbytes << " bytes of forked packet data!");
                        continue;
                    }
                    // send 1-3 beams
                    // start beam index
                    int ibeam = 0;
                    if (it.beam > 0) {
                        // find beam index
                        ibeam = -1;
                        for (int i=0; i<packet.nbeams; i++)
                            if (packet.beam_ids[i] == it.beam) {
                                ibeam = i;
                                break;
                            }
                        if (ibeam == -1) {
                            chlog("Forking: beam " << it.beam << " not found");
                            continue;
                        }
                    }
                    int nbeams = 1;
                    if (it.beam == -2) {
                        nbeams = 2;
                    } else if (it.beam == -3) {
                        nbeams = 3;
                    }
                    intensity_packet* subpacket = &(subpackets[nbeams-1]);
                    int nsub = nsubs[nbeams-1];
                    memcpy(subpacket->coarse_freq_ids,
                           packet.coarse_freq_ids,
                           nc * sizeof(uint16_t));
                    for (int i=0; i<nbeams; i++)
                        subpacket->beam_ids[i] = packet.beam_ids[ibeam + i] + it.destbeam;
                    memcpy(subpacket->scales,
                           packet.scales + ibeam * nc,
                           nbeams * nc * sizeof(float));
                    memcpy(subpacket->offsets,
                           packet.offsets + ibeam * nc,
                           nbeams * nc * sizeof(float));
                    memcpy(subpacket->data,
                           packet.data + ibeam * subpacket->data_nbytes/nbeams,
                           subpacket->data_nbytes);

                    int nsent = sendto(forking_socket, forked_packet_data[nbeams-1], nsub, 0,
                                       reinterpret_cast<struct sockaddr*>(&(it.dest)),
                                       sizeof(it.dest));
                    fs.fork_packets_sent++;
                    fs.fork_bytes_sent += nsent;
                    if (nsent == -1)
                        chlog("Failed to send forked packet data: " << strerror(errno));
                    if (nsent < nsub)
                        chlog("Sent " << nsent << " < " << nsub << " bytes of forked packet data!");
                }
            }
            if (forking_active)
                fs.finish();

            struct timeval tvnow = xgettimeofday();
            uint64_t worktime = usec_between(tva, tvnow);
            uint64_t tperiod = usec_between(tvf, tvnow);
            tvf = tvnow;
        }

	// We accumulate event counts once per udp_packet_list.
	this->_add_event_counts(assembler_thread_event_subcounts);
    }
    free(forked_packet_data1);
    free(forked_packet_data2);
    free(forked_packet_data3);
}

bool intensity_network_stream::inject_assembled_chunk(assembled_chunk* chunk) 
{
    // Find the right assembler and inject the chunk there.
    shared_ptr<assembled_chunk_ringbuf> assembler = _assembler_for_beam(chunk->beam_id);
    if (!assembler) {
        cout << "inject_assembled_chunk: no match for beam " << chunk->beam_id << endl;
        return false;
    }
    return assembler->inject_assembled_chunk(chunk);
}

// Called whenever the assembler thread exits (on all exit paths)
void intensity_network_stream::_assembler_thread_exit()
{
    this->end_stream();
    unassembled_ringbuf->end_stream();

    for (auto assembler: assemblers) {
	if (assembler)
	    assembler->end_stream(&assembler_thread_event_subcounts[0]);
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
       << "    stream mismatches: " << counts[event_type::stream_mismatch] << "\n"
       << "    assembler hits: " << counts[event_type::assembler_hit] << "\n"
       << "    assembler misses: " << counts[event_type::assembler_miss] << "\n"
       << "    assembled chunks dropped: " << counts[event_type::assembled_chunk_dropped] << "\n"
       << "    assembled chunks queued: " << counts[event_type::assembled_chunk_queued];

    cout << ss.str().c_str() << endl;
#endif
}


class CurlStringHolder {
public:
    string thestring;
};

static size_t
CurlWriteMemoryCallback(void *contents, size_t size, size_t nmemb, void *userp)
{
    size_t realsize = size * nmemb;
    CurlStringHolder* h = (CurlStringHolder*)userp;
    h->thestring += string((char*)contents, realsize);
    return realsize;
}

void intensity_network_stream::_fetch_frame0() {
    if (ini_params.frame0_url.size() == 0) {
        chlog("No 'frame0_url' set; skipping.");
        return;
    }
    CURL *curl_handle;
    CURLcode res;
    CurlStringHolder holder;
    // init the curl session
    curl_handle = curl_easy_init();
    // specify URL to get
    curl_easy_setopt(curl_handle, CURLOPT_URL, ini_params.frame0_url.c_str());
    // set timeout
    curl_easy_setopt(curl_handle, CURLOPT_TIMEOUT_MS, ini_params.frame0_timeout);
    // set received-data callback
    curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION,
                     CurlWriteMemoryCallback);
    curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, (void *)(&holder));
    // curl!
    chlog("Fetching frame0_time from " << ini_params.frame0_url);
    res = curl_easy_perform(curl_handle);
    if (res != CURLE_OK)
        throw runtime_error("ch_frb_io: fetch_frame0 failed: " + string(curl_easy_strerror(res)));
    curl_easy_cleanup(curl_handle);

    string frame0_txt = holder.thestring;
    //chlog("Received frame0 text: " << frame0_txt);
    Json::Reader frame0_reader;
    Json::Value frame0_json;
    if (!frame0_reader.parse(frame0_txt, frame0_json))
        throw runtime_error("ch_frb_io: failed to parse 'frame0' string: '" + frame0_txt + "'");

    //chlog("Parsed: " << frame0_json);
    if (!frame0_json.isObject())
        throw runtime_error("ch_frb_io: 'frame0' was not a JSON 'Object' as expected");

    string key = "frame0_nano";
    if (!frame0_json.isMember(key))
        throw runtime_error("ch_frb_io: 'frame0' did not contain key '" + key + "'");

    const Json::Value v = frame0_json[key];
    if (!v.isIntegral())
        throw runtime_error("ch_frb_io: expected 'frame0[frame0_nano]' to be integral.");
    frame0_nano = v.asUInt64();
    chlog("Found frame0_nano: " << frame0_nano);
}

}  // namespace ch_frb_io
