#ifndef _CH_FRB_IO_HPP
#define _CH_FRB_IO_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++11 support (g++ -std=c++11)"
#endif

#include <queue>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <memory>
#include <cstdint>
#include <atomic>
#include <random>
#include <thread>
#include <iostream>
#include <mutex>
#include <condition_variable>

#include <hdf5.h>

#include <arpa/inet.h>

namespace ch_frb_io {
#if 0
}; // pacify emacs c-mode
#endif

template<typename T> struct hdf5_extendable_dataset;

struct noncopyable
{
    noncopyable() { }
    noncopyable(const noncopyable &) = delete;
    noncopyable& operator=(const noncopyable &) = delete;
};

// Defined later in this file
class assembled_chunk;
class memory_slab_pool;
class output_device;

// Defined in ch_frb_io_internals.hpp
struct intensity_packet;
struct udp_packet_list;
struct udp_packet_ringbuf;
class assembled_chunk_ringbuf;

// "uptr" is a unique_ptr for memory that is allocated by
// malloc()-like calls and should therefore be freed by free().
// It uses this little custom deleter class (that is
// default-constructable, so doesn't need to specified when creating a
// uptr).  This was copy-pasted from rf_pipelines.
struct uptr_deleter {
    inline void operator()(const void *p) { std::free(const_cast<void *> (p)); }
};
template<typename T>
using uptr = std::unique_ptr<T[], uptr_deleter>;
// And we use it for blocks of memory used for assembled_chunks.
typedef uptr<uint8_t> memory_slab_t;

// -------------------------------------------------------------------------------------------------
//
// Compile-time constants
//
// FIXME many of these don't really need to be fixed at compile time, and could be runtime parameters instead.


namespace constants {
    static constexpr int cache_line_size = 64;
    
    // Number of seconds per FPGA count.  This "magic number" appears in multiple libraries
    // (ch_frb_io, rf_pipelines, ch_vdif_assembler).  FIXME: it would be better to keep it
    // in one place, to avoid any possibility of having it get out-of-sync.
    static constexpr double dt_fpga = 2.56e-6;

    // Number of "coarse" (i.e. pre-upchannelized) frequency channels.
    static constexpr int nfreq_coarse_tot = 1024;

    // For an explanation of this parameter, see class intensity_network_ostream below 
    static constexpr double default_wt_cutoff = 0.3;

    static constexpr int default_udp_port = 6677;

#ifdef __APPLE__
    // osx seems to have very small limits on socket buffer size
    static constexpr int default_socket_bufsize = 4 * 1024 * 1024;
#else
    static constexpr int default_socket_bufsize = 128 * 1024 * 1024;
#endif

    // This applies to the ring buffer between the network _output_ thread, and callers of
    // intensity_network_ostream::send_chunk().
    static constexpr int output_ringbuf_capacity = 8;

    static constexpr int nt_per_assembled_chunk = 1024;

    // These parameters don't really affect anything but appear in asserts.
    static constexpr int max_input_udp_packet_size = 9000;   // largest value the input stream will accept
    static constexpr int max_output_udp_packet_size = 8910;  // largest value the output stream will produce
    static constexpr int max_allowed_beam_id = 65535;
    static constexpr int max_allowed_nupfreq = 64;
    static constexpr int max_allowed_nt_per_packet = 1024;
    static constexpr int max_allowed_fpga_counts_per_sample = 3200;
    static constexpr double max_allowed_output_gbps = 10.0;
};


// -------------------------------------------------------------------------------------------------
//
// HDF5 file I/O
//
// Note that there are two classes here: one for reading (intensity_hdf5_file),
// and one for writing (intensity_hdf5_ofile).


struct intensity_hdf5_file : noncopyable {
    std::string filename;
    
    int nfreq;
    int npol;
    int nt_file;     // number of time samples in file (which can have gaps)
    int nt_logical;  // total number of samples in time range spanned by file, including gaps

    //
    // We currently throw an exception unless the frequencies are equally spaced and consecutive.
    // Therefore, we don't keep a full list of frequencies, just the endpoints of the range and number of samples.
    //
    // This still leaves two possibilities: the frequencies can either be ordered from lowest to highest,
    // or vice versa.  The 'frequencies_are_increasing' flag is set in the former case.  (Note that the
    // standard CHIME convention is frequencies_are_increasing=false.)
    //
    // The 'freq_lo_MHz' and 'freq_hi_MHz' fields are the lowest and highest frequencies in the entire
    // band.  Just to spell this out in detail, if frequencies_are_increasing=true, then the i-th channel
    // spans the frequency range
    //
    //   [ freq_lo + i*(freq_hi-freq_lo)/nfreq, freq_lo + (i+1)*(freq_hi-freq_lo)/nfreq ]
    //
    // and if frequences_are_increasing=false, then the i-th channel spans frequency range
    //
    //   [ freq_hi - (i+1)*(freq_hi-freq_lo)/nfreq, freq_hi - i*(freq_hi-freq_lo)/nfreq ]
    //
    // where 0 <= i <= nfreq-1.
    //
    bool frequencies_are_increasing;
    double freq_lo_MHz;
    double freq_hi_MHz;
    
    //
    // We distinguish between "file" time indices, which span the range 0 <= it < nt_file,
    // and "logical" time indices, which span the range 0 <= it < nt_logical.
    //
    // The i-th logical time sample spans the time range
    //
    //    [ time_lo + i*dt_sample, time_lo + (i+1)*dt_sample ]
    //
    double dt_sample;
    double time_lo;
    double time_hi;   // always equal to time_lo + nt_logical * dt_sample

    std::vector<double> times;                // 1D array of length nt_file
    std::vector<int> time_index_mapping;      // 1D array of length nt_file, which maps a "file" index to a "logical" index.

    // Polarization info (currently read from file but not really used)
    std::vector<std::string> polarizations;  // 1D array of length npol, each element is either "XX" or "YY"

    // 3d arrays of shape (nfreq, npol, nt_file), verbatim from the hdf5 file.
    // Rather than using them directly, you may want to use the member function get_unpolarized_intensity() below.
    std::vector<float> intensity;
    std::vector<float> weights;

    // Summary statistics
    double frac_ungapped;    // fraction of "logical" time samples which aren't in time gaps
    double frac_unmasked;    // fraction of _ungapped_ data with large weight

    // Construct from file.  If 'noisy' is true, then a one-line message will be printed when the file is read.
    explicit intensity_hdf5_file(const std::string &filename, bool noisy=true);

    //
    // Extracts a 2D array containing total intensity in time range [out_t0, out_t0+out_nt),
    // summing over polarizations.  The 'out_t0' and 'out_nt' are "logical" time indices, not
    // "file" time indices.  If this range of logical time indices contains gaps, the corresponding
    // entries of the 'out_int' and 'out_wt' arrays will be filled with zeros.
    // 
    // The 'out_int' and 'out_wt' arrays have shape (nfreq, out_nt).
    //
    // The 'out_istride' and 'out_wstride' args can be negative, if reversing the channel ordering is desired.
    // If out_istride is zero, it defaults to out_nt.  If out_wstride is zero, it defaults to out_istride.
    //
    void get_unpolarized_intensity(float *out_int, float *out_wt, int out_t0, int out_nt, int out_istride=0, int out_wstride=0) const;

    void run_unit_tests() const;
};


struct intensity_hdf5_ofile {
    std::string filename;
    double dt_sample;
    int nfreq;
    int npol;

    ssize_t curr_nt;       // current size of file (in time samples, not including gaps)
    double curr_time;      // time in seconds relative to arbitrary origin
    ssize_t curr_ipos;     // keeps track of gaps
    ssize_t initial_ipos;

    // used internally to print summary info when file is written
    double wsum;
    double wmax;

    std::unique_ptr<hdf5_extendable_dataset<double> > time_dataset;
    std::unique_ptr<hdf5_extendable_dataset<float> > intensity_dataset;
    std::unique_ptr<hdf5_extendable_dataset<float> > weights_dataset;

    //
    // The 'pol' argument is typically { "XX", "YY" }.  Note that pol.size() determines intensity_hdf5_file::npol,
    // which in turn determines the expected shape of the 'intensity' and 'weights' arguments passed to append_chunk().
    //
    // The 'freq0_MHz' and 'freq1_MHz' args should be the edges of the band, ordered the same way as the channel
    // indices.  For example in CHIME data, frequencies are ordered from highest to lowest, so we take freq0_MHz=800
    // and freq1_MHz=400.
    //
    // The optional 'ipos0' and 'time0' args are:
    //   ipos0 = index of first sample in file (in downsampled units, i.e. one sample is ~1 msec, not ~2.56 usec)
    //   time0 = arrival time of first sample in file (in seconds).
    //
    // The meaning of the 'bitshuffle' arg is:
    //   0 = no compression
    //   1 = try to compress, but if plugin fails then just write uncompressed data instead
    //   2 = try to compress, but if plugin fails then print a warning and write uncompressed data instead
    //   3 = compression mandatory
    //
    // The default nt_chunk=128 comes from ch_vdif_assembler chunk size, assuming downsampling by factor 512.
    //
    intensity_hdf5_ofile(const std::string &filename, int nfreq, const std::vector<std::string> &pol,
			 double freq0_MHz, double freq1_MHz, double dt_sample, ssize_t ipos0=0,
			 double time0=0.0, int bitshuffle=2, int nt_chunk=128);

    // Note that there is no close() member function.  The file is flushed to disk and closed when the
    // intensity_hdf5_ofile destructor is called.
    ~intensity_hdf5_ofile();
    
    //
    // Append a chunk of data, of length 'nt_chunk.
    // The 'intensity' and 'weight' arrays have shape (nfreq, npol, nt_chunk).
    // The mandatory 'chunk_ipos' arg is the index of the first sample in the chunk (in downsampled units).
    // The optional 'chunk_t0' arg is the arrival time of the first sample in the chunk (if omitted, will be inferred from 'chunk_ipos').
    //
    void append_chunk(ssize_t nt_chunk, float *intensity, float *weights, ssize_t chunk_ipos, double chunk_t0);
    void append_chunk(ssize_t nt_chunk, float *intensity, float *weights, ssize_t chunk_ipos);
};



// -------------------------------------------------------------------------------------------------
//
// intensity_network_stream: stream object which receives a packet stream from the network.
//
// (Writing a packet stream is handled by a separate class intensity_network_ostream, see below.)
//
// When the intensity_network_stream object is constructed, threads are automatically spawned to read and process
// packets.  The stream presents the incoming data to the "outside world" as a per-beam sequence 
// of regular arrays which are obtained by calling get_assembled_chunk().
//
// Reminder: normal shutdown sequence works as follows.
//  
//   - in network thread, end-of-stream packet is received (in intensity_network_stream::_network_thread_body())
//       - network thread calls intensity_network_stream::_network_thread_exit().
//       - stream state is advanced to 'stream_end_requested', this means that stream has exited but not all threads have joined.
//       - unassembled_ringbuf.end_stream() is called, which will tell the assembler thread that there are no more packets.
//
//   - in assembler thread, unassembled_ringbuf.get_packet_list() returns false.  
//       - the assembler thread loops over all beams ('assemblers') and calls assembled_chunk_ringbuf::end_stream()
//       - this sets assembled_chunk_ringbuf::doneflag
//
//   - in dedispersion thread, assembled_chunk_ringbuf::end_stream() returns an empty pointer.

struct packet_counts {
    struct timeval tv;
    double period; // in seconds
    std::unordered_map<uint64_t, uint64_t> counts;

    // Default constructor
    packet_counts();
    // Copy constructor
    packet_counts(const packet_counts& other);

    // Convert keys to IP:port format
    std::unordered_map<std::string, uint64_t> to_string() const;

    double start_time() const;
    
    void increment(const struct sockaddr_in& sender_addr, int nbytes);

    void update(const packet_counts& other);
};


class intensity_network_stream : noncopyable {
public:

    typedef std::function<void(std::vector<int> beam_id)> first_packet_listener;
    
    // The 'struct initializer' is used to construct the stream object.  A few notes:
    //
    //   - Beams are identified in the packet profile by an opaque uint16_t.  The 'beam_ids' argument
    //     should be the list of beam_ids which we expect to receive.
    //
    //   - The throw_exception_on_buffer_drop, throw_exception_on_assembler_miss flags are useful 
    //     for a unit test, but probably don't make sense otherwise.
    //
    //   - Our network protocol doesn't define any way of indicating end-of-stream.  This makes sense for the
    //     realtime search, but for testing we'd like to have a way of shutting down gracefully.  If the
    //     accept_end_of_stream_packets flag is set, then a special packet with nbeams=nupfreq=nt=0 is
    //     interpreted as an "end of stream" flag, and triggers shutdown of the network_stream.
    //
    //   - If 'nrfifreq' is initialized to a nonzero value, then memory will be allocated in each
    //     assembled_chunk for the RFI bitmask.  The "downstream" pipeline must fill the mask in each
    //     chunk.  If this does not happen quickly enough (more precisely, by the time the chunk
    //     leaves the top level of the telescoping ring buffer) then an exception will be thrown.

    struct initializer {
        int nbeams;

	std::shared_ptr<memory_slab_pool> memory_pool;
	std::vector<std::shared_ptr<output_device>> output_devices;

	int nupfreq = 0;
	int nrfifreq = 0;
	int nt_per_packet = 0;
	int fpga_counts_per_sample = 384;
	int stream_id = 0;   // only used in assembled_chunk::format_filename().

	// If 'nt_align' is set to a nonzero value, then the time sample index of the first
	// assembled_chunk in the stream must be a multiple of nt_align.  This is used in the
	// real-time server, to align all beams to the RFI and dedispersion block sizes.  Note
	// that if nt_align is enabled, then some packets may be dropped, and these will be
	// treated as assembler misses.
	int nt_align = 0;

	// If 'frame0_url' is a nonempty string, then assembler thread will retrieve frame0 info by "curling" the URL.
        std::string frame0_url = "";
        int frame0_timeout = 3000;

	// If ipaddr="0.0.0.0", then network thread will listen on all interfaces.
	std::string ipaddr = "0.0.0.0";
	int udp_port = constants::default_udp_port;

	bool force_reference_kernels = false;
	bool force_fast_kernels = false;
	bool emit_warning_on_buffer_drop = true;
	bool throw_exception_on_beam_id_mismatch = true;
	bool throw_exception_on_packet_mismatch = true;
	bool throw_exception_on_buffer_drop = false;
	bool throw_exception_on_assembler_miss = false;
	bool accept_end_of_stream_packets = true;
	bool deliberately_crash = false;  // deliberately crash dedispersion thread (for debugging purposes, obviously)

	// If nonempty, threads will be pinned to given list of cores.
	std::vector<int> network_thread_cores;
	std::vector<int> assembler_thread_cores;

	// The recv_socket_timeout determines how frequently the network thread wakes up, while blocked waiting
	// for packets.  The purpose of the periodic wakeup is to check whether intensity_network_stream::end_stream()
	// has been called, and check the timeout for flushing data to assembler threads.
	//
	// The stream_cancellation_latency_usec arg determines how frequently the network thread checks whether
	// intensity_network_stream::end_stream() has been called.  (We don't do this check in every iteration of
	// the packet read loop, since it requires acquiring a lock.)

	int max_packet_size = 9000;
	int socket_bufsize = constants::default_socket_bufsize;
	int socket_timeout_usec = 10000;                // 0.01 sec
	int stream_cancellation_latency_usec = 10000;   // 0.01 sec

        int packet_count_period_usec = 1000000; // 1 sec
        int max_packet_history_size = 3600; // keep an hour of history

	// The 'unassembled_ringbuf' is between the network thread and assembler thread.
	int unassembled_ringbuf_capacity = 16;
	int max_unassembled_packets_per_list = 16384;
	int max_unassembled_nbytes_per_list = 8 * 1024 * 1024;
	int unassembled_ringbuf_timeout_usec = 250000;   // 0.25 sec

	// The 'assembled_ringbuf' is between the assembler thread and processing threads.
	int assembled_ringbuf_capacity = 8;

	// The 'telescoping_ringbuf' stores assembled_chunks for retrieval by RPC.
	// Its capacity is a vector, whose length is the number of downsampling levels,
	// and whose elements are the number of assembled_chunks at each level.
        std::vector<int> telescoping_ringbuf_capacity;

	// A temporary hack that will go away soon.
	// Sleep for specified number of seconds, after intensity_stream starts up.
	double sleep_hack = 0.0;
    };

    // Event counts are kept in an array of the form int64_t[event_type::num_types].
    // Currently we don't do anything with the event counts besides print them at the end,
    // but they're intended to be a hook for real-time monitoring via RPC's.

    enum event_type {
	byte_received = 0,
	packet_received = 1,
	packet_good = 2,
	packet_bad = 3,
	packet_dropped = 4,            // network thread will drop packets if assembler thread runs slow and ring buffer overfills
	packet_end_of_stream = 5,
	beam_id_mismatch = 6,          // beam id in packet doesn't match any element of initializer::beam_ids
	stream_mismatch = 7,           // stream params (nupfreq, nt_per_packet, fpga_counts_per_sample) don't match values sepcified in constructor
	assembler_hit = 8,
	assembler_miss = 9,
	assembled_chunk_dropped = 10,  // assembler thread will drop assembled_chunks if processing thread runs slow
	assembled_chunk_queued = 11,
	num_types = 12                 // must be last
    };

    const initializer ini_params;

    // The largest FPGA count in a received packet.
    std::atomic<uint64_t> packet_max_fpga_seen;
    
    // It's convenient to initialize intensity_network_streams using a static factory function make(),
    // rather than having a public constructor.  Note that make() spawns network and assembler threads,
    // but won't listen for packets until start_stream() is called.  Between calling make() and start_stream(),
    // you'll want to spawn "consumer" threads which call get_assembled_chunk().

    static std::shared_ptr<intensity_network_stream> make(const initializer &ini_params);

    // High level control.
    void start_stream();         // tells network thread to start listening for packets (if stream has already started, this is not an error)
    void end_stream();           // requests stream exit (but stream will stop after a few timeouts, not immediately)
    void join_threads();         // should only be called once, does not request stream exit, blocks until network and assembler threads exit

    void reset_stream();

    void flush_end_of_stream();
    
    void wait_for_first_packet();

    std::vector<int> get_beam_ids();

    // Returns the first fpgacount of the first chunk sent downstream by
    // the given beam id.
    // Raises runtime_error if the first packet has not been received yet.
    uint64_t get_first_fpgacount();
    
    uint64_t get_frame0_nano();

    void add_first_packet_listener(first_packet_listener f);
    
    // This is the main routine called by the processing threads, to read data from one beam
    // corresponding to ini_params.beam_ids[assembler_ix].  (Note that the assembler_index
    // satisifes 0 <= assembler_ix < ini_params.beam_ids.size(), and is not a beam_id.)

    std::shared_ptr<assembled_chunk> get_assembled_chunk(int assembler_index, bool wait=true);
    
    // Can be called at any time, from any thread.  Note that the event counts returned by get_event_counts()
    // may slightly lag the real-time event counts (this behavior derives from wanting to avoid acquiring a
    // lock in every iteration of the packet read loop).
    std::vector<int64_t> get_event_counts();

    std::unordered_map<std::string, uint64_t> get_perhost_packets();

    std::vector<std::unordered_map<std::string, uint64_t> > get_statistics();

    // Retrieves chunks from one or more ring buffers.  The uint64_t
    // return value is a bitmask of l1_ringbuf_level values saying
    // where in the ringbuffer the chunk was found; this is an
    // implementation detail revealed for debugging purposes.
    //
    // If a vector of beam numbers is given, only the ring buffers for
    // those beams will be returned; otherwise the ring buffers for
    // all beams will be returned.
    std::vector< std::vector< std::pair<std::shared_ptr<assembled_chunk>, uint64_t> > >
    get_ringbuf_snapshots(const std::vector<int> &beams = std::vector<int>(),
                          uint64_t min_fpga_counts=0, uint64_t max_fpga_counts=0);

    // Searches the telescoping ring buffer for the given beam and fpgacounts start.
    // If 'toplevel' is true, then only the top level of the ring buffer is searched.
    // Returns an empty pointer iff stream has ended, and chunk is requested past end-of-stream.
    // If anything else goes wrong, an exception will be thrown.
    std::shared_ptr<assembled_chunk> find_assembled_chunk(int beam, uint64_t fpga_counts, bool toplevel=true);

    // Returns the last FPGA count processed by each of the assembler,
    // (in the same order as the "beam_ids" array), flushed downstream,
    // and retrieved by downstream callers.
    void get_max_fpga_count_seen(std::vector<uint64_t> &flushed,
                                 std::vector<uint64_t> &retrieved);
    
    // If period = 0, returns the packet rate with timestamp closest
    // to *start*.  If *start* is zero or negative, it is interpreted
    // as seconds relative to now; otherwise as gettimeofday()
    // seconds.  Zero grabs the most recent rate.  If *period* is
    // given, the rate samples overlapping [start, start+period] are
    // summed.
    std::shared_ptr<packet_counts> get_packet_rates(double start, double period);

    // Returns the set of packet rates with timestamps overlapping *start* to *end*.
    // If *start* is zero, treat as NOW.  If *start* or *end* are negative, NOW - that many seconds.
    std::vector<std::shared_ptr<packet_counts> > get_packet_rate_history(double start, double end, double period);
    
    // For debugging/testing purposes: pretend that the given
    // assembled_chunk has just arrived.  Returns true if there was
    // room in the ring buffer for the new chunk.
    bool inject_assembled_chunk(assembled_chunk* chunk);

    // For debugging/testing: pretend a packet has just arrived.
    void fake_packet_from(const struct sockaddr_in& sender, int nbytes);

    void start_forking_packets(int beam, int destbeam, const struct sockaddr_in& dest);
    void stop_forking_packets (int beam, int destbeam, const struct sockaddr_in& dest);
    void pause_forking_packets();
    void resume_forking_packets();

    // stream_to_files(): for streaming incoming data to disk.
    //
    //   'filename_pattern': see assembled_chunk::format_filename below (empty string means "streaming disabled")
    //
    //      On the DRAO backend, this will probably be one of the following two possibilities:
    //         /local/acq_data/(ACQNAME)/beam_(BEAM)/chunk_(CHUNK).msg                      (save to local SSD in node)
    //         /frb-archiver-(STREAM)/acq_data/(ACQNAME)/beam_(BEAM)/chunk_(CHUNK).msg      (save to NFS, with load-balancing)
    //
    //   'beam_ids': list of beam_ids to stream (if an empty list, this will also disable streaming)
    //
    //      This should be a subset of the beam_ids processed by the intensity_network_stream,
    //      which in turn in a subset of all beam_ids processed by the node.
    //    
    //   'priority': see write_chunk_request::priority below
    //
    // Throws an exception if anything goes wrong!  When called from an RPC thread, caller will want to
    // wrap in try..except, and use exception::what() to get the error message.

    void stream_to_files(const std::string &filename_pattern, const std::vector<int> &beam_ids, int priority, bool need_rfi);

    void get_streaming_status(std::string &filename_pattern,
                              std::vector<int> &beam_ids,
                              int &priority,
                              int &chunks_written,
                              size_t &bytes_written);

    // For debugging: print state.
    void print_state();

    ~intensity_network_stream();

protected:
    // Constant after construction, so not protected by lock
    std::vector<std::shared_ptr<assembled_chunk_ringbuf> > assemblers;

    std::vector<first_packet_listener> first_packet_listeners;

    std::vector<int> beam_ids;

    uint64_t first_fpgacount;
    
    std::map<int, std::shared_ptr<assembled_chunk_ringbuf> > beam_to_assembler;

    // Used to exchange data between the network and assembler threads
    std::unique_ptr<udp_packet_ringbuf> unassembled_ringbuf;

    // Written by network thread, read by outside thread
    // How much wall time do we spend waiting in recvfrom() vs processing?
    std::atomic<uint64_t> network_thread_waiting_usec;
    std::atomic<uint64_t> network_thread_working_usec;

    std::atomic<uint64_t> socket_queued_bytes;

    // I'm not sure how much it actually helps bottom-line performace, but it seemed like a good idea
    // to insert padding so that data accessed by different threads is in different cache lines.
    char _pad1[constants::cache_line_size];

    // Written by assembler thread, read by outside thread
    std::atomic<uint64_t> assembler_thread_waiting_usec;
    std::atomic<uint64_t> assembler_thread_working_usec;

    // Initialized to zero by constructor, set to nonzero value by assembler thread when first packet is received.
    std::atomic<uint64_t> frame0_nano;  // nanosecond time() value for fgpacount zero.

    char _pad1b[constants::cache_line_size];

    // Used only by the network thread (not protected by lock)
    //
    // Note on event counting implementation: on short timescales, the network and assembler 
    // threads accumulate event counts into "local" arrays 'network_thread_event_subcounts', 
    // 'assembler_thread_event_subcounts'.  On longer timescales, these local subcounts are 
    // accumulated into the global 'cumulative_event_counts'.  This two-level accumulation
    // scheme is designed to avoid excessive contention for the event_count_lock.

    int sockfd = -1;
    std::unique_ptr<udp_packet_list> incoming_packet_list;
    std::vector<int64_t> network_thread_event_subcounts;
    std::shared_ptr<packet_counts> network_thread_perhost_packets;
    char _pad2[constants::cache_line_size];

    // Used only by the assembler thread
    std::vector<int64_t> assembler_thread_event_subcounts;

    std::thread network_thread;
    std::thread assembler_thread;

    char _pad3[constants::cache_line_size];

    // State model.  These flags are protected by the state_lock and are set in sequence.
    // Note that the 'stream_ended_requested' flag means that the stream shutdown is imminent,
    // but doesn't mean that it has actually shut down yet, it may still be reading packets.
    // So far it hasn't been necessary to include a 'stream_ended' flag in the state model.

    std::mutex state_mutex;
    std::condition_variable cond_state_changed;
    
    bool stream_started = false;             // set asynchonously by calling start_stream()
    bool first_packet_received = false;
    bool stream_end_requested = false;       // can be set asynchronously by calling end_stream(), or by network/assembler threads on exit
    bool join_called = false;                // set by calling join_threads()
    bool threads_joined = false;             // set when both threads (network + assembler) are joined
    bool flush_end_of_stream_requested = false;
    bool stream_restart = false;

    char _pad4[constants::cache_line_size];

    std::mutex event_mutex;
    std::vector<int64_t> cumulative_event_counts;
    std::shared_ptr<packet_counts> perhost_packets;

    std::mutex packet_history_mutex;
    std::map<double, std::shared_ptr<packet_counts> > packet_history;
    
    // Streaming-related data (arguments to stream_to_files()).
    std::mutex stream_lock;
    std::string stream_filename_pattern;
    std::vector<int> stream_beam_ids;
    int stream_priority;
    bool stream_rfi_mask;
    int stream_chunks_written;
    size_t stream_bytes_written;

    struct packetfork {
        int beam;
        int destbeam;
        struct sockaddr_in dest;
    };

    std::mutex forking_mutex;
    std::vector<packetfork> forking_packets;
    int forking_socket = 0;
    std::atomic<bool> forking_paused;

    // The actual constructor is protected, so it can be a helper function 
    // for intensity_network_stream::make(), but can't be called otherwise.
    intensity_network_stream(const initializer &x);

    void _open_socket();
    void _network_flush_packets();
    void _add_event_counts(std::vector<int64_t> &event_subcounts);
    void _update_packet_rates(std::shared_ptr<packet_counts> last_packet_counts);

    std::shared_ptr<assembled_chunk_ringbuf> _assembler_for_beam(int beam_id);

    void network_thread_main();
    void assembler_thread_main();

    // Private methods called by the network thread.    
    void _network_thread_body();
    void _network_thread_one_stream();
    void _network_thread_exit();
    void _put_unassembled_packets();

    // Private methods called by the assembler thread.     
    void _assembler_thread_body();
    void _assembler_thread_exit();
    // initializes 'frame0_nano' by curling 'frame0_url', called when first packet is received.
    // NOTE that one must call curl_global_init() before, and curl_global_cleanup() after; in chime-frb-l1 we do this in the top-level main() method.
    void _fetch_frame0();
};


// -------------------------------------------------------------------------------------------------
//
// assembled_chunk
//
// This is the data structure which is returned by intensity_network_stream::get_assembled_chunk().
// It repesents a regular subarray of the intensity data in a single beam.  
//
// The data in the assembled_chunk is represented as 8-bit integers, with a coarser array of additive
// and mutiplicative offsets, but can be "decoded" to a simple floating-point array by calling
// assembled_chunk::decode().
//
// The 'offsets' and 'scales' arrays below are coarse-grained in both frequency and time, relative 
// to the 'data' array.  The level of frequency coarse-graining is the constructor argument 'nupfreq',
// and the level of time coarse-graining is the constructor-argument 'nt_per_packet'.  For no particular
// reason, the number of coarse-grained frequency channels is a compile-time constant (constants::nfreq_coarse)
// whereas the number of fine-grained time samples is also a compile-time constant (constants::nt_per_assembled_chunk).
//
// Summarizing:
//
//   fine-grained shape = (constants::nfreq_coarse * nupfreq, constants::nt_per_assembled_chunk)
//   coarse-grained shape = (constants::nfreq_coarse, constants::nt_per_assembled_chunk / nt_per_packet)
//
// A note on timestamps: there are several units of time used in different parts of the CHIME pipeline.
//
//   1 fpga count = 2.56e-6 seconds    (approx)
//   1 intensity sample = fpga_counts_per_sample * (1 fpga count)
//   1 assembled_chunk = constants::nt_per_assembled_chunk * (1 intensity sample)
//
// The 'isample' and 'ichunk' fields of the 'struct assembled_chunk' are simply defined by
//   isample = (initial fpga count) / fpga_counts_per_sample
//   ichunk = (initial fpga count) / (fpga_counts_per_sample * constants::nt_per_assembled_chunk)
//
// The 'binning' field of the 'struct assembled_chunk' is 1, 2, 4, 8, ... depending on what level
// of downsampling has been applied (i.e. the location of the assembled_chunk in the telescoping
// ring buffer).  Note that consecutive chunks src1, src2 with the same binning will satisfy
// (src2->ichunk == src1->ichunk + binning).


class assembled_chunk : noncopyable {
public:
    struct initializer {
	int beam_id = 0;
	int nupfreq = 0;
        int nrfifreq = 0;    // number of frequencies in downsampled RFI chain processing
	int nt_per_packet = 0;
	int fpga_counts_per_sample = 0;
	int binning = 1;
	int stream_id = 0;   // only used in assembled_chunk::format_filename().
	uint64_t ichunk = 0;
	bool force_reference = false;
	bool force_fast = false;

        // "ctime" in nanoseconds of FGPAcount zero
        uint64_t frame0_nano = 0;

	// If a memory slab has been preallocated from a pool, these pointers should be set.
	// Otherwise, both pointers should be empty, and the assembled_chunk constructor will allocate.
	std::shared_ptr<memory_slab_pool> pool;
        mutable memory_slab_t slab;
    };

    // Parameters specified at construction.
    const int beam_id = 0;
    const int nupfreq = 0;
    const int nrfifreq = 0;
    const int nt_per_packet = 0;
    const int fpga_counts_per_sample = 0;    // no binning factor applied here
    const int binning = 0;                   // either 1, 2, 4, 8... depending on level in telescoping ring buffer
    const int stream_id = 0;
    const uint64_t ichunk = 0;
    // "ctime" in nanoseconds of FGPAcount zero
    const uint64_t frame0_nano = 0;

    // Derived parameters.
    const int nt_coarse = 0;          // equal to (constants::nt_per_assembled_chunk / nt_per_packet)
    const int nscales = 0;            // equal to (constants::nfreq_coarse * nt_coarse)
    const int ndata = 0;              // equal to (constants::nfreq_coarse * nupfreq * constants::nt_per_assembled_chunk)
    const int nrfimaskbytes = 0;      // equal to (nrfifreq * constants::nt_per_assembled_chunk / 8)
    const uint64_t isample = 0;       // equal to ichunk * constants::nt_per_assembled_chunk
    const uint64_t fpga_begin = 0;    // equal to ichunk * constants::nt_per_assembled_chunk * fpga_counts_per_sample
    const uint64_t fpga_end = 0;      // equal to (ichunk+binning) * constants::nt_per_assembled_chunk * fpga_counts_per_sample

    // Note: you probably don't want to call the assembled_chunk constructor directly!
    // Instead use the static factory function assembed_chunk::make().
    assembled_chunk(const initializer &ini_params);
    virtual ~assembled_chunk();

    // Returns C time() (seconds since the epoch, 1970.0) of the first/last sample in this chunk.
    double time_begin() const;
    double time_end() const;

    // The following virtual member functions have default implementations,
    // which are overridden by the subclass 'fast_assembled_chunk' to be faster
    // on a CPU with the AVX2 instruction set, if certain conditions are met
    // (currently nt_per_packet==16 and nupfreq even).
    //
    // If 'prescale' is specified, then the intensity array written by decode()
    // will be multiplied by its value.  This is a temporary workaround for some
    // 16-bit overflow issues in bonsai.  (We currently don't need prescaling
    // in decode_subset(), but this could be added easily.)
    //
    // Warning (FIXME?): decode() and decode_subset() do not apply the RFI mask
    // in assembled_chunk::rfi_mask (if this exists).

    virtual void add_packet(const intensity_packet &p);
    virtual void decode(float *intensity, float *weights, int istride, int wstride, float prescale=1.0) const;
    virtual void decode_subset(float *intensity, float *weights, int t0, int nt, int istride, int wstride) const;
    virtual void downsample(const assembled_chunk *src1, const assembled_chunk *src2);   // downsamples data and RFI mask

    // Static factory functions which can return either an assembled_chunk or a fast_assembled_chunk.
    static std::unique_ptr<assembled_chunk> make(const initializer &ini_params);
    static std::shared_ptr<assembled_chunk> read_msgpack_file(const std::string& filename);

    // Note: the hdf5 file format has been phased out now..
    void write_hdf5_file(const std::string &filename);
    void write_msgpack_file(const std::string &filename, bool compress,
                            uint8_t* buffer=NULL);

    // How big can the bitshuffle-compressed data for a chunk of this size become?
    size_t max_compressed_size();

    // Performs a printf-like pattern replacement on *pattern* given the parameters of this assembled_chunk.
    // Replacements:
    //   (STREAM)  -> %01i stream_id
    //   (BEAM)    -> %04i beam_id
    //   (CHUNK)   -> %08i ichunk
    //   (NCHUNK)  -> %02i  size in chunks
    //   (BINNING) -> %02i  size in chunks
    //   (FPGA0)   -> %012i start FPGA-counts
    //   (FPGAN)   -> %08i  FPGA-counts size
    std::string format_filename(const std::string &pattern) const;

    // Utility functions currently used only for testing.
    void fill_with_copy(const std::shared_ptr<assembled_chunk> &x);
    void randomize(std::mt19937 &rng);   // also randomizes rfi_mask (if it exists)

    static ssize_t get_memory_slab_size(int nupfreq, int nt_per_packet, int nrfifreq);

    // I wanted to make the following fields protected, but msgpack doesn't like it...

    // Primary buffers.
    float *scales = nullptr;   // 2d array of shape (constants::nfreq_coarse, nt_coarse)
    float *offsets = nullptr;  // 2d array of shape (constants::nfreq_coarse, nt_coarse)
    uint8_t *data = nullptr;   // 2d array of shape (constants::nfreq_coarse * nupfreq, constants::nt_per_assembled_chunk)
    uint8_t *rfi_mask = nullptr;   // 2d array of downsampled masks, packed bitwise; (nrfifreq x constants::nt_per_assembled_chunk / 8 bits)

    // False on initialization.
    // If the RFI mask is being saved (nrfifreq > 0), it will be subsequently set to True by the processing thread.
    std::atomic<bool> has_rfi_mask;

    std::atomic<int> packets_received;
  // how many assembler misses occurred while this chunk was active_chunk0?
    std::atomic<int> packets_missed;

    // Temporary buffers used during downsampling.
    float *ds_w2 = nullptr;    // 1d array of length (nt_coarse/2)
    float *ds_data = nullptr;  // 2d array of shape (nupfreq, constants::nt_per_assembled_chunk/2)
    int *ds_mask = nullptr;    // 2d array of shape (nupfreq, constants::nt_per_assembled_chunk/2)

    // Used in the write path, to keep track of writes to disk.
    std::mutex filename_mutex;
    std::unordered_set<std::string> filename_set;
    std::unordered_map<std::string, std::string> filename_map;  // hash output_device_name -> filename

protected:
    // The array members above (scales, ..., ds_mask) are packed into a single contiguous memory slab.
    std::shared_ptr<memory_slab_pool> memory_pool;
    memory_slab_t memory_slab;

    void _check_downsample(const assembled_chunk *src1, const assembled_chunk *src2);

    // Note: destructors must call _deallocate()!  
    // Otherwise the memory_slab can't be returned to the pool.
    void _deallocate();
};


// If the CPU has the AVX2 instruction set, and for some choices of packet parameters (the precise 
// criterion is nt_per_packet == 16 and nupfreq % 2 == 0), we can speed up the assembled_chunk
// member functions using assembly language kernels.

class fast_assembled_chunk : public assembled_chunk
{
public:
    // Note: you probably don't want to call the assembled_chunk constructor directly!
    // Instead use the static factory function assembed_chunk::make().
    fast_assembled_chunk(const assembled_chunk::initializer &ini_params);
    virtual ~fast_assembled_chunk();

    // Override viruals with fast assembly language versions.
    virtual void add_packet(const intensity_packet &p) override;
    virtual void decode(float *intensity, float *weights, int istride, int wstride, float prescale=1.0) const override;
    virtual void downsample(const assembled_chunk *src1, const assembled_chunk *src2) override;

};


 
// -------------------------------------------------------------------------------------------------
//
// memory_slab_pool


class memory_slab_pool {
public:
    // The 'verbosity' parameter has the following meaning:
    //   0 = ninja-quiet
    //   1 = a little output during initialization
    //   2 = debug trace of all allocations/deallocations
    memory_slab_pool(ssize_t nbytes_per_slab, ssize_t nslabs, const std::vector<int> &allocation_cores, int verbosity=1);
    ~memory_slab_pool();

    // Returns a new slab from the pool.
    //
    // If the pool is empty, then either a null pointer is returned (wait=false),
    // or get_slab() blocks until a slab is available (wait=true).
    //
    // If zero=true, then the new slab is zeroed.

    memory_slab_t get_slab(bool zero=true, bool wait=false);
    
    // Puts a slab back in the pool.
    // Note: 'p' will be set to a null pointer after put_slab() returns.
    void put_slab(memory_slab_t &p);

    int count_slabs_available();

    const ssize_t nbytes_per_slab;
    const ssize_t nslabs;
    const int verbosity;

protected:
    std::mutex lock;
    std::condition_variable cv;

    std::vector<memory_slab_t> slabs;
    ssize_t curr_size = 0;
    ssize_t low_water_mark = 0;

    // Called by constructor, in separate thread.
    void allocate(const std::vector<int> &allocation_cores);
};


// -------------------------------------------------------------------------------------------------
//
// output_device and helper classes.
//
// FIXME(?): it would be natural to add member functions of 'class output_device' or 
// 'class output_device_pool' which return summary information, such as number of
// chunks queued for writing, total number of chunks written, etc.  (Same goes for
// the memory_slab_pool!)


// write_chunk_request: mini-struct consisting of a (chunk, filename, priority) triple.
//
// If the write_callback() virtual function is overridden, it will be called when the write
// request completes, either successfully or unsuccessfully.  Note that the callback is made
// from the i/o thread!

struct write_chunk_request {
    std::shared_ptr<assembled_chunk> chunk;
    std::string filename;
    int priority = 0;
    bool need_rfi_mask = false;

    // Called when the status of this chunk has changed --
    // due to an error, successful completion, or, eg, RFI mask added.
    virtual void status_changed(bool finished, bool success,
                                const std::string &state,
                                const std::string &error_message) { }
    virtual ~write_chunk_request() { }

    // This comparator class is used below, to make an STL priority queue.
    struct _less_than {
	bool operator()(const std::shared_ptr<write_chunk_request> &x, const std::shared_ptr<write_chunk_request> &y)
	{
	    // Compare priority first, then time.
	    if (x->priority < y->priority) return true;
	    if (x->priority > y->priority) return false;

	    // Reversed inequality is intentional here (earlier time = higher priority)
	    return x->chunk->fpga_begin > y->chunk->fpga_begin;
	}
    };
};


// output_device: corresponds to one "device" where assembled_chunks may be written.
//
// An output_device is identified by a device_name string, which is a directory such
// as '/ssd0' or /nfs1'.  When the output_device is created (with output_device::make(),
// an i/o thread is automatically spawned, which runs in the background.

class output_device : noncopyable {
public:
    struct initializer {
	std::string device_name;       // if empty string, this output_device will write all chunks
	std::string io_thread_name;    // if empty string, a default io_thread_name will be used
	std::vector<int> io_thread_allowed_cores;  // if empty vector, the io thread will be unpinned

	// Possible 'verbosity' values:
	//   0 = log nothing
	//   1 = log unsuccessful writes
	//   2 = log i/o thread startup/shutdown
	//   3 = verbose trace of all write_requests
	int verbosity = 1;
    };

    const initializer ini_params;

    // This factory function should be used to create a new output_device.
    // Returns a shared_ptr to a new output_device, and spawns a thread which
    // runs in the background, and also holds a shared_ptr.
    static std::shared_ptr<output_device> make(const initializer &ini_params);

    // Can be called by either the assembler thread, or an RPC thread.
    // Returns 'false' if request could not be queued (because end_stream() was called)
    bool enqueue_write_request(const std::shared_ptr<write_chunk_request> &req);

    // Counts the number of queued write request chunks
    int count_queued_write_requests();

    // If 'wait' is true, then end_stream() blocks until pending writes are complete.
    // If 'wait' is false, then end_stream() cancels all pending writes.
    void end_stream(bool wait);

    // Blocks until i/o thread exits.  Call end_stream() first!
    void join_thread();

    // Called (by RFI thread) to notify that the given chunk has had its
    // RFI mask filled in.
    void filled_rfi_mask(const std::shared_ptr<assembled_chunk> &chunk);

protected:
    std::thread output_thread;

    // State model.
    bool end_stream_called = false;
    bool join_thread_called = false;
    bool thread_joined = false;

    // The priority queue of write requests to be run by the output thread.
    std::priority_queue<std::shared_ptr<write_chunk_request>,
			std::vector<std::shared_ptr<write_chunk_request>>,
			write_chunk_request::_less_than> _write_reqs;

    // Write requests where need_rfi_mask is set and the rfi mask isn't
    // yet available.
    std::vector<std::shared_ptr<write_chunk_request> > _awaiting_rfi;
    
    // The state model and request queue are protected by this lock and condition variable.
    std::mutex _lock;
    std::condition_variable _cond;

    // Temporary buffer used for assembled_chunk serialization, accessed only by I/O thread
    std::unique_ptr<uint8_t[]> _buffer;

    // Constructor is protected -- use output_device::make() instead!
    output_device(const initializer &ini_params);

    // This is "main()" for the i/o thread.
    void io_thread_main();

    // Helper function called by output thread.
    // Blocks until a write request is available; returns highest-priority request.
    std::shared_ptr<write_chunk_request> pop_write_request();
};


// output_device_pool: this little helper class is a wrapper around a list of
// output_devices with different device_names.

class output_device_pool {
public:
    // Note: length-zero vector is OK!
    output_device_pool(const std::vector<std::shared_ptr<output_device>> &streams);

    // Sends the write request to the appropriate output_device (based on filename).
    // Returns 'false' if request could not be queued.
    bool enqueue_write_request(const std::shared_ptr<write_chunk_request> &req);
    
    // Loops over output_devices, and calls either end_stream() or join_thread().
    void end_streams(bool wait);
    void join_threads();

protected:
    std::vector<std::shared_ptr<output_device>> streams;
    std::vector<std::string> device_names;
};


// -------------------------------------------------------------------------------------------------
//
// intensity_network_ostream: this class is used to packetize intensity data and send
// it over the network.
//
// The stream object presents the "oustide world" with an object which is fed regular
// chunks of data via a member function send_chunk().  Under the hood, send_chunk() is
// encoding each chunk as multiple UDP packets, which are handed off to a network
// thread running in the background.


class intensity_network_ostream : noncopyable {
public:

    // The 'struct initializer' is used to construct the ostream object.  A few notes:
    //
    //   - It's convenient for the input to intensity_network_ostream to be a pair of
    //     floating-point arrays (intensity, weights).  Since the packet protocol allows
    //     a boolean mask but not floating-point weights, we simply mask all samples
    //     whose weight is below some cutoff.  This is the 'wt_cutoff' parameter.
    //
    //   - Each UDP packet is a 4D logical array of shape 
    //        (nbeams, nfreq_coarse_per_packet, nupfreq, nt_per_packet). 
    //     Each chunk passed to send_chunk() is a 4D logical array of shape
    //        (nbeams, coarse_freq_ids.size(), nupfreq, nt_per_chunk)
    //     and will generally correspond to multiple packets.
    //
    //   - If throttle=false, then UDP packets will be written as quickly as possible.
    //     If throttle=true, then the packet-throttling logic works as follows:
    //
    //        - If target_gbps > 0.0, then packets will be transmitted at the
    //          specified rate in Gbps.
    //
    //        - If target_gbps = 0, then the packet transmit rate will be inferred
    //          from the value of 'fpga_counts_per_sample', assuming 2.56 usec per
    //          fpga count.  (This is the default.)

    struct initializer {
	std::string dstname;
	std::vector<int> beam_ids;
	std::vector<int> coarse_freq_ids;   // will usually be [ 0, 1, ..., 1023 ], but reordering is allowed

	int nupfreq = 0;
	int nt_per_chunk = 0;
	int nfreq_coarse_per_packet = 0;
	int nt_per_packet = 0;
	int fpga_counts_per_sample = 0;
	float wt_cutoff = constants::default_wt_cutoff;
        int bind_port = 0; // 0: don't bind; send from randomly assigned port
        std::string bind_ip = "0.0.0.0";

	bool throttle = true;
	double target_gbps = 0.0;

	bool is_blocking = true;
	bool emit_warning_on_buffer_drop = true;
	bool throw_exception_on_buffer_drop = false;
	bool send_end_of_stream_packets = true;
	bool print_status_at_end = true;
    };

    const initializer ini_params;
    
    // Some of these fields are redundant with fields in the ini_params, but the redundancy is convenient.
    const int nbeams;
    const int nfreq_coarse_per_packet;
    const int nfreq_coarse_per_chunk;
    const int nupfreq;
    const int nt_per_packet;
    const int nt_per_chunk;
    const int nbytes_per_packet;
    const int npackets_per_chunk;
    const int nbytes_per_chunk;
    const int elts_per_chunk;
    const uint64_t fpga_counts_per_sample;
    const uint64_t fpga_counts_per_packet;
    const uint64_t fpga_counts_per_chunk;

    // Note: this->target_gpbs is not identical to ini_params.target_gpbs, since there is a case 
    // (ini_params.throttle=true and ini_params.target_gbps=0) where ini_params.target_gbps is zero, 
    // but this->target_gbps has a nonzero value inferred from 'fpga_counts_per_sample'.

    double target_gbps = 0.0;

    // It's convenient to initialize intensity_network_ostreams using a static factory function make(),
    // rather than having a public constructor.  Note that make() spawns a network thread which runs
    // in the background.

    static std::shared_ptr<intensity_network_ostream> make(const initializer &ini_params);

    // Data is added to the network stream by calling send_chunk().  This routine packetizes
    // the data and puts the packets in a thread-safe ring buffer, for the network thread to send.
    //
    // The 'intensity' and 'weights' arrays have logical shape
    //   (nbeams, nfreq_coarse_per_chunk, nupfreq, nt_per_chunk).
    //
    // The 'istride/wstride' args are tmemory offsets between time series in the intensity/weights arrays.
    // To write this out explicitly, the intensity sample with logical indices (b,fcoarse,u,t) has memory location
    //   intensity + ((b + fcoarse*nupfreq) * nupfreq + u) * stride + t

    void send_chunk(const float *intensity, int istride, const float *weights, int wstride, uint64_t fpga_count);

    void get_statistics(int64_t& curr_timestamp,
                        int64_t& npackets_sent,
                        int64_t& nbytes_sent);
    
    // Called when there is no more data; sends end-of-stream packets.
    void end_stream(bool join_network_thread);

    ~intensity_network_ostream();

    // This is a helper function called by send_chunk(), but we make it public so that the unit tests can call it.
    void _encode_chunk(const float *intensity, int istride, const float *weights, int wstride, uint64_t fpga_count, const std::unique_ptr<udp_packet_list> &out);

    void print_status(std::ostream &os = std::cout);

    bool is_sending();

protected:
    std::vector<uint16_t> beam_ids_16bit;
    std::vector<uint16_t> coarse_freq_ids_16bit;

    int sockfd = -1;
    std::string hostname;
    uint16_t udp_port = constants::default_udp_port;
    
    std::mutex statistics_lock;
    int64_t curr_timestamp = 0;    // microseconds between first packet and most recent packet
    int64_t npackets_sent = 0;
    int64_t nbytes_sent = 0;

    // State model.
    std::thread network_thread;
    std::mutex state_lock;
    std::condition_variable cond_state_changed;
    bool network_thread_started = false;
    bool network_thread_joined = false;

    // Ring buffer used to exchange packets with network thread
    std::unique_ptr<udp_packet_ringbuf> ringbuf;
    
    // Temp buffer for packet encoding
    std::unique_ptr<udp_packet_list> tmp_packet_list;

    // The actual constructor is protected, so it can be a helper function 
    // for intensity_network_ostream::make(), but can't be called otherwise.
    intensity_network_ostream(const initializer &ini_params);

    void _network_thread_main();
    
    void _network_thread_body();

    // For testing purposes (eg, can create a subclass that randomly drops packets), a wrapper on the underlying packet send() function.
    virtual ssize_t _send(int socket, const uint8_t* packet, int nbytes, int flags);

    void _open_socket();
    void _send_end_of_stream_packets();
};


// -------------------------------------------------------------------------------------------------
//
// Miscellaneous
    

// Used in RPC's which return chunks, to indicate where the chunk was found.
// Note: implementation of assembled_chunk_ringbuf::get_ringbuf_snapshot() assumes
// that L1RB_LEVELn == 2^n, so be careful when modifying this!
enum l1_ringbuf_level {
    L1RB_DOWNSTREAM = 1,
    L1RB_LEVEL1 = 2,
    L1RB_LEVEL2 = 4,
    L1RB_LEVEL3 = 8,
    L1RB_LEVEL4 = 0x10,
    // queued for writing in the L1 RPC system
    L1RB_WRITEQUEUE = 0x100,
};


extern void pin_thread_to_cores(const std::vector<int> &core_list);


// Utility routine: converts a string to type T (only a few T's are defined; see lexical_cast.cpp)
// Returns true on success, false on failure
template<typename T> extern bool lexical_cast(const std::string &x, T &ret);

// Also defined in lexical_cast.cpp (for the same values of T)
template<typename T> extern const char *typestr();

// Version of lexical_cast() which throws exception on failure.
template<typename T> inline T lexical_cast(const std::string &x, const char *name="string")
{
    T ret;
    if (lexical_cast(x, ret))
	return ret;
    throw std::runtime_error("couldn't convert " + std::string(name) + "='" + x + "' to " + typestr<T>());
}

// Unit tests
extern void test_lexical_cast();
extern void test_packet_offsets(std::mt19937 &rng);
extern void test_avx2_kernels(std::mt19937 &rng);
extern void peek_at_unpack_kernel();


}  // namespace ch_frb_io

#endif // _CH_FRB_IO_HPP
