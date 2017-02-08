#ifndef _CH_FRB_IO_HPP
#define _CH_FRB_IO_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++11 support (g++ -std=c++11)"
#endif

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <atomic>
#include <random>
#include <thread>
#include <hdf5.h>

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
struct assembled_chunk;

// Defined in ch_frb_io_internals.hpp
struct intensity_packet;
struct udp_packet_list;
struct udp_packet_ringbuf;
class assembled_chunk_ringbuf;


// -------------------------------------------------------------------------------------------------
//
// Compile-time constants
// Note: many of these don't really need to be fixed at compile time, and could be runtime parameters instead.


namespace constants {
    static constexpr int cache_line_size = 64;

    // Number of "coarse" (i.e. pre-upchannelized) frequency channels.
    static constexpr int nfreq_coarse_tot = 1024;

    // For an explanation of this parameter, see class intensity_network_ostream below 
    static constexpr double default_wt_cutoff = 0.3;

    // Network parameters.
    //
    // The recv_socket_timeout determines how frequently the network thread wakes up, while blocked waiting
    // for packets.  The purpose of the periodic wakeup is to check whether intensity_network_stream::end_stream()
    // has been called, and check the timeout for flushing data to assembler threads.
    //
    // The max packet size params have been chosen around ~9KB, as appropriate for 1 Gbps ethernet
    // with jumbo frames and no fragmentation.
    //
    // The stream_cancellation_latency_usec arg determines how frequently the network thread checks whether
    // intensity_network_stream::end_stream() has been called.  (We don't do this check in every iteration of
    // the packet read loop, since it requires acquiring a lock.)

    static constexpr int default_udp_port = 10252;
    static constexpr int max_input_udp_packet_size = 9000;   // largest value the input stream will accept
    static constexpr int max_output_udp_packet_size = 8910;  // largest value the output stream will produce
    static constexpr int recv_socket_timeout_usec = 10000;   // 0.01 sec
    static constexpr int stream_cancellation_latency_usec = 10000;    // 0.01 sec
    static constexpr double default_gbps = 1.0;              // default transmit rate for output stream

#ifdef __APPLE__
    // osx seems to have very small limits on socket buffer size
    static constexpr int recv_socket_bufsize = 4 * 1024 * 1024;
    static constexpr int send_socket_bufsize = 4 * 1024 * 1024;
#else
    static constexpr int recv_socket_bufsize = 128 * 1024 * 1024;
    static constexpr int send_socket_bufsize = 128 * 1024 * 1024;
#endif

    // This applies to the ring buffer between the network _output_ thread, and callers of
    // intensity_network_ostream::send_chunk().
    static constexpr int output_ringbuf_capacity = 8;

    // Ring buffer between network _input_ thread and assembler thread.
    static constexpr int unassembled_ringbuf_capacity = 16;
    static constexpr int max_unassembled_packets_per_list = 16384;
    static constexpr int max_unassembled_nbytes_per_list = 8 * 1024 * 1024;
    static constexpr int unassembled_ringbuf_timeout_usec = 250000;   // 0.25 sec

    // Ring buffers between assembler thread and processing threads.
    static constexpr int assembled_ringbuf_capacity = 8;
    static constexpr int assembled_ringbuf_nlevels = 4;

    static constexpr int nt_per_assembled_chunk = 1024;

    // These parameters don't really affect anything but appear in range-checking asserts.
    static constexpr int max_allowed_beam_id = 65535;
    static constexpr int max_allowed_nupfreq = 64;
    static constexpr int max_allowed_nt_per_packet = 1024;
    static constexpr int max_allowed_fpga_counts_per_sample = 3200;
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
    // The 'out_stride' arg can be negative, if reversing the channel ordering is desired.
    // If out_stride is zero, it defaults to out_nt.
    //
    void get_unpolarized_intensity(float *out_int, float *out_wt, int out_t0, int out_nt, int out_stride=0) const;

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
    //   ipos0 = index of first sample in file (in downsampled units, i.e. one sample is ~1.3 msec, not ~2.5 usec)
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


class intensity_network_stream : noncopyable {
public:
    
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

    struct initializer {
	std::vector<int> beam_ids;
	int udp_port = constants::default_udp_port;
	bool mandate_reference_kernels = false;
	bool mandate_fast_kernels = false;
	bool emit_warning_on_buffer_drop = true;
	bool throw_exception_on_beam_id_mismatch = true;
	bool throw_exception_on_first_packet_mismatch = true;
	bool throw_exception_on_buffer_drop = false;
	bool throw_exception_on_assembler_miss = false;
	bool accept_end_of_stream_packets = true;
        std::vector<int> ringbuf_n;
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
	first_packet_mismatch = 7,     // stream params (nupfreq, nt_per_packet, fpga_counts_per_sample) don't match first packet received
	assembler_hit = 8,
	assembler_miss = 9,
	assembled_chunk_dropped = 10,  // assembler thread will drop assembled_chunks if processing thread runs slow
	assembled_chunk_queued = 11,
	num_types = 12                 // must be last
    };
    
    // It's convenient to initialize intensity_network_streams using a static factory function make(),
    // rather than having a public constructor.  Note that make() spawns network and assembler threads,
    // but won't listen for packets until start_stream() is called.  Between calling make() and start_stream(),
    // you'll want to spawn "consumer" threads which call get_assembled_chunk().

    static std::shared_ptr<intensity_network_stream> make(const initializer &ini_params);

    // High level control.
    void start_stream();         // tells network thread to start listening for packets
    void end_stream();           // requests stream exit (but stream will stop after a few timeouts, not immediately)
    void join_threads();         // should only be called once, does not request stream exit, blocks until network and assembler threads exit

    // This is the main routine called by the processing threads, to read data from one beam
    // corresponding to ini_params.beam_ids[assembler_ix].  (Note that the assembler_index
    // satisifes 0 <= assembler_ix < ini_params.beam_ids.size(), and is not a beam_id.)

    std::shared_ptr<assembled_chunk> get_assembled_chunk(int assembler_index,
                                                         bool wait=true);

    // Will block until first packet is received, or stream ends.  Returns true if packet was received.
    bool get_first_packet_params(int &nupfreq, int &nt_per_packet, uint64_t &fpga_counts_per_sample, uint64_t &fpga_count);
    
    // Can be called at any time, from any thread.  Note that the event counts returned by get_event_counts()
    // may slightly lag the real-time event counts (this behavior derives from wanting to avoid acquiring a
    // lock in every iteration of the packet read loop).

    initializer get_initializer();
    std::vector<int64_t> get_event_counts();

    std::unordered_map<std::string, uint64_t> get_perhost_packets();
    std::unordered_map<uint64_t, uint64_t> get_perhost_packets_raw();

    std::vector<std::unordered_map<std::string, uint64_t> > get_statistics();

    // Retrieves chunks from one or more ring buffers.  The uint64_t
    // return value is a bitmask of l1_ringbuf_level values saying
    // where in the ringbuffer the chunk was found; this is an
    // implementation detail revealed for debugging purposes
    std::vector< std::vector< std::pair<std::shared_ptr<assembled_chunk>, uint64_t> > >
    get_ringbuf_snapshots(std::vector<uint64_t> &beams,
                          uint64_t min_fpga_counts=0, uint64_t max_fpga_counts=0);

    // For debugging/testing purposes: pretend that the given
    // assembled_chunk has just arrived.  Returns true if there was
    // room in the ring buffer for the new chunk.
    bool inject_assembled_chunk(assembled_chunk* chunk);

    // For debugging/testing: stream data to disk.  Filename pattern: see assembled_chunk::format_filename.  Empty string to turn off streaming.
    void stream_to_files(const std::string& filename_pattern);

    // For debugging: print state.
    void print_state();

    ~intensity_network_stream();

protected:
    // Constant after construction, so not protected by lock
    const initializer ini_params;
    const int nassemblers = 0;

    // This is initialized by the assembler thread before it sets 'first_packet_received' flag.
    // Therefore, other threads can access it without a lock, but should wait for this flag to be set (which does
    // require a lock).  There is a corner case where the vector is still length-zero after the flag gets set.
    // This happens if the stream was asynchronously cancelled before receiving the first packet.

    std::vector<std::unique_ptr<assembled_chunk_ringbuf> > assemblers;
    
    // These fields are initialized from the first packet received ("fp_" stands for "first packet").
    // They are initialized by the assembler thread, which then advances the state model to "first_packet_received".
    // They are not protected by a lock!  This is OK as long as non-assembler threads access them read-only, and
    // only after checking the first_packet_received flag (which does require a lock).

    uint16_t fp_nupfreq = 0;
    uint16_t fp_nt_per_packet = 0;
    uint16_t fp_fpga_counts_per_sample = 0;
    uint64_t fp_fpga_count = 0;

    // Used to exchange data between the network and assembler threads
    std::unique_ptr<udp_packet_ringbuf> unassembled_ringbuf;

    // Written by network thread, read by outside thread
    // How much wall time do we spend waiting in recvfrom() vs processing?
    std::atomic_uint64_t network_thread_recv_usec;
    std::atomic_uint64_t network_thread_processing_usec;

    // I'm not sure how much it actually helps bottom-line performace, but it seemed like a good idea
    // to insert padding so that data accessed by different threads is in different cache lines.
    char _pad1[constants::cache_line_size];

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
    std::unordered_map<uint64_t, uint64_t> network_thread_perhost_packets;
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

    pthread_mutex_t state_lock;
    pthread_cond_t cond_state_changed;       // threads wait here for state to change

    bool stream_started = false;             // set asynchonously by calling start_stream()
    bool first_packet_received = false;      // set by network thread
    bool assemblers_initialized = false;     // set by assembler thread
    bool stream_end_requested = false;       // can be set asynchronously by calling end_stream(), or by network/assembler threads on exit
    bool join_called = false;                // set by calling join_threads()
    char _pad4[constants::cache_line_size];

    pthread_mutex_t event_lock;
    std::vector<int64_t> cumulative_event_counts;
    std::unordered_map<uint64_t, uint64_t> perhost_packets;

    std::string stream_filename;

    // The actual constructor is protected, so it can be a helper function 
    // for intensity_network_stream::make(), but can't be called otherwise.
    intensity_network_stream(const initializer &x);

    void _wait_for_assemblers_initialized(bool prelocked=false);

    void _open_socket();
    void _network_flush_packets();
    void _add_event_counts(std::vector<int64_t> &event_subcounts);

    void network_thread_main();
    void assembler_thread_main();

    // Private methods called by the network thread.    
    void _network_thread_body();
    void _network_thread_exit();
    void _put_unassembled_packets();

    // Private methods called by the assembler thread.     
    void _assembler_thread_body();
    void _assembler_thread_exit();
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


struct assembled_chunk : noncopyable {
    // Stream parameters, specified at construction.
    const int beam_id = 0;
    const int nupfreq = 0;
    const int nt_per_packet = 0;
    int fpga_counts_per_sample = 0;

    // More parameters which are constant after construction.
    const int nt_coarse = 0;   // equal to (constants::nt_per_assembled_chunk / nt_per_packet)
    const int nscales = 0;     // equal to (constants::nfreq_coarse * nt_coarse)
    const int ndata = 0;       // equal to (constants::nfreq_coarse * nupfreq * constants::nt_per_assembled_chunk)

    // Chunks are indexed by 'ichunk', which differs by 1 in adjacent chunks.
    //
    // Reminder: there are several units of time used in different parts of the CHIME backend!
    //   1 assembled_chunk = constants::nt_per_assembled_chunk * (1 intensity sample)
    //   1 intensity sample = assembled_chunk::fpga_counts_per_sample * (1 fpga count)
    //   1 fpga count = 2.56e-6 seconds    (approx)

    uint64_t ichunk = 0;
    uint64_t isample = 0;   // always equal to ichunk * constants::nt_per_assembled_chunk

    // How many original samples have been binned into each sample in this chunk.  Used in the L1 telescoping ring buffer, where binning = 2, 4, 8.
    int binning = 1;

    float *scales = nullptr;   // shape (constants::nfreq_coarse, nt_coarse)
    float *offsets = nullptr;  // shape (constants::nfreq_coarse, nt_coarse)
    uint8_t *data = nullptr;   // shape (constants::nfreq_coarse, nupfreq, constants::nt_per_assembled_chunk)

    bool msgpack_bitshuffle = false;

    assembled_chunk(int beam_id, int nupfreq, int nt_per_packet, int fpga_counts_per_sample, uint64_t ichunk);
    virtual ~assembled_chunk();

    // Performs a printf-like pattern replacement on *pattern* given the parameters of this assembled_chunk.
    // Replacements:
    //   (BEAM)    -> %04i beam_id
    //   (CHUNK)   -> %08i ichunk
    //   (NCHUNK)  -> %02i  size in chunks
    //   (BINNING) -> %02i  size in chunks
    //   (FPGA0)   -> %012i start FPGA-counts
    //   (FPGAN)   -> %08i  FPGA-counts size
    std::string format_filename(const std::string &pattern) const;

    // the first fpga-counts sample in this chunk
    uint64_t fpgacounts_begin() const { return isample * fpga_counts_per_sample / binning; }

    // the number of fpga-counts in this chunk
    uint64_t fpgacounts_N() const { return (constants::nt_per_assembled_chunk * fpga_counts_per_sample); }

    // the last fpga-counts sample in this chunk + 1
    uint64_t fpgacounts_end() const { return fpgacounts_begin() + fpgacounts_N(); }

    // These are virtual so that subclasses can be written with optimized implementations 
    // for specific parameter choices (e.g. full CHIME nt_per_packet=16)
    virtual void add_packet(const intensity_packet &p);
    virtual void decode(float *intensity, float *weights, int stride) const;

    virtual void decode_subset(float *intensity, float *weights,
                               int t0, int nt, int stride) const;

    // Downsamples (in time) the two given source chunks, writing into the
    // given destination chunk.
    // *dest* can be equal to *src1* to downsample in-place.
    // If *dest* is NULL, a new assembled_chunk will be allocated.
    static assembled_chunk* downsample(assembled_chunk* dest,
                                       const assembled_chunk* src1,
                                       const assembled_chunk* src2);

    // Static factory function which returns either the assembled_chunk base class, or the fast_assembled_chunk
    // subclass (see below), based on the packet parameters.
    static std::unique_ptr<assembled_chunk> make(int beam_id, int nupfreq, int nt_per_packet, int fpga_counts_per_sample, uint64_t ichunk, bool force_reference=false, bool force_fast=false);

    static std::shared_ptr<assembled_chunk> read_msgpack_file(const std::string& filename);

    // Utility functions currently used only for testing.
    void fill_with_copy(const std::shared_ptr<assembled_chunk> &x);
    void randomize(std::mt19937 &rng);

    // HDF5 file output
    void write_hdf5_file(const std::string &filename);

    // msgpack file output
    void write_msgpack_file(const std::string &filename);

};

// For some choices of packet parameters (the precise criterion is nt_per_packet == 16 and
// nupfreq % 2 == 0) we can speed up assembled_chunk::add_packet() and assembled_chunk::decode()
// using assembly language kernels.

struct fast_assembled_chunk : public assembled_chunk
{
    // Constructor throws exception unless nt_per_packet == 16 and (nupfreq % 2) == 0.
    fast_assembled_chunk(int beam_id, int nupfreq, int nt_per_packet, int fpga_counts_per_sample, uint64_t ichunk);

    // Overrides assembled_chunk::add_packet() and assembled_chunk::decode() with fast assembly language versions.
    virtual void add_packet(const intensity_packet &p) override;
    virtual void decode(float *intensity, float *weights, int stride) const override;
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
    //   - The target_gbps arg enables throttling of the output to some target bandwidth in Gbps.

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
	double target_gbps = constants::default_gbps;   // if 0.0, then data will be written as quickly as possible!
        int bind_port = 0; // 0: don't bind; send from randomly assigned port

	bool is_blocking = true;
	bool emit_warning_on_buffer_drop = true;
	bool throw_exception_on_buffer_drop = false;
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
    const double target_gbps;

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
    // The 'stride' arg is the memory offset between time series in the arrays.
    // To write this out explicitly, the intensity sample with logical indices (b,fcoarse,u,t) has memory location
    //   intensity + ((b + fcoarse*nupfreq) * nupfreq + u) * stride + t

    void send_chunk(const float *intensity, const float *weights, int stride, uint64_t fpga_count);

    void get_statistics(int64_t& curr_timestamp,
                        int64_t& npackets_sent,
                        int64_t& nbytes_sent);
    
    // Called when there is no more data; sends end-of-stream packets.
    void end_stream(bool join_network_thread);

    ~intensity_network_ostream();

    // This is a helper function called by send_chunk(), but we make it public so that the unit tests can call it.
    void _encode_chunk(const float *intensity, const float *weights, int stride, uint64_t fpga_count, const std::unique_ptr<udp_packet_list> &out);
    
protected:
    std::vector<uint16_t> beam_ids_16bit;
    std::vector<uint16_t> coarse_freq_ids_16bit;

    int sockfd = -1;
    std::string hostname;
    uint16_t udp_port = constants::default_udp_port;
    
    pthread_mutex_t statistics_lock;
    int64_t curr_timestamp = 0;    // microseconds between first packet and most recent packet
    int64_t npackets_sent = 0;
    int64_t nbytes_sent = 0;

    // State model.
    pthread_t network_thread;
    pthread_mutex_t state_lock;
    pthread_cond_t cond_state_changed;
    bool network_thread_started = false;
    bool network_thread_joined = false;

    // Ring buffer used to exchange packets with network thread
    std::unique_ptr<udp_packet_ringbuf> ringbuf;
    
    // Temp buffer for packet encoding
    std::unique_ptr<udp_packet_list> tmp_packet_list;

    // The actual constructor is protected, so it can be a helper function 
    // for intensity_network_ostream::make(), but can't be called otherwise.
    intensity_network_ostream(const initializer &ini_params);

    static void *network_pthread_main(void *opaque_args);

    void _network_thread_body();

    // For testing purposes (eg, can create a subclass that randomly drops packets), a wrapper on the underlying packet send() function.
    virtual ssize_t _send(int socket, const uint8_t* packet, int nbytes, int flags);

    void _open_socket();
    void _announce_end_of_stream();
};


// -------------------------------------------------------------------------------------------------
//
// Miscellaneous


// Utility routine: converts a string to type T (only a few T's are defined; see lexical_cast.cpp)
// Returns true on success, false on failure
template<typename T> extern bool lexical_cast(const std::string &x, T &ret);

// Also defined in lexical_cast.cpp (for the same values of T)
template<typename T> extern const char *typestr();

// Version of lexical_cast() which throws exception on failure.
template<typename T> inline T lexical_cast(const std::string &x)
{
    T ret;
    if (lexical_cast(x, ret))
	return ret;
    throw std::runtime_error("couldn't convert string '" + x + "' to " + typestr<T>());
}

// Unit tests
extern void test_lexical_cast();
extern void test_packet_offsets(std::mt19937 &rng);
extern void test_avx2_kernels(std::mt19937 &rng);
extern void peek_at_unpack_kernel();


}  // namespace ch_frb_io

#endif // _CH_FRB_IO_HPP
