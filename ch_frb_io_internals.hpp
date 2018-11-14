#ifndef _CH_FRB_IO_INTERNALS_HPP
#define _CH_FRB_IO_INTERNALS_HPP

#if (__cplusplus < 201103) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
#error "This source file needs to be compiled with C++11 support (g++ -std=c++11)"
#endif

#include <cmath>
#include <cstring>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <memory>

#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>

#include "ch_frb_io.hpp"


// This will compile to a "hint" for CPU branch prediction.  We use it mainly for error detection 
// in critical paths (e.g. packet parsing, where there are many possible ways a packet can be bad,
// leading to many unlikely branches in the code).  I found that it gives a few-percent speedup
// if used consistently, so it's not important, but seems worth doing anyway since it's so easy.

#ifndef _unlikely
#define _unlikely(cond)  (__builtin_expect(cond,0))
#endif


// ch_assert(): assert-like macro, but throws an exception instead of terminating the process.

#ifndef ch_assert
#define ch_assert(cond) ch_assert2(cond, __LINE__)
#define ch_assert2(cond,line) \
    do { \
        if (_unlikely(!(cond))) { \
	    const char *msg = "ch_frb_io: assertion '" __STRING(cond) "' failed (" __FILE__ ":" __STRING(line) ")\n"; \
	    throw std::runtime_error(msg); \
	} \
    } while (0)
#endif


namespace ch_frb_io {
#if 0
}; // pacify emacs c-mode
#endif


// -------------------------------------------------------------------------------------------------
//
// struct intensity_packet: a lightweight struct representing one UDP packet.
// See L0_L1_packet.hpp for a more verbose description of the packet format.


struct intensity_packet {
    // "Header fields".   These 24 bytes should have the same ordering and byte count as the 
    // "on-wire" packet, since we use memcpy(..., 24) to initialize them from the raw packet data.
    uint32_t  protocol_version;
    int16_t   data_nbytes;
    uint16_t  fpga_counts_per_sample;
    uint64_t  fpga_count;
    uint16_t  nbeams;
    uint16_t  nfreq_coarse;
    uint16_t  nupfreq;
    uint16_t  ntsamp;

    // "Pointer" fields
    uint16_t  *beam_ids;          // 1D array of length nbeams
    uint16_t  *coarse_freq_ids;   // 1D array of length nfreq_coarse
    float     *scales;            // 2D array of shape (nbeam, nfreq_coarse)
    float     *offsets;           // 2D array of shape (nbeam, nfreq_coarse)
    uint8_t   *data;              // array of shape (nbeam, nfreq_coarse, nupfreq, ntsamp)


    static inline int header_size(int nbeams, int nfreq_coarse)
    {
	return 24 + 2*nbeams + 2*nfreq_coarse + 8*nbeams*nfreq_coarse;
    }

    static inline int packet_size(int nbeams, int nfreq_coarse, int nupfreq, int nt_per_packet)
    {
	int data_size = nbeams * nfreq_coarse * nupfreq * nt_per_packet;
	return header_size(nbeams, nfreq_coarse) + data_size;
    }


    // Initializes a 'struct intensity_packet' from raw packet data.  The "pointer" fields of the
    // struct intensity_packet are initialized to pointers into the 'src' buffer, so the caller is
    // responsible for ensuring that this buffer doesn't get freed while the struct intensity_packet 
    // is in scope.
    //
    // Does a bunch of sanity checks and returns 'true' if packet is good, 'false' if bad.
    // (See extended comment in intensity_packet.cpp for a complete list of checks performed.)

    bool decode(const uint8_t *src, int src_nbytes);

    
    // Encodes a floating-point array of intensities into raw packet data, before sending packet.
    // The semantics of encode() aren't very intuitive, so we document them carefully here!
    //
    //    - Caller should initialize the "header" fields of the struct intensity packet.
    //
    //    - Caller should initialize the pointer fields 'beam_ids' and 'coarse_freq_ids' to 
    //      point to arrays of appropriate size.
    //
    //    - Caller should initialize the 'intensity' and 'weights' arrays to point to logical arrays 
    //      of shape (nbeams, nfreq_coarse, nupfreq, ntsamp).  This is the data that will be encoded 
    //      into the packet.  The istride/wstride arguments are defined so that the intensity array element with 
    //      logical indices (b,f,u,t) has memory location
    //
    //          intensity + b*beam_stride + (f*nupfreq+u)*freq_stride + t
    //
    //      and likewise for the weights.
    //
    //    - Since the binary packet format doesn't support weights but does support a boolean mask,
    //      encode() simply masksdata whose weight is below the 'wt_cutoff' argument.
    //
    //    - Caller must ensure that 'dst' points to a large enough buffer to encode the packet.
    //      For example, in intensity_network_ostream.cpp, we compute the buffer size ahead of
    //      time using intensity_packet::packet_size().
    //
    //    - Caller doesn't need to initialize the pointer fields 'scales', 'offsets', 'data'.
    //      These pointers are initialized in encode(), to point into the appropriate locations
    //      in the 'dst' buffer.  The actual scales and offsets are computed in encode() based
    //      on the mean and variance of the data.  The scales are chosen so that the intensities
    //      are masked if they deviate from the mean by approx 5 sigma.
    //
    // Returns size of the encoded packet in bytes.
    //
    // Caveat emptor: encode() doesn't do any argument checking at all, so it's easy to segfault
    // by calling it wrong!

    int encode(uint8_t *dst, 
	       const float *intensity, int beam_istride, int freq_istride, 
	       const float *weights, int beam_wstride, int freq_wstride, 
	       float wt_cutoff);


    // Currently used only for debugging
    int find_coarse_freq_id(int id) const;
    bool contains_coarse_freq_id(int id) const;
};


// -------------------------------------------------------------------------------------------------
//
// udp_packet_list: a buffer containing opaque UDP packets.
// udp_packet_ringbuf: a thread-safe ring buffer for exchanging udp_packet_lists between threads.


struct udp_packet_list {
    // Capacity of buffer
    const int max_npackets;
    const int max_nbytes;  // summed over all packets

    // Current buffer size
    int curr_npackets = 0;
    int curr_nbytes = 0;   // summed over all packets
    bool is_full = false;

    // Packets are concatenated into the 'buf' array, and off_buf[i] stores the offset
    // of the i-th packet relative to the start of the buffer.  It's convenient to set
    // the sentinel value
    //   off_buf[curr_npackets] = curr_nbytes
    // so that the i-th packet always has size (off_buf[i+1] - off_buf[i])

    std::unique_ptr<uint8_t[]> buf;   // points to an array of length (max_nbytes + max_packet_size).
    std::unique_ptr<int[]> off_buf;   // points to an array of length (max_npackets + 1).

    // Bare pointers.
    uint8_t *data_start = nullptr;    // points to &buf[0]
    uint8_t *data_end = nullptr;      // points to &buf[curr_nbytes]
    int *packet_offsets = nullptr;    // points to &off_buf[0].  Note that packet_offsets[npackets] is always equal to 'nbytes'.

    udp_packet_list(int max_npackets, int max_nbytes);

    // Accessors (not range-checked)
    inline uint8_t *get_packet_data(int i)  { return data_start + packet_offsets[i]; }
    inline int get_packet_nbytes(int i)     { return packet_offsets[i+1] - packet_offsets[i]; }

    // To add a packet, we copy its data to the udp_packet_list::data_end pointer, then call add_packet()
    // to update the rest of the udp_packet_list fields consistently.
    void add_packet(int packet_nbytes);

    // Doesn't deallocate buffers or change the max_* fields, but sets the current packet count to zero.
    void reset();
};


// High-level comment: the get/put methods of udp_packet_ringbuf have been designed so that a
// fixed pool of udp_packet_lists is recycled throughout the lifetime of the ringbuf, rather than
// having buffers which are continually freed and allocated.  This is to avoid the page-faulting
// cost of Linux malloc.

struct udp_packet_ringbuf : noncopyable {
    // Specified at construction, used when new udp_packet_list objects are allocated
    const int ringbuf_capacity;
    const int max_npackets_per_list;
    const int max_nbytes_per_list;

    pthread_mutex_t lock;
    pthread_cond_t cond_packets_added;
    pthread_cond_t cond_packets_removed;
    bool stream_ended = false;

    int ringbuf_size = 0;
    int ringbuf_pos = 0;
    std::vector<std::unique_ptr<udp_packet_list> > ringbuf;

    udp_packet_ringbuf(int ringbuf_capacity, int max_npackets_per_list, int max_nbytes_per_list);
    ~udp_packet_ringbuf();
    
    // Note!  The pointer 'p' is _swapped_ with an empty udp_packet_list from the ring buffer.
    // In other words, when put_packet_list() returns, the argument 'p' points to an empty udp_packet_list.
    // Returns true on success, returns false if packets were dropped due to full ring buffer.
    // Throws an exception if called after end-of-stream.
    bool put_packet_list(std::unique_ptr<udp_packet_list> &p, bool is_blocking);

    void get_size(int* currsize, int* maxsize);

    // Note!  The pointer 'p' is _swapped_ with the udp_packet_list which is extracted from the ring buffer.
    // In other words, when get_packet_list() returns, the original udp_packet_list will be "recycled" (rather than freed).
    // Returns true on success (possibly after blocking), returns false if ring buffer is empty and stream has ended.
    bool get_packet_list(std::unique_ptr<udp_packet_list> &p);

    // Called by producer thread, when stream has ended.
    void end_stream();

    // No longer used, but I left it in anyway.
    bool is_alive();
};


// -------------------------------------------------------------------------------------------------
//
// assembled_chunk_ringbuf
//
// A thread-safe data structure for exchanging assembled_chunks between the assembler and processing threads.  
// It also manages the "active" assembled_chunks, which are being filled with data as new packets arrive.
// There is one assembled_chunk_ringbuf for each beam.

class assembled_chunk_ringbuf : noncopyable,
                                public std::enable_shared_from_this<assembled_chunk_ringbuf> {
public:
    // When an assembled_chunk is put in the downstream ringbuf
    std::atomic<uint64_t> max_fpga_flushed;
    // When an assembled_chunk is retrieved from the downstream ringbuf
    // by, eg, rf_pipelines::chime_network_stream
    std::atomic<uint64_t> max_fpga_retrieved;
    
    assembled_chunk_ringbuf(const intensity_network_stream::initializer &ini_params, int beam_id, int stream_id);

    ~assembled_chunk_ringbuf();

    // Called by assembler thread, to "assemble" an intensity_packet into the appropriate assembled_chunk.
    // The length-(intensity_network_stream::event_type::num_types) event_counts array is incremented 
    // in the appropriate indices.
    //
    // The packet must have nbeams==1.  If the actual packet stream has multiple beams, then each packet
    // is split into multiple packets (in a zero-copy way, see intensity_network_stream.cpp) which are
    // passed to different assembled_chunk_ringbufs.
    //
    // The 'event_counts' array will be modified without acquiring any locks.  This is OK because the
    // assembler thread passes an event_subcounts array which is updated on a per-packet basis, and
    // accumulated into the global event_counts on a per-udp_packet_list basis (with locks acquired!).
    //
    // This routine is nonblocking.  If it would need to block, then it drops data somewhere and
    // increments the appropriate event_count.  (Depending which flags are set in ini_params, it may
    // also print a warning or throw an exception.)
    //
    // Warning: only safe to call from assembler thread!

    void put_unassembled_packet(const intensity_packet &packet, int64_t *event_counts);
    
    // Called by the assembler thread, when it exits.
    // Moves any remaining active chunks into the ring buffer, sets 'doneflag', initializes 'final_fpga'.
    void end_stream(int64_t *event_counts);

    void set_frame0(uint64_t frame0_nano);
    
    // Debugging: inject the given chunk
    bool inject_assembled_chunk(assembled_chunk* chunk);

    // For debugging/testing: stream data to disk.  
    //   'filename pattern': see assembled_chunk::format_filename (empty string to turn off streaming)
    //   'priority': see write_chunk_request::priority
    void stream_to_files(const std::string &filename_pattern, int priority,
                         bool need_rfi);

    // Callback from the write thread when a stream_to_files() write finishes
    void chunk_streamed(const std::string &filename);

    void get_streamed_chunks(int &achunks, size_t &abytes);

    // Debugging: print state
    void print_state();

    // Called by processing threads, via intensity_network_stream::get_assembled_chunk().
    // Returns the next assembled_chunk from the ring buffer, blocking if necessary to wait for data.
    // If the ring buffer is empty and end_stream() has been called, it returns an empty pointer
    // to indicate end-of-stream.
    std::shared_ptr<assembled_chunk> get_assembled_chunk(bool wait=true);

    // Find an assembled_chunk with the given fpgacounts start time, if it exists in the ring buffer.
    // Called by processing threads, in order to fill the RFI mask.
    // Returns an empty pointer iff stream has ended, and chunk is requested past end-of-stream.
    // If anything else goes wrong, an exception will be thrown.
    std::shared_ptr<assembled_chunk> find_assembled_chunk(uint64_t fpga_counts, bool top_level_only=false);
                                                          

    // The return value is a vector of (chunk, where) pairs, where 'where' is of type enum l1_ringbuf_level (defined in ch_frb_io.hpp)
    std::vector<std::pair<std::shared_ptr<assembled_chunk>, uint64_t>> get_ringbuf_snapshot(uint64_t min_fpga_counts=0, uint64_t max_fpga_counts=0);

    // Returns stats about the ring buffer, for the get_statistics RPC.
    //  *ringbuf_fpga_next* is the FPGA-counts of the next chunk that will be delivered to get_assembled_chunk().
    //  *ringbuf_n_ready* is the number of chunks available to be consumed by get_assembled_chunk().
    //  *ringbuf_capacity* is the maximum number of chunks that can be held in the ring buffer.
    //  *ringbuf_nelements* counts the number of valid chunks, including old chunks that have already been consumed by get_assembled_chunk.
    //  *ringbuf_fpga_min* is the smallest FPGA-counts number available in the ring buffer (including ones that have already been consumed by get_assembled_chunk().)
    //  *ringbuf_fpga_max* is the largest FPGA-counts number available in the ring buffer (including ones that have already been consumed by get_assembled_chunk().).  This includes the number of FPGA samples in the chunks.
    //
    // If *level* is specified, then returns *capacity*, *nelements*,
    // *fpga_min* and *fpga_max* for one level of the ringbuffer;
    // level=1 is the original intensity data, level=2 is
    // downsampled-by-2, level=3 is downsampled x 4, etc.
    void get_ringbuf_size(uint64_t* ringbuf_fpga_next,
                          uint64_t* ringbuf_n_ready,
                          uint64_t* ringbuf_capacity,
                          uint64_t* ringbuf_nelements,
                          uint64_t* ringbuf_fpga_min,
                          uint64_t* ringbuf_fpga_max,
                          int level=0);

protected:
    const intensity_network_stream::initializer ini_params;
    const int beam_id;
    const int stream_id;   // only used in assembled_chunk::format_filename().

    uint64_t frame0_nano; // nanosecond time() value for fgpacount zero
    
    output_device_pool output_devices;

    // Set to 'true' in the first call to put_unassembled_packet().
    bool first_packet_received = false;

    // Helper function called in assembler thread, to add a new assembled_chunk to the ring buffer.
    // Resets 'chunk' to a null pointer.
    // Warning: only safe to call from assembler thread.
    bool _put_assembled_chunk(std::unique_ptr<assembled_chunk> &chunk, int64_t *event_counts);

    // For debugging.
    // Warning: only safe to call from assembler thread.
    void _check_invariants();
    
    // Helper function: allocates new assembled chunk (from a memory_slab_pool, if one has been
    // specified in ini_params::memory_pool).  For now, an exception is thrown if the allocation
    // fails.  (FIXME: add code to recover gracefully.)
    std::unique_ptr<assembled_chunk> _make_assembled_chunk(uint64_t ichunk, int binning, bool zero=true);

    // The "active" chunks are in the process of being filled with data as packets arrive.
    // Currently we take the active window to be two assembled_chunks long, but this could be generalized.
    // When an active chunk is finished, it is added to the ring buffer.
    // Note: the active_chunk pointers are not protected by a lock, but are only accessed by the assembler thread.
    std::unique_ptr<assembled_chunk> active_chunk0;
    std::unique_ptr<assembled_chunk> active_chunk1;

    // Not sure if this really affects bottom-line performance, but thought it would be a good idea
    // to ensure that the "assembler-only" and "shared" fields were on different cache lines.
    char pad[constants::cache_line_size];

    // All fields below are protected by the lock
    pthread_mutex_t lock;

    // Processing thread waits here if the ring buffer is empty.
    pthread_cond_t cond_assembled_chunks_added;
    
    // Telescoping ring buffer.
    // All ringbuf* vectors have length num_downsampling_levels.
    // ringbuf[i] is a vector of length ringbuf_capacity[i]

    int num_downsampling_levels;
    std::vector<int> ringbuf_pos;
    std::vector<int> ringbuf_size;
    std::vector<int> ringbuf_capacity;
    std::vector<std::vector<std::shared_ptr<assembled_chunk>>> ringbuf;

    int downstream_pos;      // Position of "downstream" thread in ringbuf[0]
    int downstream_bufsize;  // Buffering capacity (in assembled_chunks) between assembler and downstream.

    inline std::shared_ptr<assembled_chunk> &ringbuf_entry(int ids, int ipos)
    {
	return ringbuf[ids][ipos % ringbuf_capacity[ids]];
    }

    // Are we streaming data to disk?  (Note: these fields require the lock for either read or write access.)
    std::string stream_pattern;
    int stream_priority = 0;
    bool stream_rfi_mask = false;
    int stream_chunks_written = 0;
    size_t stream_bytes_written = 0;

    bool doneflag = false;
    uint64_t final_fpga = 0;   // last fpga count which has an assembled_chunk, only initialized when 'doneflag' is set to true.
};


// -------------------------------------------------------------------------------------------------
//
// Some internal building blocks for downsampling kernels.  
// These are externally callable for use in unit tests.


extern void ds_slow_kernel1(float *out_data, int *out_mask, const uint8_t *in_data, 
			    const float *in_offsets, const float *in_scales, float *out_count, 
			    float *out_mean, int nupfreq, int nt_per_chunk, int nt_per_packet);

extern void ds_slow_kernel2(const float *in_data, const int *in_mask, float *w0, float *w1, 
			    float *w2, int nupfreq, int nt_per_chunk, int nt_per_packet);

extern void ds_slow_kernel3(uint8_t *out, const float *data, const int *mask, const float *enc_off, 
			    const float *enc_scal, int nupfreq, int nt_per_chunk, int nt_per_packet);


// -------------------------------------------------------------------------------------------------
//
// Miscelleanous inlines


inline bool is_power_of_two(int n)
{
    if (n <= 0)
	throw std::runtime_error("ch_frb_io: internal error: is_power_of_two() received argument <= 0");
    return (n & (n-1)) == 0;
}

inline int round_down_to_power_of_two(int n)
{
    if (n <= 0)
	throw std::runtime_error("ch_frb_io: internal error: is_power_of_two() received argument <= 0");
    return 1 << (int)log2(n+0.5);
}

template<typename T> 
inline bool equal3(T x, T y, T z)
{
    return (x == y) && (y == z);
}

inline double uniform_rand(std::mt19937 &rng)
{
    return std::uniform_real_distribution<>()(rng);
}

inline double uniform_rand(std::mt19937 &rng, double lo, double hi)
{
    return lo + (hi-lo) * uniform_rand(rng);
}

template<typename T> 
inline void uniform_rand(std::mt19937 &rng, T *p, int n)
{
    for (int i = 0; i < n; i++)
	p[i] = uniform_rand(rng);
}

template<typename T> 
inline void uniform_rand(std::mt19937 &rng, T *p, int n, double lo, double hi)
{
    for (int i = 0; i < n; i++)
	p[i] = uniform_rand(rng, lo, hi);
}

inline int randint(std::mt19937 &rng, int lo, int hi)
{
    return std::uniform_int_distribution<>(lo,hi-1)(rng);   // note hi-1 here!
}

// Returns true if string 's1' is a prefix of 's2'.
inline bool is_prefix(const std::string &s1, const std::string &s2)
{
    const char *p1 = s1.c_str();
    const char *p2 = s2.c_str();

    for(;;) {
	if (*p1 == 0)
            return true;
	if (*p1 != *p2)
	    return false;
	p1++;
	p2++;
    }
}

inline bool file_exists(const std::string &filename)
{
    struct stat s;

    int err = stat(filename.c_str(), &s);
    if (err >= 0)
        return true;
    if (errno == ENOENT)
        return false;

    throw std::runtime_error(filename + ": " + strerror(errno));
}

inline void hard_link(const std::string &src_filename, const std::string &dst_filename)
{
    int err = link(src_filename.c_str(), dst_filename.c_str());
    if (err < 0)
	throw std::runtime_error("hard_link() failed (" + src_filename + " -> " + dst_filename + ": " + strerror(errno));
}

template<typename T> inline T sum(const std::vector<T> &v)
{
    T ret = T(0);
    for (unsigned int i = 0; i < v.size(); i++)
	ret += v[i];
    return ret;
}

template<typename T> inline T prod(const std::vector<T> &v)
{
    T ret = (T)1;
    for (unsigned int i = 0; i < v.size(); i++)
	ret *= v[i];
    return ret;
}

// returns string representation of a vector
template<typename T> inline std::string vstr(const T *buf, int n, int stride=1)
{
    std::stringstream ss;
    ss << "[";
    for (int i = 0; i < n; i++)
	ss << " " << buf[i*stride];
    ss << " ]";
    return ss.str();
}

template<typename T> inline std::string vstr(const std::vector<T> &buf)
{
    return vstr(&buf[0], buf.size());
}


template<typename T>
inline bool vcontains(const std::vector<T> &v, T x)
{
    for (size_t i = 0; i < v.size(); i++) {
	if (v[i] == x)
	    return true;
    }

    return false;
}


template<typename T>
inline T *aligned_alloc(size_t nelts)
{
    if (nelts == 0)
	return NULL;

    // align to 64-byte cache lines
    void *p = NULL;
    if (posix_memalign(&p, 64, nelts * sizeof(T)) != 0)
	throw std::runtime_error("couldn't allocate memory");

    memset(p, 0, nelts * sizeof(T));
    return reinterpret_cast<T *> (p);
}

template<typename T>
inline std::unique_ptr<T[]> aligned_unique_ptr(size_t nelts)
{
    return std::unique_ptr<T[]> (aligned_alloc<T> (nelts));
}

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&& ...args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}


inline struct timeval xgettimeofday()
{
    struct timeval tv;

    int err = gettimeofday(&tv, NULL);
    if (_unlikely(err))
	throw std::runtime_error("gettimeofday failed");

    return tv;
}

inline int64_t usec_between(const struct timeval &tv1, const struct timeval &tv2)
{
    return 1000000 * int64_t(tv2.tv_sec - tv1.tv_sec) + int64_t(tv2.tv_usec - tv1.tv_usec);
}

inline void xusleep(useconds_t usec)
{
    int err = usleep(usec);
    if (err)
	throw std::runtime_error("usleep failed");
}


// -------------------------------------------------------------------------------------------------
//
// HDF5 wrappers (these are generally useful since the libhdf5 C/C++ api is so clunky)


template<typename T> inline hid_t hdf5_type();

// Reference: https://www.hdfgroup.org/HDF5/doc/H5.user/Datatypes.html
template<> inline hid_t hdf5_type<int>()            { return H5T_NATIVE_INT; }
template<> inline hid_t hdf5_type<unsigned long>() { return H5T_NATIVE_ULONG; }
template<> inline hid_t hdf5_type<unsigned long long>() { return H5T_NATIVE_ULLONG; }
template<> inline hid_t hdf5_type<float>()          { return H5T_NATIVE_FLOAT; }
template<> inline hid_t hdf5_type<double>()         { return H5T_NATIVE_DOUBLE; }
template<> inline hid_t hdf5_type<unsigned char>()  { return H5T_NATIVE_UCHAR; }


struct hdf5_file : noncopyable {
    std::string filename;
    hid_t file_id;

    // If write=false, the file is opened read-only, and an exception is thrown if it doesn't exist.
    // If write=true, the file is opened for writing.  If the file already exists, it will either be clobbered
    // or an exception will be thrown, depending on the value of 'clobber'.
    hdf5_file(const std::string &filename, bool write=false, bool clobber=true);
    ~hdf5_file();
};


struct hdf5_group : noncopyable {
    std::string filename;
    std::string group_name;
    hid_t group_id;

    // If create=true, the group will be created if it doesn't exist.
    hdf5_group(const hdf5_file &f, const std::string &group_name, bool create=false);
    ~hdf5_group();

    bool has_attribute(const std::string &attr_name) const;
    bool has_dataset(const std::string &dataset_name) const;

    void get_attribute_shape(const std::string &attr_name, std::vector<hsize_t> &shape) const;
    void get_dataset_shape(const std::string &attr_name, std::vector<hsize_t> &shape) const;

    // Read scalar attribute
    template<typename T> T read_attribute(const std::string &attr_name) const
    {
	T ret;
	this->_read_attribute(attr_name, hdf5_type<T>(), reinterpret_cast<void *> (&ret), std::vector<hsize_t>());
	return ret;
    }

    // Write scalar attribute
    template<typename T> void write_attribute(const std::string &attr_name, const T &x)
    {
	this->_write_attribute(attr_name, hdf5_type<T>(), reinterpret_cast<const void *> (&x), std::vector<hsize_t>());
    }

    // Write 1D vector attribute
    template<typename T> void write_attribute(const std::string &attr_name, const std::vector<T> &x)
    {
	std::vector<hsize_t> shape(1, x.size());
	this->_write_attribute(attr_name, hdf5_type<T>(), reinterpret_cast<const void *> (&x[0]), shape);
    }

    // Read multidimensional dataset
    template<typename T> void read_dataset(const std::string &dataset_name, T *out, const std::vector<hsize_t> &expected_shape) const
    {
	this->_read_dataset(dataset_name, hdf5_type<T>(), reinterpret_cast<void *> (out), expected_shape);
    }
    
    // Write multidimensional dataset
    template<typename T> void write_dataset(const std::string &dataset_name, const T *data, const std::vector<hsize_t> &shape)
    {
	this->_write_dataset(dataset_name, hdf5_type<T>(), reinterpret_cast<const void *> (data), shape);
    }

    // This interface is intended for small string-valued datasets.
    void write_string_dataset(const std::string &dataset_name, const std::vector<std::string> &data, const std::vector<hsize_t> &shape);
    void read_string_dataset(const std::string &dataset_name, std::vector<std::string> &data, const std::vector<hsize_t> &expected_shape) const;

    // Helpers
    void _get_attribute_shape(const std::string &attr_name, hid_t attr_id, std::vector<hsize_t> &shape) const;
    void _read_attribute(const std::string &attr_name, hid_t hdf5_type, void *out, const std::vector<hsize_t> &expected_shape) const;
    void _write_attribute(const std::string &attr_name, hid_t hdf5_type, const void *data, const std::vector<hsize_t> &shape);
    void _get_dataset_shape(const std::string &dataset_name, hid_t dataset_id, std::vector<hsize_t> &shape) const;
    void _check_dataset_shape(const std::string &dataset_name, hid_t dataset_id, const std::vector<hsize_t> &expected_shape) const;
    void _read_dataset(const std::string &dataset_name, hid_t hdf5_type, void *out, const std::vector<hsize_t> &expected_shape) const;
    void _write_dataset(const std::string &dataset_name, hid_t hdf5_type, const void *data, const std::vector<hsize_t> &shape);
};


// This class isn't intended to be used directly; use the wrapper hdf5_extendable_dataset<T> below
struct _hdf5_extendable_dataset : noncopyable {
    std::string filename;
    std::string group_name;
    std::string dataset_name;
    std::string full_name;
    std::vector<hsize_t> curr_shape;
    int axis;

    hid_t type;
    hid_t dataset_id;

    _hdf5_extendable_dataset(const hdf5_group &g, const std::string &dataset_name, 
			     const std::vector<hsize_t> &chunk_shape, int axis, hid_t type, int bitshuffle);

    ~_hdf5_extendable_dataset();

    void write(const void *data, const std::vector<hsize_t> &shape);
};


template<typename T>
struct hdf5_extendable_dataset {
    _hdf5_extendable_dataset base;

    //
    // The 'bitshuffle' argument has the following meaning:
    //   0 = no compression
    //   1 = try to compress, but if plugin fails then just write uncompressed data instead
    //   2 = try to compress, but if plugin fails then print a warning and write uncompressed data instead
    //   3 = compression mandatory
    //
    // List of all filter_ids: https://www.hdfgroup.org/services/contributions.html
    // Note that the compile-time constant 'bitshuffle_id' (=32008) is defined above.
    //
    hdf5_extendable_dataset(const hdf5_group &g, const std::string &dataset_name, const std::vector<hsize_t> &chunk_shape, int axis, int bitshuffle=0) :
	base(g, dataset_name, chunk_shape, axis, hdf5_type<T>(), bitshuffle)
    { }

    void write(const T *data, const std::vector<hsize_t> &shape)
    {
	base.write(data, shape);
    }

};


}  // namespace ch_frb_io

#endif // _CH_FRB_IO_INTERNALS_HPP
