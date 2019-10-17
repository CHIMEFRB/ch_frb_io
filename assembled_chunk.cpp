#include <sys/types.h>
#include <sys/stat.h>
#include <cstdio>
#include <iostream>
#include <stdarg.h> // for va_start/va_end
#include <unistd.h>
#include <stdio.h>
#include <cinttypes>  // PRIu64
#include <immintrin.h>
#include <msgpack/fbuffer.hpp>
#include "assembled_chunk_msgpack.hpp"
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


// -------------------------------------------------------------------------------------------------
//
// Constructor and allocation logic.


// This helper function's purpose in life is to compute
//   nt_coarse = constants::nt_per_assembled_chunk / nt_per_packet
// in a way which throws a C++ exception (rather than a hardware exception)
// if nt_per_packet == 0.

static int _nt_c(int nt_per_packet)
{
    int nt_f = constants::nt_per_assembled_chunk;

    if (nt_per_packet <= 0)
	throw runtime_error("ch_frb_io: assembled_chunk constructor expects nt_per_packet > 0");
    if (nt_f % nt_per_packet != 0)
	throw runtime_error("ch_frb_io: assembled_chunk constructor expects nt_per_packet to be a divisor of " + to_string(nt_f));
    return nt_f / nt_per_packet;
}


// A helper class describing the layout of assembled_chunk::memory_slab.
struct memory_slab_layout {
    const int nfreq_c;
    const int nfreq_f;
    const int nt_f;
    const int nt_c;

    // All nb_* fields are array sizes in bytes.
    const int nb_data;
    const int nb_scales;
    const int nb_offsets;
    const int nb_rfimask;
    const int nb_ds_data;
    const int nb_ds_mask;
    const int nb_ds_w2;
    const int nb_detrend_t;
    const int nb_detrend_f;

    // All ib_* fields are array offsets within the memory_slab, in bytes
    const int ib_data;
    const int ib_scales;
    const int ib_offsets;
    const int ib_rfimask;
    const int ib_ds_data;
    const int ib_ds_mask;
    const int ib_ds_w2;
    const int ib_detrend_t;
    const int ib_detrend_f;
    const int slab_size;

    static int align(int nbytes) { return ((nbytes+63)/64) * 64; }

    memory_slab_layout(int nupfreq, int nt_per_packet, int nrfifreq,
                       int n_detrend_t, int n_detrend_f) :
	nfreq_c(constants::nfreq_coarse_tot),
	nfreq_f(constants::nfreq_coarse_tot * nupfreq),
	nt_f(constants::nt_per_assembled_chunk),
	nt_c(_nt_c(nt_per_packet)),
	nb_data(nfreq_f * nt_f),
	nb_scales(nfreq_c * nt_c * sizeof(float)),
	nb_offsets(nfreq_c * nt_c * sizeof(float)),
        nb_rfimask(nrfifreq * nt_f / 8 * sizeof(uint8_t)),
	nb_ds_data(nupfreq * (nt_f/2) * sizeof(float)),
	nb_ds_mask(nupfreq * (nt_f/2) * sizeof(int)),
	nb_ds_w2((nt_c/2) * sizeof(float)),
        nb_detrend_t(nt_f * n_detrend_t * sizeof(float)),
        nb_detrend_f(nfreq_f * n_detrend_f * sizeof(float)),
	ib_data(0),
	ib_scales(align(ib_data + nb_data)),
	ib_offsets(align(ib_scales + nb_scales)),
        ib_rfimask(align(ib_offsets + nb_offsets)),
	ib_ds_data(align(ib_rfimask + nb_rfimask)),
	ib_ds_mask(align(ib_ds_data + nb_ds_data)),
	ib_ds_w2(align(ib_ds_mask + nb_ds_mask)),
        ib_detrend_t(align(ib_ds_w2 + nb_ds_w2)),
        ib_detrend_f(align(ib_detrend_t + nb_detrend_t)),
	slab_size(align(ib_detrend_f + nb_detrend_f))
    { }
};


assembled_chunk::assembled_chunk(const assembled_chunk::initializer &ini_params) :
    beam_id(ini_params.beam_id), 
    nupfreq(ini_params.nupfreq),
    nrfifreq(ini_params.nrfifreq),
    nt_per_packet(ini_params.nt_per_packet),
    fpga_counts_per_sample(ini_params.fpga_counts_per_sample), 
    binning(ini_params.binning),
    stream_id(ini_params.stream_id),
    ichunk(ini_params.ichunk),
    frame0_nano(ini_params.frame0_nano),
    nt_coarse(_nt_c(nt_per_packet)),
    nscales(constants::nfreq_coarse_tot * nt_coarse),
    ndata(constants::nfreq_coarse_tot * nupfreq * constants::nt_per_assembled_chunk),
    nrfimaskbytes(nrfifreq * constants::nt_per_assembled_chunk / 8),
    n_detrend_t(ini_params.n_detrend_t),
    n_detrend_f(ini_params.n_detrend_f),
    isample(ichunk * constants::nt_per_assembled_chunk),
    fpga_begin(ichunk * constants::nt_per_assembled_chunk * fpga_counts_per_sample),
    fpga_end((ichunk+binning) * constants::nt_per_assembled_chunk * fpga_counts_per_sample),
    has_rfi_mask(false),
    has_detrend_t(false),
    has_detrend_f(false),
    packets_received(0)
{
    if ((beam_id < 0) || (beam_id > constants::max_allowed_beam_id))
	throw runtime_error("assembled_chunk constructor: bad 'beam_id' argument");
    if ((nupfreq <= 0) || (nupfreq > constants::max_allowed_nupfreq))
	throw runtime_error("assembled_chunk constructor: bad 'nupfreq' argument");
    if (nrfifreq < 0)
	throw runtime_error("assembled_chunk constructor: bad 'nrfifreq' argument");
    if ((nrfifreq > 0) && ((constants::nfreq_coarse_tot * nupfreq) % nrfifreq))
	throw runtime_error("assembled_chunk constructor: bad 'nrfifreq' argument");
    if ((nt_per_packet <= 0) || !is_power_of_two(nt_per_packet) || (nt_per_packet > constants::nt_per_assembled_chunk))
	throw runtime_error("assembled_chunk constructor: bad 'nt_per_packet' argument");
    if ((fpga_counts_per_sample <= 0) || (fpga_counts_per_sample > constants::max_allowed_fpga_counts_per_sample))
	throw runtime_error("assembled_chunk constructor: bad 'fpga_counts_per_sample' argument");
    if ((binning <= 0) || !is_power_of_two(binning) || (binning > 8))
	throw runtime_error("assembled_chunk constructor: bad 'binning' argument");
    if ((stream_id < 0) || (stream_id > 9))
	throw runtime_error("assembled_chunk constructor: bad 'stream_id' argument");

    uint64_t ichunk_max = UINT64_MAX / uint64_t(constants::nt_per_assembled_chunk * fpga_counts_per_sample);
    if (ichunk > ichunk_max)
	throw runtime_error("assembled_chunk constructor: bad 'ichunk' argument");

    memory_slab_layout mc(nupfreq, nt_per_packet, nrfifreq,
                          n_detrend_t, n_detrend_f);

    if (ini_params.pool) {
	if (!ini_params.slab)
	    throw runtime_error("assembled_chunk constructor: 'pool' is a nonempty pointer, but 'slab' is empty");

	if (ini_params.pool->nbytes_per_slab < mc.slab_size) {
	    throw runtime_error("assembled_chunk constructor: memory_slab_pool::nbytes_per_slab (=" 
				+ to_string(ini_params.pool->nbytes_per_slab) 
				+ ") is less than required slab size (="
				+ to_string(mc.slab_size) + ")");
	}

	this->memory_pool = ini_params.pool;
	this->memory_slab.swap(ini_params.slab);
    }
    else {
	if (ini_params.slab)
	    throw runtime_error("assembled_chunk constructor: 'pool' is an empty pointer, but 'slab' is nonempty");

	uint8_t *p = aligned_alloc<uint8_t> (mc.slab_size);
	this->memory_slab = memory_slab_t(p);
    }

    this->data = memory_slab.get() + mc.ib_data;
    this->scales = reinterpret_cast<float *> (memory_slab.get() + mc.ib_scales);
    this->offsets = reinterpret_cast<float *> (memory_slab.get() + mc.ib_offsets);
    this->rfi_mask = nrfifreq ? reinterpret_cast<uint8_t *> (memory_slab.get() + mc.ib_rfimask) : nullptr;
    this->detrend_params_t = n_detrend_t ? reinterpret_cast<float *> (memory_slab.get() + mc.ib_detrend_t) : nullptr;
    this->detrend_params_f = n_detrend_f ? reinterpret_cast<float *> (memory_slab.get() + mc.ib_detrend_f) : nullptr;
    this->ds_data = reinterpret_cast<float *> (memory_slab.get() + mc.ib_ds_data);
    this->ds_mask = reinterpret_cast<int *> (memory_slab.get() + mc.ib_ds_mask);
    this->ds_w2 = reinterpret_cast<float *> (memory_slab.get() + mc.ib_ds_w2);
}


assembled_chunk::~assembled_chunk()
{
    this->_deallocate();
}


// Helper function for destructors.
void assembled_chunk::_deallocate()
{
    if (memory_pool) {
	memory_pool->put_slab(memory_slab);
	memory_pool = shared_ptr<memory_slab_pool> ();
    }
    
    // Shouldn't be necessary, but guards against accidental pointer reuse after free().
    this->scales = nullptr;
    this->offsets = nullptr;
    this->data = nullptr;
    this->rfi_mask = nullptr;
    this->detrend_params_t = nullptr;
    this->detrend_params_f = nullptr;
    this->ds_w2 = nullptr;
    this->ds_data = nullptr;
    this->ds_mask = nullptr;
}

    
// Static member function
ssize_t assembled_chunk::get_memory_slab_size(int nupfreq, int nt_per_packet, int nrfifreq, int n_detrend_t, int n_detrend_f)
{
    if ((nupfreq <= 0) || (nupfreq > constants::max_allowed_nupfreq))
	throw runtime_error("assembled_chunk::get_memory_slab_size(): bad 'nupfreq' argument");
    if ((nt_per_packet <= 0) || !is_power_of_two(nt_per_packet) || (nt_per_packet > constants::nt_per_assembled_chunk))
	throw runtime_error("assembled_chunk::get_memory_slab_size(): bad 'nt_per_packet' argument");
    if (nrfifreq < 0)
	throw runtime_error("assembled_chunk::get_memory_slab_size(): bad 'nrfifreq' argument");
    if ((nrfifreq > 0) && ((constants::nfreq_coarse_tot * nupfreq) % nrfifreq))
	throw runtime_error("assembled_chunk::get_memory_slab_size(): bad 'nrfifreq' argument");

    memory_slab_layout mc(nupfreq, nt_per_packet, nrfifreq, n_detrend_t, n_detrend_f);
    return mc.slab_size;
}


// -------------------------------------------------------------------------------------------------

// Nanoseconds per FPGA count
static const uint64_t fpga_nano = 2560;

double assembled_chunk::time_begin() const {
    return (frame0_nano + fpga_counts_per_sample * fpga_nano * fpga_begin) * 1e-9;
}

double assembled_chunk::time_end() const {
    return (frame0_nano + fpga_counts_per_sample * fpga_nano * fpga_end) * 1e-9;}

static string 
__attribute__ ((format(printf,1,2)))
stringprintf(const char* format, ...) {
    va_list lst;
    /* Yarrrr vasprintf is not in the C++ standard.
    int rtn;
     va_start(lst, format);
     char* strp = NULL;
     rtn = vasprintf(strp, format, lst);
     if (rtn == -1)
     throw runtime_error("stringprintf failed: " + string(strerror(errno)));
     va_end(lst);
     string s(strp);
     free(strp);
     */
    char temps[256];
    va_start(lst, format);
    // truncates if length > size of 'temps'
    if (vsnprintf(temps, sizeof(temps), format, lst) < 0)
        throw runtime_error("stringprintf failed: " + string(strerror(errno)));
    va_end(lst);
    return string(temps);
}

// Replaces all instances of the string "from" to the string "to" in
// input string "input".
static string replaceAll(const string &input, const string &from, const string &to) {
    string s = input;
    size_t i;
    while ((i = s.find(from)) != std::string::npos)
        s.replace(i, from.length(), to);
    return s;
}

string assembled_chunk::format_filename(const string &pattern) const {
    //   (STREAM)  -> %01i stream_id
    //   (BEAM)    -> %04i beam_id
    //   (CHUNK)   -> %08i ichunk
    //   (NCHUNK)  -> %02i  size in chunks
    //   (BINNING) -> %02i  size in chunks
    //   (FPGA0)   -> %012i start FPGA-counts
    //   (FPGAN)   -> %08i  FPGA-counts size
    string s = pattern;
    s = replaceAll(s, "(STREAM)",  stringprintf("%01i",        stream_id));
    s = replaceAll(s, "(BEAM)",    stringprintf("%04i",        beam_id));
    s = replaceAll(s, "(CHUNK)",   stringprintf("%08"  PRIu64, ichunk));
    s = replaceAll(s, "(NCHUNK)",  stringprintf("%02i",        binning));
    s = replaceAll(s, "(BINNING)", stringprintf("%02i",        binning));
    s = replaceAll(s, "(FPGA0)",   stringprintf("%012" PRIu64, fpga_begin));
    s = replaceAll(s, "(FPGAN)",   stringprintf("%08"  PRIu64, fpga_end - fpga_begin));
    return s;
}

void assembled_chunk::fill_with_copy(const shared_ptr<assembled_chunk> &x)
{
    if (!x)
	throw runtime_error("assembled_chunk::fill_with_copy() called with empty pointer");
    if ((this->nupfreq != x->nupfreq) || (this->nt_per_packet != x->nt_per_packet))
	throw runtime_error("assembled_chunk::fill_with_copy() called on non-conformable chunks");
    if (this->nrfifreq != x->nrfifreq)
	throw runtime_error("assembled_chunk::fill_with_copy() called on non-conformable chunks (nrfifreq)");

    if (x.get() == this)
	return;

    memcpy(this->data, x->data, ndata);
    memcpy(this->scales, x->scales, nscales * sizeof(float));
    memcpy(this->offsets, x->offsets, nscales * sizeof(float));
    memcpy(this->rfi_mask, x->rfi_mask, nrfimaskbytes);
}


// Used in unit tests.
// If nrfifreq > 0, also randomizes the RFI mask and sets the 'has_rfi_mask' flag.
void assembled_chunk::randomize(std::mt19937 &rng)
{
    for (int i = 0; i < ndata; i++) {
	// Assign ~10% probability to 0x00 or 0xff
	int x = randint(rng, -25, 281);
	x = max(x, 0);
	x = min(x, 255);
	this->data[i] = uint8_t(x);
    }

    uniform_rand(rng, this->scales, nscales);
    uniform_rand(rng, this->offsets, nscales);

    if (nrfifreq > 0) {
	for (int i = 0; i < nrfimaskbytes; i++)
	    this->rfi_mask[i] = randint(rng, 0, 256);
	this->has_rfi_mask = true;
    }
}


// virtual member function; any changes made here should be reflected in override fast_assembled_chunk::add_packet().
void assembled_chunk::add_packet(const intensity_packet &packet)
{
    uint64_t packet_t0 = packet.fpga_count / uint64_t(fpga_counts_per_sample);

    // Offset relative to beginning of packet
    uint64_t t0 = packet_t0 - isample;
    
    // The runtime checks in intensity_network_stream::_process_packet() should
    // ensure that the following checks are redundant.  I decided to include the 
    // redundant checks here in the "generic" assembled_chunk::add_packet(), but 
    // omit them in fast_assembled_chunk::add_packet().

    bool bad = ((packet.nbeams != 1) ||
		(packet.nupfreq != this->nupfreq) ||
		(packet.ntsamp != this->nt_per_packet) ||
		(packet.fpga_counts_per_sample != this->fpga_counts_per_sample) ||
		(packet.fpga_count % (fpga_counts_per_sample * nt_per_packet)) ||
		(packet.beam_ids[0] != this->beam_id) ||
		(packet_t0 < isample) ||
		(packet_t0 + nt_per_packet > isample + constants::nt_per_assembled_chunk));

    if (_unlikely(bad))
	throw runtime_error("ch_frb_io: internal error in assembled_chunk::add_packet()");

    for (int f = 0; f < packet.nfreq_coarse; f++) {
	int coarse_freq_id = packet.coarse_freq_ids[f];

	this->scales[coarse_freq_id*nt_coarse + (t0/nt_per_packet)] = packet.scales[f];
	this->offsets[coarse_freq_id*nt_coarse + (t0/nt_per_packet)] = packet.offsets[f];

	for (int u = 0; u < nupfreq; u++) {
	    memcpy(data + (coarse_freq_id*nupfreq + u) * constants::nt_per_assembled_chunk + t0, 
		   packet.data + (f*nupfreq + u) * nt_per_packet,
		   nt_per_packet);
	}
    }
    this->packets_received++;
}


// virtual member function; any changes made here should be reflected in override fast_assembled_chunk::decode().
// Note: decode() and decode_subset() do not apply the RFI mask in assembled_chunk::rfi_mask (if this exists).

void assembled_chunk::decode(float *intensity, float *weights, int istride, int wstride, float prescale) const
{
    if (!intensity || !weights)
	throw runtime_error("ch_frb_io: null pointer passed to assembled_chunk::decode()");	
    if (istride < constants::nt_per_assembled_chunk)
	throw runtime_error("ch_frb_io: bad istride passed to assembled_chunk::decode()");
    if (wstride < constants::nt_per_assembled_chunk)
	throw runtime_error("ch_frb_io: bad wstride passed to assembled_chunk::decode()");

    for (int if_coarse = 0; if_coarse < constants::nfreq_coarse_tot; if_coarse++) {
	const float *scales_f = this->scales + if_coarse * nt_coarse;
	const float *offsets_f = this->offsets + if_coarse * nt_coarse;
	
	for (int if_fine = if_coarse*nupfreq; if_fine < (if_coarse+1)*nupfreq; if_fine++) {
	    const uint8_t *src_f = this->data + if_fine * constants::nt_per_assembled_chunk;
	    float *int_f = intensity + if_fine * istride;
	    float *wt_f = weights + if_fine * wstride;

	    for (int it_coarse = 0; it_coarse < nt_coarse; it_coarse++) {
		float scale = scales_f[it_coarse] * prescale;
		float offset = offsets_f[it_coarse] * prescale;
		
		for (int it_fine = it_coarse*nt_per_packet; it_fine < (it_coarse+1)*nt_per_packet; it_fine++) {
		    float x = float(src_f[it_fine]);
		    int_f[it_fine] = scale*x + offset;
		    wt_f[it_fine] = ((x==0) || (x==255)) ? 0.0 : 1.0;
		}
	    }
	}
    }
}

void assembled_chunk::decode_subset(float *intensity, float *weights, int t0, int NT, int istride, int wstride) const 
{                                    
    if (!intensity || !weights)
	throw runtime_error("ch_frb_io: null pointer passed to assembled_chunk::decode_subset()");
    if (istride < NT)
	throw runtime_error("ch_frb_io: bad istride passed to assembled_chunk::decode_subset()");
    if (wstride < NT)
	throw runtime_error("ch_frb_io: bad wstride passed to assembled_chunk::decode_subset()");
    if ((t0 < 0) || (NT < 0) || (t0 + NT > constants::nt_per_assembled_chunk))
	throw runtime_error("ch_frb_io: bad (t0,NT) passed to assembled_chunk::decode_subset()");

    for (int if_coarse = 0; if_coarse < constants::nfreq_coarse_tot; if_coarse++) {
	const float * scales_f = this->scales  + if_coarse * nt_coarse;
	const float *offsets_f = this->offsets + if_coarse * nt_coarse;

	for (int if_fine = if_coarse*nupfreq; if_fine < (if_coarse+1)*nupfreq; if_fine++) {
	    const uint8_t *src_f = this->data + if_fine * constants::nt_per_assembled_chunk;
	    float *int_f = intensity + if_fine * istride;
	    float * wt_f = weights   + if_fine * wstride;

            for (int i=0; i<NT; i++) {
                int it = t0 + i;
                int it_coarse = it / nt_per_packet;

		float scale  =  scales_f[it_coarse];
		float offset = offsets_f[it_coarse];

                float x = float(src_f[it]);
                int_f[i] = scale*x + offset;
                wt_f [i] = ((x==0) || (x==255)) ? 0.0 : 1.0;
	    }
	}
    }
}


static void ds_slow_kernel(uint8_t *out_data, float *out_offsets, float *out_scales, const uint8_t *in_data, 
			   const float *in_offsets, const float *in_scales, float *tmp_data, int *tmp_mask, 
			   float *tmp_scales, int nupfreq, int nt_per_chunk, int nt_per_packet)
{
    ds_slow_kernel1(tmp_data, tmp_mask, in_data, in_offsets, in_scales, out_offsets, out_scales, nupfreq, nt_per_chunk, nt_per_packet);
    ds_slow_kernel2(tmp_data, tmp_mask, out_offsets, out_scales, tmp_scales, nupfreq, nt_per_chunk, nt_per_packet);
    ds_slow_kernel3(out_data, tmp_data, tmp_mask, out_offsets, tmp_scales, nupfreq, nt_per_chunk, nt_per_packet);
}


// Reads 'nbits_in' bits (not bytes!), writes (nbits_in/2) bits.
// Assumes nbits_in is a multiple of 16.

static void ds_slow_kernel_rfimask(uint8_t *dst, const uint8_t *src, int nbits_in)
{
    memset(dst, 0, nbits_in/16);

    for (int i = 0; i < nbits_in/2; i++) {
	int idst = i / 8;
	int bdst = i % 8;
	int isrc = (2*i) / 8;
	int bsrc = (2*i) % 8;
	
	uint8_t x = src[isrc];
	uint8_t y = (1 << bsrc) | (1 << (bsrc+1));
	
	if ((x & y) == y)
	    dst[idst] |= (1 << bdst);
    }
}


// Checks validity of call to assembled_chunk::downsample().
void assembled_chunk::_check_downsample(const assembled_chunk *src1, const assembled_chunk *src2)
{
    if (!src1 || !src2)
	throw runtime_error("ch_frb_io: null pointer passed to assembled_chunk::downsample()");
    if (this->has_rfi_mask)
	throw runtime_error("ch_frb_io: assembled_chunk::downsample() called, and has_rfi_mask=true in destination chunk");

    if (src1->binning != src2->binning)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample(): mismatched binning");
    if (src2->ichunk != src1->ichunk + src2->binning)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample(): source ichunks are not consecutive");

    if (this->binning != 2 * src1->binning)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample(): wrong value of this->binning");
    if (this->ichunk != src1->ichunk)
	throw runtime_error("ch_frb_io: assembled_chunk::downsample(): wrong value of this->ichunk");

    if (!equal3(this->beam_id, src1->beam_id, src2->beam_id))
	throw runtime_error("ch_frb_io: assembled_chunk::downsample(): mismatched beam_id");
    if (!equal3(this->nupfreq, src1->nupfreq, src2->nupfreq))
        throw runtime_error("ch_frb_io: assembled_chunk::downsample(): mismatched nupfreq");
    if (!equal3(this->nrfifreq, src1->nrfifreq, src2->nrfifreq))
        throw runtime_error("ch_frb_io: assembled_chunk::downsample(): mismatched nrfifreq");
    if (!equal3(this->nt_per_packet, src1->nt_per_packet, src2->nt_per_packet))
        throw runtime_error("ch_frb_io: assembled_chunk::downsample(): mismatched nt_per_packet");

    if ((nrfifreq > 0) && (!src1->rfi_mask || !src2->rfi_mask || !src1->has_rfi_mask || !src2->has_rfi_mask))
        throw runtime_error("ch_frb_io: assembled_chunk::downsample(): RFI masks were not initialized as expected, maybe your ring buffer is too small?");
}

// This is the slow, reference version of assembled_chunk::downsample().
//
// In production, it is overridden by the fast, assembly-language-kernelized version
// in fast_assembled_chunk::downsample().  (See avx2_kernels.cpp)

void assembled_chunk::downsample(const assembled_chunk *src1, const assembled_chunk *src2)
{
    this->_check_downsample(src1, src2);

    int nfreq_c = constants::nfreq_coarse_tot;
    int nt_f = constants::nt_per_assembled_chunk;
    int nt_c = nt_f / nt_per_packet;
 
    for (int ifreq_c = 0; ifreq_c < nfreq_c; ifreq_c++) {
	int ifreq_f = ifreq_c * nupfreq;

	ds_slow_kernel(this->data + ifreq_f * nt_f,
		       this->offsets + ifreq_c * nt_c,
		       this->scales + ifreq_c * nt_c,
		       src1->data + ifreq_f * nt_f,
		       src1->offsets + ifreq_c * nt_c,
		       src1->scales + ifreq_c * nt_c,
		       this->ds_data, this->ds_mask, this->ds_w2,
		       nupfreq, nt_f, nt_per_packet);

	ds_slow_kernel(this->data + (ifreq_f * nt_f) + (nt_f/2),
		       this->offsets + (ifreq_c * nt_c) + (nt_c/2),
		       this->scales + (ifreq_c * nt_c) + (nt_c/2),
		       src2->data + ifreq_f * nt_f,
		       src2->offsets + ifreq_c * nt_c,
		       src2->scales + ifreq_c * nt_c,
		       this->ds_data, this->ds_mask, this->ds_w2,
		       nupfreq, nt_f, nt_per_packet);
    }

    if (nrfifreq <= 0)
	return;

    // Downsample RFI mask.
    // Note: if we get here, then _check_downsample() has already checked that src1, src2
    // have initialized RFI masks with the same value of 'nfreq'.

    for (int ifreq = 0; ifreq < nrfifreq; ifreq++) {
	ds_slow_kernel_rfimask(this->rfi_mask + ifreq * (nt_f/8),
			       src1->rfi_mask + ifreq * (nt_f/8),
			       nt_f);

	ds_slow_kernel_rfimask(this->rfi_mask + ifreq * (nt_f/8) + (nt_f/16),
			       src2->rfi_mask + ifreq * (nt_f/8),
			       nt_f);
    }

    this->has_rfi_mask = true;
}


unique_ptr<assembled_chunk> assembled_chunk::make(const assembled_chunk::initializer &ini_params)
{
    bool fast_kernel_exists = (ini_params.nt_per_packet == 16) && (ini_params.nupfreq % 2 == 0);

    if (ini_params.force_reference && ini_params.force_fast)
        throw runtime_error("ch_frb_io: assembled_chunk::make(): both force_reference and force_fast flags were set");
    if (ini_params.force_fast && !fast_kernel_exists)
	throw runtime_error("ch_frb_io: assembled_chunk::make(): force_fast flag was set, but conditions for a fast kernel (nt_per_packet=16 and nupfreq even) were not met");

#ifdef __AVX2__
    if (fast_kernel_exists && !ini_params.force_reference)
	return make_unique<fast_assembled_chunk> (ini_params);
#else
    if (ini_params.force_fast)
        throw runtime_error("ch_frb_io: assembled_chunk::make(): force_fast flag set on a machine without AVX2");
#endif

    return make_unique<assembled_chunk> (ini_params);
}


void assembled_chunk::write_hdf5_file(const string &filename)
{
    if (this->rfi_mask != nullptr)
	throw runtime_error("ch_frb_io: assembled_chunk hdf5 format does not implement RFI mask yet");
    
    bool write = true;
    bool clobber = true;
    hdf5_file f(filename, write, clobber);

    string chunkname = "/assembled-chunk-beam" + to_string(beam_id)
        + "-ichunk" + to_string(ichunk);
    bool create = true;
    hdf5_group g_chunk(f, chunkname, create);

    // Header
    g_chunk.write_attribute("beam_id", this->beam_id);
    g_chunk.write_attribute("nupfreq", this->nupfreq);
    g_chunk.write_attribute("nt_per_packet", this->nt_per_packet);
    g_chunk.write_attribute("fpga_counts_per_sample", this->fpga_counts_per_sample);
    g_chunk.write_attribute("nt_coarse", this->nt_coarse);
    g_chunk.write_attribute("nscales", this->nscales);
    g_chunk.write_attribute("ndata", this->ndata);
    g_chunk.write_attribute("ichunk", this->ichunk);
    g_chunk.write_attribute("isample", this->isample);

    // Offset & scale vectors
    vector<hsize_t> scaleshape = { (hsize_t)constants::nfreq_coarse_tot,
                                   (hsize_t)this->nt_coarse };
    g_chunk.write_dataset("scales",  this->scales,  scaleshape);
    g_chunk.write_dataset("offsets", this->offsets, scaleshape);

    // Raw data
    int bitshuffle = 0;
    vector<hsize_t> datashape = {
        (hsize_t)constants::nfreq_coarse_tot,
        (hsize_t)nupfreq,
        (hsize_t)constants::nt_per_assembled_chunk };
    unique_ptr<hdf5_extendable_dataset<uint8_t> > data_dataset =
        make_unique<hdf5_extendable_dataset<uint8_t> >(g_chunk, "data", datashape, 2, bitshuffle);
    data_dataset->write(this->data, datashape);
    // close
    data_dataset = unique_ptr<hdf5_extendable_dataset<uint8_t> > ();
}

void assembled_chunk::write_msgpack_file(const string &filename, bool compress, uint8_t *buffer)
{
    if ((this->rfi_mask != nullptr) && (!this->has_rfi_mask))
	throw runtime_error("ch_frb_io::assembled_chunk::write_msgpack_file() called on chunk whose RFI mask has not been initialized yet");
    
    char tempfilename[filename.size() + 10];
    sprintf(tempfilename, "%s.tmpXXXXXX", filename.c_str());
    int fd = mkstemp(tempfilename);
    if (fd == -1)
        throw runtime_error("ch_frb_io: failed to create temp file for " + filename + " for writing an assembled_chunk in msgpack format: " + strerror(errno));
    FILE* f = fdopen(fd, "w+");
    if (!f)
        throw runtime_error("ch_frb_io: failed to open temp file " + string(tempfilename) + " for writing an assembled_chunk in msgpack format: " + strerror(errno));
    // msgpack buffer that will write to file "f"
    msgpack::fbuffer fb(f);
    // Construct a shared_ptr from this, carefully
    shared_ptr<assembled_chunk> shthis(shared_ptr<assembled_chunk>(), this);
    msgpack::packer<msgpack::fbuffer> packer(fb);
    // The real deal: in assembled_chunk_msgpack.hpp
    pack_assembled_chunk(packer, shthis, compress, buffer);

    if (fclose(f))
        throw runtime_error("ch_frb_io: failed to close assembled_chunk msgpack temp file " + string(tempfilename) + ": " + string(strerror(errno)));

    // Set permissions -- mkstemp creates the file with mode 0600
    // the umask() call sets the mask to the given value and returns the previous
    // value -- thus the double-umask() call to leave it unchanged.
    mode_t mask = umask(0);
    umask(mask);
    if (chmod(tempfilename, 0666 & ~mask))
      throw runtime_error("ch_frb_io: failed to chmod temp file in writing assembled_chunk msgpack file: " + string(tempfilename) + ": " + string(strerror(errno)));

    if (rename(tempfilename, filename.c_str()))
        throw runtime_error("ch_frb_io: failed to rename temp file in writing assembled_chunk msgpack file: " + string(tempfilename) + ": " + string(strerror(errno)));
}

/*
// We could do this using msgpack::sbuffer, but if we're just going to write to a file, no real point first serializing into memory and then writing to disk.
 int assembled_chunk::serialize(uint8_t* buf, int capacity, bool compress) {
 }
 */

size_t assembled_chunk::max_compressed_size() {
    return bshuf_compress_lz4_bound(this->ndata, 1, 0);
}

shared_ptr<assembled_chunk> assembled_chunk::read_msgpack_file(const string &filename)
{
    struct stat st;
    if (stat(filename.c_str(), &st))
        throw runtime_error("ch_frb_io: failed to stat file " + filename + " for reading an assembled_chunk in msgpack format: " + strerror(errno));

    size_t len = st.st_size;
    FILE* f = fopen(filename.c_str(), "r");
    if (!f)
        throw runtime_error("ch_frb_io: failed to open file " + filename + " for reading an assembled_chunk in msgpack format: " + strerror(errno));

    unique_ptr<char> fdata(new char[len]);

    size_t nr = fread(fdata.get(), 1, len, f);
    if (nr != len)
        throw runtime_error("ch_frb_io: failed to read " + to_string(len) + " from file " + filename + " for reading an assembled_chunk in msgpack format: " + strerror(errno));
    fclose(f);

    msgpack::object_handle oh = msgpack::unpack(fdata.get(), len);
    msgpack::object obj = oh.get();
    shared_ptr<assembled_chunk> ch;
    obj.convert(ch);
    return ch;
}


// -------------------------------------------------------------------------------------------------
//
// Downsampling kernels
//
// The logic in assembled_chunk::downsample() is split into three kernels ds_slow_kernel[123]().
//
// This splitting parallels the fast AVX2 kernels, and allows the fast kernels to be unit-tested
// individually.  Otherwise, there's no reason to call the ds_slow_kernel*() functions directly,
// they are just helpers for assembled_chunk::downsample().


// ds_slow_kernel1(): downsample-and-decode kernel
//
// Assumes nt (=nt_per_chunk) is a multiple of (2*nt_per_packet).
// Assumes nt_per_packet is even.
//
// Reads (nupfreq, nt) input data values (uint8).
// Reads (nt/16) offset values (float32).
// Reads (nt/16) scale values (float32).
// Writes (nupfreq, nt/2) output data values (float32).
// Writes (nupfreq, nt/2) output mask values (int32).
// Writes (nt/32) count values (float32).
// Writes (nt/32) mean values (float32).
//
// Reminder: packets are decoded as
//   (intensity) = scale * (8-bit value) + offset


template<typename T> inline T square(T x) { return x*x; }


void ds_slow_kernel1(float *out_data, int *out_mask, const uint8_t *in_data, 
		     const float *in_offsets, const float *in_scales, float *out_count, 
		     float *out_mean, int nupfreq, int nt_per_chunk, int nt_per_packet)
{
    ch_assert(nt_per_chunk > 0);
    ch_assert(nt_per_packet > 0);
    ch_assert(nt_per_chunk % (2*nt_per_packet) == 0);
    ch_assert(nt_per_packet % 2 == 0);

    int np = nt_per_chunk / nt_per_packet;
    int nq = nt_per_packet / 2;

    memset(out_data, 0, nupfreq * (nt_per_chunk/2) * sizeof(float));
    memset(out_mask, 0, nupfreq * (nt_per_chunk/2) * sizeof(int));
    memset(out_count, 0, (np/2) * sizeof(float));
    memset(out_mean, 0, (np/2) * sizeof(float));

    for (int iupfreq = 0; iupfreq < nupfreq; iupfreq++) {
	for (int p = 0; p < np; p++) {
	    float offset = in_offsets[p];
	    float scale = 0.5 * in_scales[p];  // Note factor 0.5 here

	    for (int q = 0; q < nq; q++) {
		uint8_t d0 = in_data[iupfreq*nt_per_chunk + p*nt_per_packet + 2*q];
		uint8_t d1 = in_data[iupfreq*nt_per_chunk + p*nt_per_packet + 2*q+1];
		
		if ((d0 == 0x00) || (d0 == 0xff) || (d1 == 0x00) || (d1 == 0xff))
		    continue;

		float x = scale * (float(d0) + float(d1)) + offset;

		out_data[iupfreq*(nt_per_chunk/2) + p*nq + q] = x;
		out_mask[iupfreq*(nt_per_chunk/2) + p*nq + q] = -1;
		out_count[p/2] += 1.0;
		out_mean[p/2] += x;
	    }
	}
    }

    for (int i = 0; i < (np/2); i++)
	out_mean[i] /= max(out_count[i], 0.5f);
}


// ds_slow_kernel2(): computes variance, offsets, scales
//
// The 'in_data' and 'in_mask' arrays have shape (nupfreq,nt/2).  (Not shape (nupfreq,nt)!)
//
// The w0,w1,w2 arrays have length (nt / (2*nt_per_packet)).
//
//   - On input:
//       w0 = count
//       w1 = mean
//
//   - On output:
//       w0 = offset
//       w1 = dec_scale
//       w2 = enc_scale
//
// Reminder: packets are decoded and encoded as follows
//   (decoded intensity) = (dec_scale) * (8-bit value) + offset
//   (encoded 8-bit value) = (enc_scale) * (intensity - offset)

void ds_slow_kernel2(const float *in_data, const int *in_mask, float *w0, float *w1, float *w2, int nupfreq, int nt_per_chunk, int nt_per_packet)
{
    ch_assert(nt_per_chunk > 0);
    ch_assert(nt_per_packet > 0);
    ch_assert(nt_per_chunk % (2*nt_per_packet) == 0);

    int nt2 = nt_per_chunk/2;
    int np = nt_per_chunk / (2*nt_per_packet);

    memset(w2, 0, np * sizeof(float));

    for (int iupfreq = 0; iupfreq < nupfreq; iupfreq++) {
	for (int p = 0; p < np; p++) {
	    float mean = w1[p];
	    float var = 0.0;
	    
	    for (int i = p*nt_per_packet; i < (p+1)*nt_per_packet; i++)
		if (in_mask[iupfreq*nt2+i])
		    var += square(in_data[iupfreq*nt2+i] - mean);

	    w2[p] += var;
	}
    }

    for (int p = 0; p < np; p++) {
	float mean = w1[p];
	float var = w2[p] / max(w0[p],1.0f);

	var += 1.0e-10 * square(mean);
	
	float rms = (var > 0.0) ? sqrt(var) : 1.0;
	float dec_scale = 0.04 * rms;
	float enc_scale = 1.0 / dec_scale;
	float offset = mean - 128. * dec_scale;

	w0[p] = offset;
	w1[p] = dec_scale;
	w2[p] = enc_scale;
    }
}


// ds_slow_kernel3(): slow encode kernel.
//
// Note that the arrays here have (nt/2) time samples, not nt time samples! (where nt=nt_per_chunk)
//
// Reads (nupfreq,nt/2) data values (float32)
// Reads (nupfreq,nt/2) mask values (int32)
// Reads (nt/(2*nt_per_packet)) enc_offset values (float32)
// Reads (nt/(2*nt_per_packet)) enc_scale values (float32)
// Writes (nupfreq,nt/2) data values (uint8_t).

void ds_slow_kernel3(uint8_t *out, const float *data, const int *mask, const float *enc_off, const float *enc_scal, int nupfreq, int nt_per_chunk, int nt_per_packet)
{
    ch_assert(nt_per_chunk > 0);
    ch_assert(nt_per_packet > 0);
    ch_assert(nt_per_chunk % (2*nt_per_packet) == 0);

    int nt2 = nt_per_chunk / 2;
    int np = nt2 / nt_per_packet;

    for (int iupfreq = 0; iupfreq < nupfreq; iupfreq++) {
	for (int p = 0; p < np; p++) {
	    float off = enc_off[p];
	    float scal = enc_scal[p];

	    for (int i = p*nt_per_packet; i < (p+1)*nt_per_packet; i++) {
		float x = data[iupfreq*nt2 + i];
		int m = mask[iupfreq*nt2 + i];

		x = scal * (x - off);
		x = max(x, 0.0f);
		x = min(x, 255.0f);

		if (m != -1)
		    x = 0.0;

		x = roundf(x);
		out[iupfreq*nt_per_chunk + i] = uint8_t(x);
	    }
	}
    }
}


}  // namespace ch_frb_io
