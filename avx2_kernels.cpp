// This file contains assembly language kernels for packet assembly and decoding.
// These kernels improve bottom-line performance of the network front end by ~30%.  
// It might be possible to improve further, using memory access optimizations such
// as streaming writes and aligned loads/stores (reference: chapter 7 of the Intel 
// optimization manual).
//
// I didn't bother writing assembly language kernels for packet _encoding_, since
// this code is used only for testing and performance isn't as critical.
//
// The assembly language kernels make assumptions on the packet parameters:
// nt_per_packet must be equal to 16, and nupfreq must be even.  According to
// a recent email from Andre, these assumptions may not hold up, in which case
// we may need to write more kernels!
//
// The fast kernels are called using the following mechanism.  The member functions
// assembled_chunk::add_packet(), assembled_chunk::decode() are virtual.  When an
// assembled_chunk is allocated, we test whether the packet parameters (nt_per_packet,
// nupfreq) satisfy the constraints which permit fast kernels.  If so, then we allocate
// a special subclass fast_assembled_chunk which overrides the virtuals and calls the
// fast kernels.  Thus, this source file is responsible for defining the functions
// fast_assembled_chunk::add_packet() and fast_assembled_chunk::decode().

#include <cassert>
#include <iomanip>
#include <algorithm>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif

#ifndef __AVX2__

// If compiling on a machine without the AVX2 instruction set, we include some placeholder routines

fast_assembled_chunk::fast_assembled_chunk(int beam_id_, int nupfreq_, int nt_per_packet_, int fpga_counts_per_sample_, uint64_t ichunk_) :
    assembled_chunk(beam_id_, nupfreq_, nt_per_packet_, fpga_counts_per_sample_, ichunk_)
{
    throw runtime_error("ch_frb_io: internal error: fast_assembled_chunk constructor called on a non-AVX2 machine");
}

void fast_assembled_chunk::add_packet(const intensity_packet &packet)
{
    throw runtime_error("ch_frb_io: internal error: fast_assembled_chunk::add_packet() called on a non-AVX2 machine");
}

void fast_assembled_chunk::decode(float *intensity, float *weights, int stride) const
{
    throw runtime_error("ch_frb_io: internal error: fast_assembled_chunk::decode() called on a non-AVX2 machine");
}

void test_avx2_kernels(std::mt19937 &rng)
{
    cerr << "test_avx2_kernels(): this machine does not have the AVX2 instruction set, nothing to do\n";
}

#else  // __AVX2__

#include <immintrin.h>

// -------------------------------------------------------------------------------------------------
//
// Utils


// Given 128-bit SIMD register which holds 4 single-precision floats, extract the N-th float
template<unsigned int N>
inline float _extract128(__m128 x)
{
    // _mm_extract_ps() returns int instead of float?!
    union { int i; float x; } u;
    u.i = _mm_extract_ps(x, N);
    return u.x;
}

// Given 256-bit SIMD register which holds 8 single-precision floats, extract the N-th float
template<unsigned int N>
inline float _extract256(__m256 x)
{
    __m128 x2 = _mm256_extractf128_ps(x, N/4);
    return _extract128<N%4> (x2);
}

// These functions return string representations of SIMD registers, for debugging:
//   _vstr8_partial:  256-bit SIMD register which holds 32 8-bit integers
//   _vstr32_partial: 256-bit SIMD register which holds 8 32-bit integers
//   _vstrps_partial: 256-bit SIMD register which holds 8 single-precision floats

template<unsigned int N> inline void _vstr8_partial(stringstream &ss, __m256i x, bool hexflag);
template<unsigned int N> inline void _vstr32_partial(stringstream &ss, __m256i x, bool hexflag);
template<unsigned int N> inline void _vstrps_partial(stringstream &ss, __m256 x);

template<> inline void _vstr8_partial<0>(stringstream &ss, __m256i x, bool hex) { return; }
template<> inline void _vstr32_partial<0>(stringstream &ss, __m256i x, bool hex) { return; }
template<> inline void _vstrps_partial<0>(stringstream &ss, __m256 x) { return; }

template<unsigned int N> 
inline void _vstr8_partial(stringstream &ss, __m256i x, bool hexflag) 
{
    _vstr8_partial<N-1>(ss, x, hexflag);
    if (hexflag)
	ss << " " << setfill('0') << setw(2) << hex << uint32_t(uint8_t(_mm256_extract_epi8(x,N-1)));
    else
	ss << " " << int32_t(_mm256_extract_epi8(x,N-1));
}

template<unsigned int N>
inline void _vstr32_partial(stringstream &ss, __m256i x, bool hexflag) 
{
    _vstr32_partial<N-1>(ss, x, hexflag);
    if (hexflag)
	ss << " " << setfill('0') << setw(8) << hex << uint32_t(_mm256_extract_epi32(x,N-1));
    else
	ss << " " << _mm256_extract_epi32(x,N-1);
}

template<unsigned int N>
inline void _vstrps_partial(stringstream &ss, __m256 x)
{
    _vstrps_partial<N-1>(ss, x);
    ss << " " << _extract256<N-1>(x);
}


inline string _vstr8(__m256i x, bool hexflag=true)
{
    stringstream ss;
    ss << "[";
    _vstr8_partial<32> (ss, x, hexflag);
    ss << " ]";
    return ss.str();
}

inline string _vstr32(__m256i x, bool hexflag=false)
{
    stringstream ss;
    ss << "[";
    _vstr32_partial<8> (ss, x, hexflag);
    ss << " ]";
    return ss.str();
}

inline string _vstr(__m256 x)
{
    stringstream ss;
    ss << "[";
    _vstrps_partial<8> (ss, x);
    ss << " ]";
    return ss.str();
}


// -------------------------------------------------------------------------------------------------
//
// add_packet_kernel()
//
// Just copies src -> dst, where both 'src' and 'dst' are logical 2D arrays of shape (nupfreq, 16).
// The src array has frequency stride 16, as appropriate for a "close-packed" array.
// The dst array has frequency stride nt_per_assembled_chunk, as appopriate for an assembled_chunk subarray.
// The kernel assumes nt_per_packet=16, and nupfreq is even.


inline void _add_packet_kernel(uint8_t *dst, const uint8_t *src, int nupfreq)
{
    constexpr int s = constants::nt_per_assembled_chunk;

    for (int i = 0; i < nupfreq; i += 2) {
	__m256i x = _mm256_loadu_si256((const __m256i *) (src + 16*i));
	__m128i x0 = _mm256_extractf128_si256(x, 0);
	__m128i x1 = _mm256_extractf128_si256(x, 1);
	
	_mm_storeu_si128((__m128i *) (dst + i*s), x0);
	_mm_storeu_si128((__m128i *) (dst + (i+1)*s), x1);
    }
}


// -------------------------------------------------------------------------------------------------
//
// Decode kernels
//
// FIXME: these can be made faster for sure!


// Input: a 256-bit SIMD register which holds 32 8-bit unsigned integers [ x0, ..., x31 ].
// Output: four 256-bit SIMD registers which each hold 8 32-bit unsigned integers
//   out0 = [ x0, ..., x7 ]
//   out1 = [ x8, ..., x15 ]
//   out2 = [ x16, ..., x23 ]
//   out3 = [ x24, ..., x31 ]

inline void _decode_unpack(__m256i &out0, __m256i &out1, __m256i &out2, __m256i &out3, __m256i x)
{
    static const __m256i ctl0 = _mm256_set_epi8(15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0,
						15,11,7,3,14,10,6,2,13,9,5,1,12,8,4,0);

    // 4-by-4 transpose within each 128-bit lane
    x = _mm256_shuffle_epi8(x, ctl0);

    __m256i y0 = _mm256_and_si256(x, _mm256_set1_epi32(0xff));
    __m256i y1 = _mm256_and_si256(x, _mm256_set1_epi32(0xff00));
    __m256i y2 = _mm256_and_si256(x, _mm256_set1_epi32(0xff0000));
    __m256i y3 = _mm256_and_si256(x, _mm256_set1_epi32(0xff000000));

    y1 = _mm256_srli_epi32(y1, 8);
    y2 = _mm256_srli_epi32(y2, 16);
    y3 = _mm256_srli_epi32(y3, 24);

    out0 = _mm256_permute2f128_si256(y0, y1, 0x20);
    out1 = _mm256_permute2f128_si256(y2, y3, 0x20);
    out2 = _mm256_permute2f128_si256(y0, y1, 0x31);
    out3 = _mm256_permute2f128_si256(y2, y3, 0x31);
}


// Input: a 256-bit register which holds 8 32-bit unsigned integers
// Output: 8 single-precision floating point values.
// Each output value is 1.0 if the corresponding input satisfies 1 <= x <= 254, or 0.0 otherwise.

inline void _decode_weights(float *wtp, __m256i x, __m256i i0, __m256i i254, __m256 f0, __m256 f1)
{
    __m256i gt0 = _mm256_cmpgt_epi32(x, i0);
    __m256i gt254 = _mm256_cmpgt_epi32(x, i254);
    __m256i valid = _mm256_andnot_si256(gt254, gt0);
    __m256  wt = _mm256_blendv_ps(f0, f1, (__m256)valid);

    _mm256_storeu_ps(wtp, wt);
}

// Decodes 32 bytes of 8-bit intensity data, obtaining 32 floating-point intensities and
// 32 floating-point weights.  The kernel assumes nt_per_packet = 16.  Thus there are two
// relevant scales 'scale0', 'scale1' and two offsets 'offset0', 'offset1'.
//
// 'data' is a 256-bit SIMD register holding the input data.
//
// scale0, scale1, offset0, offset1 are scalar values, redundantly packed eight times into
// a 256-bit SIMD register.

inline void _decode_kernel32(float *intp, float *wtp, __m256i data, __m256 scale0, __m256 scale1, __m256 offset0, __m256 offset1)
{
    __m256i in0, in1, in2, in3;
    _decode_unpack(in0, in1, in2, in3, data);
    
    _mm256_storeu_ps(intp, scale0 * _mm256_cvtepi32_ps(in0) + offset0);
    _mm256_storeu_ps(intp+8, scale0 * _mm256_cvtepi32_ps(in1) + offset0);
    _mm256_storeu_ps(intp+16, scale1 * _mm256_cvtepi32_ps(in2) + offset1);
    _mm256_storeu_ps(intp+24, scale1 * _mm256_cvtepi32_ps(in3) + offset1);

    __m256i i0 = _mm256_set1_epi32(0);
    __m256i i254 = _mm256_set1_epi32(254);
    __m256 f0 = _mm256_set1_ps(0.0);
    __m256 f1 = _mm256_set1_ps(1.0);

    _decode_weights(wtp, in0, i0, i254, f0, f1);
    _decode_weights(wtp+8, in1, i0, i254, f0, f1);
    _decode_weights(wtp+16, in2, i0, i254, f0, f1);
    _decode_weights(wtp+24, in3, i0, i254, f0, f1);
}

// Decodes 128 bytes of 8-bit intensity data.

inline void _decode_kernel128(float *intp, float *wtp, const uint8_t *data, const float *scalep, const float *offsetp)
{
    __m256 scale = _mm256_loadu_ps(scalep);
    __m256 offset = _mm256_loadu_ps(offsetp);
    __m256 scale0, offset0;

    scale0 = _mm256_permute2f128_ps(scale, scale, 0x00);
    offset0 = _mm256_permute2f128_ps(offset, offset, 0x00);

    _decode_kernel32(intp, wtp, 
		     _mm256_loadu_si256((const __m256i *) (data)),
		     _mm256_shuffle_ps(scale0, scale0, 0x00), 
		     _mm256_shuffle_ps(scale0, scale0, 0x55),
		     _mm256_shuffle_ps(offset0, offset0, 0x00), 
		     _mm256_shuffle_ps(offset0, offset0, 0x55));

    _decode_kernel32(intp + 32, wtp + 32, 
		     _mm256_loadu_si256((const __m256i *) (data + 32)),
		     _mm256_shuffle_ps(scale0, scale0, 0xaa), 
		     _mm256_shuffle_ps(scale0, scale0, 0xff),
		     _mm256_shuffle_ps(offset0, offset0, 0xaa), 
		     _mm256_shuffle_ps(offset0, offset0, 0xff));


    scale0 = _mm256_permute2f128_ps(scale, scale, 0x11);
    offset0 = _mm256_permute2f128_ps(offset, offset, 0x11);

    _decode_kernel32(intp + 64, wtp + 64,
		     _mm256_loadu_si256((const __m256i *) (data + 64)),
		     _mm256_shuffle_ps(scale0, scale0, 0x00), 
		     _mm256_shuffle_ps(scale0, scale0, 0x55),
		     _mm256_shuffle_ps(offset0, offset0, 0x00), 
		     _mm256_shuffle_ps(offset0, offset0, 0x55));

    _decode_kernel32(intp + 96, wtp + 96, 
		     _mm256_loadu_si256((const __m256i *) (data + 96)),
		     _mm256_shuffle_ps(scale0, scale0, 0xaa), 
		     _mm256_shuffle_ps(scale0, scale0, 0xff),
		     _mm256_shuffle_ps(offset0, offset0, 0xaa), 
		     _mm256_shuffle_ps(offset0, offset0, 0xff));
}


// Decodes one "row" (i.e. nt_per_assembled_chunk time values) of 8-bit intensity data.

inline void _decode_kernel(float *intp, float *wtp, const uint8_t *data, const float *scalep, const float *offsetp)
{
    static_assert(constants::nt_per_assembled_chunk % 128 == 0, "_decode_kernel() assumes nt_per_assembled_chunk divisible by 128");

    constexpr int n = constants::nt_per_assembled_chunk / 128;

    for (int i = 0; i < n; i++)
	_decode_kernel128(intp + i*128, wtp + i*128, data + i*128, scalep + i*8, offsetp + i*8);
}


// -------------------------------------------------------------------------------------------------
//
// class fast_assembled_chunk


fast_assembled_chunk::fast_assembled_chunk(int beam_id_, int nupfreq_, int nt_per_packet_, int fpga_counts_per_sample_, uint64_t ichunk_) :
    assembled_chunk(beam_id_, nupfreq_, nt_per_packet_, fpga_counts_per_sample_, ichunk_)
{
    if (nt_per_packet_ != 16)
	throw runtime_error("ch_frb_io: internal error: fast_assembled_chunk constructor called with nt_per_packet != 16");
    if (nupfreq % 2 != 0)
	throw runtime_error("ch_frb_io: internal error: fast_assembled_chunk constructor called with odd value of nupfreq");
}


// virtual override
void fast_assembled_chunk::add_packet(const intensity_packet &packet)
{
    // Offset relative to beginning of packet
    uint64_t t0 = packet.fpga_count / uint64_t(fpga_counts_per_sample) - isample;

    for (int f = 0; f < packet.nfreq_coarse; f++) {
	int coarse_freq_id = packet.coarse_freq_ids[f];

	int d = coarse_freq_id*nt_coarse + (t0/nt_per_packet);
	this->scales[d] = packet.scales[f];
	this->offsets[d] = packet.offsets[f];

	uint8_t *dst = data + coarse_freq_id * nupfreq * constants::nt_per_assembled_chunk + t0;
	const uint8_t *src = packet.data + f * nupfreq * 16;

	_add_packet_kernel(dst, src, nupfreq);
    }
}


// virtual override
void fast_assembled_chunk::decode(float *intensity, float *weights, int stride) const
{
    if (!intensity || !weights)
	throw runtime_error("ch_frb_io: null pointer passed to fast_assembled_chunk::decode()");
    if (stride < constants::nt_per_assembled_chunk)
	throw runtime_error("ch_frb_io: bad stride passed to fast_assembled_chunk::decode()");

    for (int if_coarse = 0; if_coarse < constants::nfreq_coarse_tot; if_coarse++) {
	const float *scales_f = this->scales + if_coarse * nt_coarse;
	const float *offsets_f = this->offsets + if_coarse * nt_coarse;
	
	for (int if_fine = if_coarse*nupfreq; if_fine < (if_coarse+1)*nupfreq; if_fine++) {
	    const uint8_t *src_f = this->data + if_fine * constants::nt_per_assembled_chunk;
	    float *int_f = intensity + if_fine * stride;
	    float *wt_f = weights + if_fine * stride;

	    _decode_kernel(int_f, wt_f, src_f, scales_f, offsets_f);
	}
    }    
}


// -------------------------------------------------------------------------------------------------
//
// Unit testing


// helper function used by test_avx2_kernels()
static vector<float> randvec(std::mt19937 &rng, ssize_t n)
{
    vector<float> ret(n);
    uniform_rand(rng, &ret[0], n);
    return ret;
}


void peek_at_unpack_kernel()
{
    uint8_t v[32];
    for (int i = 0; i < 32; i++)
	v[i] = i+1;

    __m256i x = _mm256_loadu_si256((const __m256i *) v);

    __m256i y0, y1, y2, y3;
    _decode_unpack(y0, y1, y2, y3, x);

    cout << _vstr8(x,false) << endl
	 << _vstr32(y0,false) << endl
	 << _vstr32(y1,false) << endl
	 << _vstr32(y2,false) << endl
	 << _vstr32(y3,false) << endl;
}


void test_avx2_kernels(std::mt19937 &rng)
{
    cerr << "test_avx2_kernels()";

    // Required by fast decode kernel
    const int nt_per_packet = 16;

    // Initialized arbitrarily, since decode() doesn't use them.
    const int beam_id = 0;
    const int fpga_counts_per_sample = 384;
    const uint64_t ichunk = 0;

    for (int iouter = 0; iouter < 128; iouter++) {
	cerr << ".";

	// Randomized in every iteration
	const int nupfreq = 2 * randint(rng, 1, 9);
	const int nfreq_coarse_per_packet = 1 << randint(rng, 0, 6);
	const int stride = randint(rng, constants::nt_per_assembled_chunk, constants::nt_per_assembled_chunk + 16);

	// Set up intensity_packet

	vector<uint16_t> beam_ids = { beam_id };
	vector<uint16_t> coarse_freq_ids(nfreq_coarse_per_packet, 0);
	vector<float> scales(nfreq_coarse_per_packet, 0.0);
	vector<float> offsets(nfreq_coarse_per_packet, 0.0);
	vector<uint8_t> packet_data(nfreq_coarse_per_packet * nupfreq * nt_per_packet, 0);

	intensity_packet p;
	p.protocol_version = 1;
	p.data_nbytes = nfreq_coarse_per_packet * nupfreq * nt_per_packet;
	p.fpga_counts_per_sample = fpga_counts_per_sample;
	p.fpga_count = 0;
	p.nbeams = 1;
	p.nfreq_coarse = nfreq_coarse_per_packet;
	p.nupfreq = nupfreq;
	p.ntsamp = nt_per_packet;
	p.beam_ids = &beam_ids[0];
	p.coarse_freq_ids = &coarse_freq_ids[0];
	p.scales = &scales[0];
	p.offsets = &offsets[0];
	p.data = &packet_data[0];

	// Some auxiliary data which is useful when simulating random packets
	
	const int nt_coarse = constants::nt_per_assembled_chunk / nt_per_packet;
	const int npacket_max = (constants::nfreq_coarse_tot / nfreq_coarse_per_packet) * nt_coarse;
	const int npackets = randint(rng, npacket_max/2, npacket_max+1);

	vector<int> coarse_freq_ids_remaining(nt_coarse, constants::nfreq_coarse_tot);
	vector<uint16_t> coarse_freq_id_pool(nt_coarse * constants::nfreq_coarse_tot);

	for (int i = 0; i < nt_coarse; i++) {
	    for (int j = 0; j < constants::nfreq_coarse_tot; j++)
		coarse_freq_id_pool[i*constants::nfreq_coarse_tot + j] = j;

	    std::shuffle(coarse_freq_id_pool.begin() + i*constants::nfreq_coarse_tot,
			 coarse_freq_id_pool.begin() + (i+1)*constants::nfreq_coarse_tot,
			 rng);
	}

	// Test 1: Equivalence of assembled_chunk::add_packet() and fast_assembled_chunk::add_packet()
	
	auto chunk0 = make_shared<assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);
	auto chunk1 = make_shared<fast_assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);

	for (int ipacket = 0; ipacket < npackets; ipacket++) {
	    // Simulate random packet
	    int it_coarse;

	    do {
		it_coarse = randint(rng, 0, nt_coarse);
	    } while (coarse_freq_ids_remaining[it_coarse] < nfreq_coarse_per_packet);

	    coarse_freq_ids_remaining[it_coarse] -= nfreq_coarse_per_packet;

	    const uint16_t *src_ids = &coarse_freq_id_pool[0] + it_coarse*constants::nfreq_coarse_tot + coarse_freq_ids_remaining[it_coarse];
	    memcpy(p.coarse_freq_ids, src_ids, nfreq_coarse_per_packet * sizeof(uint16_t));
	    
	    uniform_rand(rng, p.scales, nfreq_coarse_per_packet);
	    uniform_rand(rng, p.offsets, nfreq_coarse_per_packet);

	    for (unsigned int i = 0; i < packet_data.size(); i++) {
		// Assign ~10% probability to 0x00 or 0xff
		int x = randint(rng, -25, 281);
		x = max(x, 0);
		x = min(x, 255);
		p.data[i] = uint8_t(x);
	    }

	    chunk0->add_packet(p);
	    chunk1->add_packet(p);
	}

	for (int i = 0; i < chunk0->nscales; i++) {
	    assert(chunk0->scales[i] == chunk1->scales[i]);
	    assert(chunk0->offsets[i] == chunk1->offsets[i]);
	}

	for (int i = 0; i < chunk0->ndata; i++)
	    assert(chunk0->data[i] == chunk1->data[i]);

	// Test 2: Equivalence of assembled_chunk::decode() and fast_assembled_chunk::decode()

	chunk0->randomize(rng);
	chunk1->fill_with_copy(chunk0);

	int nfreq_fine = constants::nfreq_coarse_tot * nupfreq;
	vector<float> intensity0 = randvec(rng, nfreq_fine * stride);
	vector<float> intensity1 = randvec(rng, nfreq_fine * stride);
	vector<float> weights0 = randvec(rng, nfreq_fine * stride);
	vector<float> weights1 = randvec(rng, nfreq_fine * stride);

	chunk0->decode(&intensity0[0], &weights0[0], stride);
	chunk1->decode(&intensity1[0], &weights1[0], stride);

	for (int ifreq = 0; ifreq < nfreq_fine; ifreq++) {
	    for (int it = 0; it < constants::nt_per_assembled_chunk; it++) {
		int i = ifreq*stride + it;
		int j = ifreq*constants::nt_per_assembled_chunk + it;

		if (fabs(intensity0[i] - intensity1[i]) > 1.0e-5) {
		    cerr << "\n " << nupfreq << " " << ifreq << " " << it 
			 << " " << int32_t(chunk0->data[j]) << " " << int32_t(chunk1->data[j])
			 << " " << intensity0[i] << " " << intensity1[i] << endl;
		    throw runtime_error("test_avx2_kernels: intensity mismatch");
		}

		if (weights0[i] != weights1[i]) {
		    cerr << "\n " << nupfreq << " " << ifreq << " " << it 
			 << " " << int32_t(chunk0->data[j]) << " " << int32_t(chunk1->data[j])
			 << " " << weights0[i] << " " << weights1[i] << endl;
		    throw runtime_error("test_avx2_kernels: weights mismatch");
		}
	    }
	}
    }

    cerr << "success\n";
}

#endif  // __AVX2__

}  // namespace ch_frb_io
