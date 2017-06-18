// This file contains assembly language kernels for packet assembly and decoding.
// These kernels improve bottom-line performance of the network front end by ~30%.  
// It might be possible to improve further, using memory access optimizations such
// as streaming writes and aligned loads/stores (reference: chapter 7 of the Intel 
// optimization manual).
//
// I didn't bother writing assembly language kernels for packet _encoding_, since
// this code is used only for testing, and performance isn't as critical.
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
// Fast downsampling kernels start here, and continue for the next ~500 lines!


// _extract_16_8: a helper class for the downsampling kernels.
//
// Takes a uint16[32] array (represented as two __m256i's), extracts the low 8 bits from each element, 
// and outputs a uint8[32] array (represented at __m256i).

struct _extract_16_8 {
    const __m256i ctl0;
    const __m256i ctl1;

    _extract_16_8() :
	ctl0(_mm256_set_epi8(0x0e, 0x0c, 0x0a, 0x08, 0x06, 0x04, 0x02, 0x00,
			     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			     0x0e, 0x0c, 0x0a, 0x08, 0x06, 0x04, 0x02, 0x00,
			     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff)),
	ctl1(_mm256_set_epi8(0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			     0x0e, 0x0c, 0x0a, 0x08, 0x06, 0x04, 0x02, 0x00,
			     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			     0x0e, 0x0c, 0x0a, 0x08, 0x06, 0x04, 0x02, 0x00))
    { }

    inline __m256i extract(__m256i x, __m256i y)
    {
	x = _mm256_shuffle_epi8(x, ctl0);    // [ 0 x0 0 x1 ]
	y = _mm256_shuffle_epi8(y, ctl1);    // [ y0 0 y1 0 ]

	__m256i z = _mm256_blend_epi32(x, y, 0x33);  // [ y0 x0 y1 x1 ].  Note 0x33 = (00110011)_2

	// To finish the kernel, we just need to output [ x0 x1 y0 y1 ]

	__m256i u = _mm256_shuffle_epi32(z, 0x4e);          // [ x0 y0 x1 y1 ].  Note 0x4e = (1032)_4
	__m256i v = _mm256_permute2x128_si256(z, z, 0x01);  // [ y1 x1 y0 x0 ]
	
	return _mm256_blend_epi32(u, v, 0x3c);  // [ x0 x1 y0 y1 ]
    }
};


// _extract_32_8: a helper class for the downsampling kernels.
//
// Takes a uint32[32] array (represented as four __m256i's), extracts the low 8 bits from each 
// element, and outputs a uint8[32] array (represented as __m256i).

struct _extract_32_8 {
    const __m256i ctl0;
    const __m256i ctl1;
    unsigned int saved_rounding_mode;

    _extract_32_8() :
	ctl0(_mm256_set_epi8(0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			     0x0c, 0x08, 0x04, 0x00, 0xff, 0xff, 0xff, 0xff,
			     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			     0xff, 0xff, 0xff, 0xff, 0x0c, 0x08, 0x04, 0x00)),
	ctl1(_mm256_set_epi8(0x0c, 0x08, 0x04, 0x00, 0xff, 0xff, 0xff, 0xff,
			     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			     0xff, 0xff, 0xff, 0xff, 0x0c, 0x08, 0x04, 0x00,
			     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff))
    { 
	this->saved_rounding_mode = _MM_GET_ROUNDING_MODE();
	_MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
    }

    ~_extract_32_8()
    {
	_MM_SET_ROUNDING_MODE(saved_rounding_mode);
    }

    inline __m256i extract(__m256 a, __m256 b, __m256 c, __m256 d)
    {
	__m256 f0 = _mm256_setzero_ps();
	__m256 f255 = _mm256_set1_ps(255.0);

	a = _mm256_min_ps(_mm256_max_ps(a, f0), f255);
	b = _mm256_min_ps(_mm256_max_ps(b, f0), f255);
	c = _mm256_min_ps(_mm256_max_ps(c, f0), f255);
	d = _mm256_min_ps(_mm256_max_ps(d, f0), f255);
	
	__m256i ai = _mm256_cvtps_epi32(a);
	__m256i bi = _mm256_cvtps_epi32(b);
	__m256i ci = _mm256_cvtps_epi32(c);
	__m256i di = _mm256_cvtps_epi32(d);

	ai = _mm256_shuffle_epi8(ai, ctl0);  // [ a0 0 0 0 0 a1 0 0 ]
	bi = _mm256_shuffle_epi8(bi, ctl1);  // [ 0 0 b0 0 0 0 0 b1 ]
	ci = _mm256_shuffle_epi8(ci, ctl0);  // [ c0 0 0 0 0 c1 0 0 ]
	di = _mm256_shuffle_epi8(di, ctl1);  // [ 0 0 d0 0 0 0 0 d1 ]

	ai = _mm256_blend_epi32(ai, bi, 0xcc);    // [ a0 0 b0 0 0 a1 0 b1 ].  Note 0xcc = (11001100)_2
	ci = _mm256_blend_epi32(ci, di, 0xcc);    // [ c0 0 d0 0 0 c1 0 d1 ].
	
	__m256i t = _mm256_blend_epi32(ai, ci, 0xf0);   // [ a0 0 b0 0 0 c1 0 d1 ]
	__m256i u = _mm256_blend_epi32(ai, ci, 0x0f);   // [ c0 0 d0 0 0 a1 0 b1 ]

	u = _mm256_permute2x128_si256(u, u, 0x01);      // [ 0 a1 0 b1 c0 0 d0 0 ]
	
	return  _mm256_blend_epi32(t, u, 0x5a);  // [ a0 a1 b0 b1 c0 c1 d0 d1 ].  Note 0x5a = (01011010)_2
    }
};


// incremental_downsample8: a helper class for the downsampling kernels.
//
// Incrementally "absorbs" a float32[64], represented as eight __m256's, 
// downsamples float32[8], and returns the result as an __m256.

struct incremental_downsample8
{
    __m256 x0;
    __m256 x1;
    __m256 x2;

    static inline __m256 f(__m256 a, __m256 b) { return _mm256_shuffle_ps(a, b, 0x88) + _mm256_shuffle_ps(a, b, 0xdd); }

    // After add<7>(), the downsampled vector is 'x0'.
    template<int N> inline void add(__m256 x);

    // These only make sense to call after add<7>().
    inline void store(float *p)  { _mm256_store_ps(p, x0); }
    inline void update(float *p) { _mm256_store_ps(p, x0 + _mm256_load_ps(p)); }
};


template<> inline void incremental_downsample8::add<0> (__m256 x)  { x0 = x; }
template<> inline void incremental_downsample8::add<1> (__m256 x)  { x1 = f(x0,x); }
template<> inline void incremental_downsample8::add<2> (__m256 x)  { x0 = x; }
template<> inline void incremental_downsample8::add<3> (__m256 x)  { x2 = f(x1,f(x0,x)); }
template<> inline void incremental_downsample8::add<4> (__m256 x)  { x0 = x; }
template<> inline void incremental_downsample8::add<5> (__m256 x)  { x1 = f(x0,x); }
template<> inline void incremental_downsample8::add<6> (__m256 x)  { x0 = x; }

template<> inline void incremental_downsample8::add<7> (__m256 x)
{ 
    x0 = f(x1,f(x0,x));
    x0 = _mm256_blend_ps(x2, x0, 0xf0) + _mm256_permute2f128_ps(x2, x0, 0x21);
}


// incremental_upsample8: a helper class for the downsampling kernels.
//
// Loads a float32[8] from memory, and allows caller to request a constant
// __m256 containing any of the 8 floats (repeated 8 times).

struct incremental_upsample8 {
    __m256 x0;  // [x0 x0]
    __m256 x1;  // [x1 x1]
    
    incremental_upsample8(const float *p)
    {
	__m256 x01 = _mm256_loadu_ps(p);                      // [x0 x1]
	__m256 x10 = _mm256_permute2f128_ps(x01, x01, 0x01);  // [x1 x0]

	x0 = _mm256_blend_ps(x01, x10, 0xf0);  // [x0 x0]
	x1 = _mm256_blend_ps(x01, x10, 0x0f);  // [x1 x1]
    }

    template<int N> inline __m256 get();
};


template<> inline __m256 incremental_upsample8::get<0>() { return _mm256_permute_ps(x0, 0x00); }  // (0000)_4
template<> inline __m256 incremental_upsample8::get<1>() { return _mm256_permute_ps(x0, 0x55); }  // (1111)_4
template<> inline __m256 incremental_upsample8::get<2>() { return _mm256_permute_ps(x0, 0xaa); }  // (2222)_4
template<> inline __m256 incremental_upsample8::get<3>() { return _mm256_permute_ps(x0, 0xff); }  // (3333)_4
template<> inline __m256 incremental_upsample8::get<4>() { return _mm256_permute_ps(x1, 0x00); }
template<> inline __m256 incremental_upsample8::get<5>() { return _mm256_permute_ps(x1, 0x55); }
template<> inline __m256 incremental_upsample8::get<6>() { return _mm256_permute_ps(x1, 0xaa); }
template<> inline __m256 incremental_upsample8::get<7>() { return _mm256_permute_ps(x1, 0xff); }


// Helper function called by _kernel1a().
inline __m256 _unpack_16bit_data(__m128i x, __m256 offset, __m256 scale)
{
    __m256i y = _mm256_cvtepi16_epi32(x);
    __m256 z = _mm256_cvtepi32_ps(y);
    return scale * z + offset;
}


// Decode kernel.
// Reads 32 input data values (uint8).
// Reads 2 offset values, by calling offset.get<2N> and offset.get<2N+1>
// Reads 2 scale values, by calling scales.get<2N> and scales.get<2N+1>
// Writes 1 logical count value, by calling out_count.add<N>
// Writes 1 logical mean value, by calling out_mean.add<N>

template<int N>
inline void _ds_kernel1a(float *out_data, int *out_mask, const uint8_t *in_data,
			 incremental_upsample8 &in_offset, incremental_upsample8 &in_scale,
			 incremental_downsample8 &out_count, incremental_downsample8 &out_mean)
{
    static constexpr int N2 = (2*N) % 8;

    __m256i all_ones = _mm256_set1_epi8(0xff);

    // 8-bit input data
    __m256i d8 = _mm256_loadu_si256((const __m256i *) in_data);

    // 8-bit mask (actually the mask complement: 0x00 means "OK", 0xff means "invalid".
    __m256i m8a = _mm256_cmpeq_epi8(d8, _mm256_setzero_si256());
    __m256i m8b = _mm256_cmpeq_epi8(d8, all_ones);
    __m256i m8c = _mm256_or_si256(m8a, m8b);
    
    // 16-bit output data
    __m256i d16a = _mm256_and_si256(d8, _mm256_set1_epi16(0xff));
    __m256i d16b = _mm256_srli_epi16(d8, 8);
    __m256i d16 = _mm256_add_epi16(d16a, d16b);

    // 16-bit mask
    __m256i m16 = m8c;
    m16 = _mm256_or_si256(m16, _mm256_slli_epi16(m8c,8));  // shift left
    m16 = _mm256_or_si256(m16, _mm256_srli_epi16(m8c,8));  // shift right
    m16 = _mm256_xor_si256(m16, all_ones);

    // Unpack output data to pair (float32[8], float32[8])
    __m256 x0 = _unpack_16bit_data(_mm256_extracti128_si256(d16,0), in_offset.get<N2>(), in_scale.get<N2>());
    __m256 x1 = _unpack_16bit_data(_mm256_extracti128_si256(d16,1), in_offset.get<N2+1>(), in_scale.get<N2+1>());

    // Unpack 16-bit output mask to pair (int32[8], int32[8])
    __m256i m0 = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(m16,0));
    __m256i m1 = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(m16,1));

    // May as well apply mask here, since we need to do it before computing 'out_mean' anyway.
    x0 = _mm256_and_ps(x0, _mm256_castsi256_ps(m0));
    x1 = _mm256_and_ps(x1, _mm256_castsi256_ps(m1));

    // Write out_data
    _mm256_storeu_ps(out_data, x0);
    _mm256_storeu_ps(out_data+8, x1);

    // Write out_mask
    _mm256_storeu_si256((__m256i *) (out_mask), m0);
    _mm256_storeu_si256((__m256i *) (out_mask+8), m1);

    // Write out_count.
    __m256i count = _mm256_add_epi32(m0, m1);
    count = _mm256_sub_epi32(_mm256_setzero_si256(), count);
    out_count.add<N> (_mm256_cvtepi32_ps(count));

    // Write out_mean
    out_mean.add<N> (x0+x1);
}


// Reads 256 input data values (uint8).
// Reads 16 scale values (float32)
// Reads 16 offset values (float32)
// Writes 128 output data values (float32, scale/offset applied).
// Writes 128 output mask values (int32).
// Updates (not writes!) 8 mean values
// Updates (not writes!) 8 count values

inline void _ds_kernel1b(float *out_data, int *out_mask, const uint8_t *in_data, 
			 const float *in_offset, const float *in_scale, 
			 float *update_count, float *update_mean)
{
    incremental_downsample8 count;
    incremental_downsample8 mean;

    incremental_upsample8 offset0(in_offset);
    incremental_upsample8 scale0(in_scale);

    _ds_kernel1a<0> (out_data, out_mask, in_data, offset0, scale0, count, mean);
    _ds_kernel1a<1> (out_data+16, out_mask+16, in_data+32, offset0, scale0, count, mean);
    _ds_kernel1a<2> (out_data+32, out_mask+32, in_data+64, offset0, scale0, count, mean);
    _ds_kernel1a<3> (out_data+48, out_mask+48, in_data+96, offset0, scale0, count, mean);

    incremental_upsample8 offset1(in_offset+8);
    incremental_upsample8 scale1(in_scale+8);

    _ds_kernel1a<4> (out_data+64, out_mask+64, in_data+128, offset1, scale1, count, mean);
    _ds_kernel1a<5> (out_data+80, out_mask+80, in_data+160, offset1, scale1, count, mean);
    _ds_kernel1a<6> (out_data+96, out_mask+96, in_data+192, offset1, scale1, count, mean);
    _ds_kernel1a<7> (out_data+112, out_mask+112, in_data+224, offset1, scale1, count, mean);

    count.update(update_count);
    mean.update(update_mean);
}


// Fast decode kernel.
// Assumes nt (=nt_per_assembled_chunk) is a multiple of 256.
// Assumes nt_per_packet = 16.
//
// Reads (nupfreq, nt) input data values (uint8).
// Reads (nt/16) offset values (float32).
// Reads (nt/16) scale values (float32).
// Writes (nupfreq, nt/2) output data values (float32).
// Writes (nupfreq, nt/2) output mask values (int32).
// Writes (nt/32) count values (float32).
// Writes (nt/32) mean values (float32).
//
// Note that values of the 'out_data' array are undefined (but guaranteed to be non-NaN)
// if the corresponding 'out_mask' entry is zero.

inline void _ds_kernel1(float *out_data, int *out_mask, const uint8_t *in_data, const float *in_offsets, 
			const float *in_scales, float *out_count, float *out_mean, int nupfreq, int nt)
{
    // The first pass computes out_count = sum(w_i), but
    // uses out_mean as a temp buffer to store sum(w_i d_i).

    memset(out_count, 0, (nt/32) * sizeof(float));
    memset(out_mean, 0, (nt/32) * sizeof(float));

    for (int iupfreq = 0; iupfreq < nupfreq; iupfreq++) {
	for (int i = 0; i < (nt/256); i++)
	    _ds_kernel1b(out_data + iupfreq*(nt/2) + 128*i,
			 out_mask + iupfreq*(nt/2) + 128*i,
			 in_data + iupfreq*nt + 256*i,
			 in_offsets + 16*i,
			 in_scales + 16*i,
			 out_count + 8*i,
			 out_mean + 8*i);
    }

    // The second pass computes out_mean = sum(w_i d_i) / sum(w_i).

    __m256 one = _mm256_set1_ps(1.0);

    for (int i = 0; i < (nt/32); i += 8) {
	__m256 num = _mm256_loadu_ps(out_mean + i);
	__m256 den = _mm256_loadu_ps(out_count + i);
	__m256 mean = num / _mm256_max_ps(den, one);

	_mm256_storeu_ps(out_mean + i, mean);
    }
}


// Reads 16 data values (float32)
// Reads 16 mask values (int32)
// Calls in_mean::get<N>() to get the mean.
// Calls out_var::add<N>() with sum w d^2

template<int N>
inline void _ds_kernel2a(const float *data, const int *mask, incremental_upsample8 &in_mean, incremental_downsample8 &out_var)
{
    __m256 x0 = _mm256_loadu_ps(data);
    __m256 x1 = _mm256_loadu_ps(data+8);

    // Subtract offset.
    __m256 offset = in_mean.get<N>();
    x0 -= offset;
    x1 -= offset;

    // Apply mask
    const float *fmask = (const float *) (mask);
    x0 = _mm256_and_ps(x0, _mm256_loadu_ps(fmask));
    x1 = _mm256_and_ps(x1, _mm256_loadu_ps(fmask+8));
    
    out_var.add<N> (x0*x0 + x1*x1);
}

// Reads 128 data values (float32)
// Reads 128 mask values (int32)
// Reads 8 mean values (float32)
// Updates (not writes!) 8 var values.

inline void _ds_kernel2b(const float *in_data, const int *in_mask, const float *in_mean, float *update_var)
{
    incremental_upsample8 mean(in_mean);
    incremental_downsample8 var;
    
    _ds_kernel2a<0> (in_data, in_mask, mean, var);
    _ds_kernel2a<1> (in_data+16, in_mask+16, mean, var);
    _ds_kernel2a<2> (in_data+32, in_mask+32, mean, var);
    _ds_kernel2a<3> (in_data+48, in_mask+48, mean, var);
    _ds_kernel2a<4> (in_data+64, in_mask+64, mean, var);
    _ds_kernel2a<5> (in_data+80, in_mask+80, mean, var);
    _ds_kernel2a<6> (in_data+96, in_mask+96, mean, var);
    _ds_kernel2a<7> (in_data+112, in_mask+112, mean, var);
    
    var.update(update_var);
}



// kernel2: computes offsets, scales
//
// The 'in_data' and 'in_mask' arrays have shape (nupfreq,nt/2).  (Not shape (nupfreq,nt)!)
//
// The w0,w1,w2 arrays have length (nt/32).
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
// Assumes nt is divisible by 256, and nt_per_packet=16.
//
// Reminder: packets are decoded and encoded as follows
//   (decoded intensity) = (dec_scale) * (8-bit value) + offset
//   (encoded 8-bit value) = (enc_scale) * (intensity - offset)


inline void _ds_kernel2(const float *in_data, const int *in_mask, float *w0, float *w1, float *w2, int nupfreq, int nt)
{
    memset(w2, 0, (nt/32) * sizeof(float));

    for (int iupfreq = 0; iupfreq < nupfreq; iupfreq++) {
	for (int i = 0; i < (nt/256); i++)
	    _ds_kernel2b(in_data + iupfreq*(nt/2) + 128*i,
			 in_mask + iupfreq*(nt/2) + 128*i,
			 w1 + 8*i,
			 w2 + 8*i);
    }

    __m256 c1 = _mm256_set1_ps(1.0);
    __m256 c004 = _mm256_set1_ps(0.04);
    __m256 c128 = _mm256_set1_ps(128.0);
    __m256 eps = _mm256_set1_ps(1.0e-10);

    for (int i = 0; i < nt/32; i += 8) {
	__m256 count = _mm256_loadu_ps(w0+i);
	__m256 mean = _mm256_loadu_ps(w1+i);
	__m256 var = _mm256_loadu_ps(w2+i) / _mm256_max_ps(count,c1);

	var += eps * mean * mean;

	__m256 rms = _mm256_sqrt_ps(var);
	__m256 mask = _mm256_cmp_ps(rms, _mm256_setzero_ps(), _CMP_GT_OQ);

	rms = _mm256_blendv_ps(c1, rms, mask);

	__m256 dec_scale = c004 * rms;
	__m256 enc_scale = c1 / dec_scale;
	__m256 offset = mean - c128 * dec_scale;

	_mm256_storeu_ps(w0+i, offset);
	_mm256_storeu_ps(w1+i, dec_scale);
	_mm256_storeu_ps(w2+i, enc_scale);
    }
}


// Encode kernel.
// Reads 32 data values (float32)
// Reads 32 mask values (float32)
// Calls enc_off::get<2N>() and enc_off::get<2N+1> to get encode offsets.
// Calls enc_scal::get<2N>() and enc_scal::get<2N+1> to get encode scale.
// Writes 32 output values (uint8)

template<int N>
inline void _ds_kernel3a(uint8_t *out, const float *data, const int *mask, incremental_upsample8 &enc_off, incremental_upsample8 &enc_scal, _extract_32_8 &ex)
{
    __m256 x0 = _mm256_loadu_ps(data);
    __m256 x1 = _mm256_loadu_ps(data+8);
    __m256 x2 = _mm256_loadu_ps(data+16);
    __m256 x3 = _mm256_loadu_ps(data+24);

    __m256 off = enc_off.get<2*N>();
    __m256 scal = enc_scal.get<2*N>();
    
    x0 = (x0 - off) * scal;
    x1 = (x1 - off) * scal;

    off = enc_off.get<2*N+1>();
    scal = enc_scal.get<2*N+1>();

    x2 = (x2 - off) * scal;
    x3 = (x3 - off) * scal;

    const float *fmask = (const float *) mask;

    x0 = _mm256_and_ps(x0, _mm256_loadu_ps(fmask));
    x1 = _mm256_and_ps(x1, _mm256_loadu_ps(fmask+8));
    x2 = _mm256_and_ps(x2, _mm256_loadu_ps(fmask+16));
    x3 = _mm256_and_ps(x3, _mm256_loadu_ps(fmask+24));

    __m256i y = ex.extract(x0, x1, x2, x3);

    _mm256_storeu_si256((__m256i *) out, y);
}


// Encode kernel.
// Reads 128 data values (float32)
// Reads 128 mask values (float32)
// Reads 8 offset values from 'w1' (float32)
// Reads 8 inverse_scale values from 'w2' (float32)
// Writes 128 output values (int8)

inline void _ds_kernel3b(uint8_t *out, const float *data, const int *mask, const float *enc_off, const float *enc_scal, _extract_32_8 &ex)
{
    incremental_upsample8 w1(enc_off);
    incremental_upsample8 w2(enc_scal);

    _ds_kernel3a<0> (out, data, mask, w1, w2, ex);
    _ds_kernel3a<1> (out+32, data+32, mask+32, w1, w2, ex);
    _ds_kernel3a<2> (out+64, data+64, mask+64, w1, w2, ex);
    _ds_kernel3a<3> (out+96, data+96, mask+96, w1, w2, ex);
}


// Fast encode kernel.
// Note that the arrays here have (nt/2) time samples, not nt time samples!
// Assumes nt is a multiple of 256, and that 'out' has stride nt (not stride nt/2).
// Assumes nt_per_packet = 16.
//
// Reads (nupfreq,nt/2) data values (float32)
// Reads (nupfreq,nt/2) mask values (int32)
// Reads (nt/32) enc_offset values (float32)
// Reads (nt/32) enc_scale values (float32)
// Writes (nupfreq,nt/2) data values (uint8_t).

inline void _ds_kernel3(uint8_t *out, const float *data, const int *mask, const float *enc_off, const float *enc_scal, int nupfreq, int nt)
{
    _extract_32_8 ex;

    for (int iupfreq = 0; iupfreq < nupfreq; iupfreq++) {
	for (int i = 0; i < (nt/256); i++) {
	    _ds_kernel3b(out + iupfreq*nt + 128*i, 
			 data + iupfreq*(nt/2) + 128*i, 
			 mask + iupfreq*(nt/2) + 128*i, 
			 enc_off + 8*i,
			 enc_scal + 8*i,
			 ex);
	}
    }
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
    
    // assumed in downsampling kenrnel
    if (constants::nt_per_assembled_chunk % 256 != 0)
	throw runtime_error("ch_frb_io: internal error: fast_assembled_chunk constructor called nt_coarse not divisible by 256");
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


// Virtual override.
// There is some code duplication between this and assembled_chunk::downsample(), so be sure to make changes in sync.
void fast_assembled_chunk::downsample(const assembled_chunk *src1, const assembled_chunk *src2)
{
    int nfreq_c = constants::nfreq_coarse_tot;
    int nt_f = constants::nt_per_assembled_chunk;
    int nt_c = nt_f / nt_per_packet;

    if (!src1 || !src2)
	throw runtime_error("ch_frb_io: null pointer in assembled_chunk::downsample()");
    if (this == src2)
	throw runtime_error("ch_frb_io: assembled_chunk::downsample: 'this' and 'src2' pointers cannot be equal");

    if (!equal3(this->beam_id, src1->beam_id, src2->beam_id))
	throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched beam_id");
    if (!equal3(this->nupfreq, src1->nupfreq, src2->nupfreq))
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched nupfreq");
    if (!equal3(this->nt_per_packet, src1->nt_per_packet, src2->nt_per_packet))
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched nupfreq");

    if (src1->binning != src2->binning)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched binning");
    if (src2->ichunk != src1->ichunk + src2->binning)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: source ichunks are not consecutive");

    for (int ifreq_c = 0; ifreq_c < nfreq_c; ifreq_c++) {
	int ifreq_f = ifreq_c * nupfreq;

	float *out_offsets = this->offsets + (ifreq_c * nt_c);
	float *out_scales = this->scales + (ifreq_c * nt_c);

	_ds_kernel1(ds_data, ds_mask,
		    src1->data + ifreq_f * nt_f,
		    src1->offsets + ifreq_c * nt_c,
		    src1->offsets + ifreq_c * nt_c,
		    out_offsets, out_scales, nupfreq, nt_f);

	_ds_kernel2(ds_data, ds_mask, out_offsets, out_scales, ds_w2, nupfreq, nt_f);

	_ds_kernel3(this->data + (ifreq_f * nt_f),
		    ds_data, ds_mask, out_offsets, ds_w2, nupfreq, nt_f);

	out_offsets = this->offsets + (ifreq_c * nt_c) + (nt_c/2);
	out_scales = this->scales + (ifreq_c * nt_c) + (nt_c/2);

	_ds_kernel1(ds_data, ds_mask,
		    src2->data + ifreq_f * nt_f,
		    src2->offsets + ifreq_c * nt_c,
		    src2->offsets + ifreq_c * nt_c,
		    out_offsets, out_scales, nupfreq, nt_f);

	_ds_kernel2(ds_data, ds_mask, out_offsets, out_scales, ds_w2, nupfreq, nt_f);

	_ds_kernel3(this->data + (ifreq_f * nt_f) + (nt_f/2),
		    ds_data, ds_mask, out_offsets, ds_w2, nupfreq, nt_f);
    }

    this->binning = src1->binning + 1;
    this->ichunk = src1->ichunk;
    this->isample = src1->isample;

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


static void test_avx2_ds_kernel1(std::mt19937 &rng, int nupfreq, int nt_f)
{
    ch_assert(nt_f % 32 == 0);

    int nt_per_packet = 16;
    int nt_c = nt_f / nt_per_packet;

    vector<uint8_t> in_data(nupfreq * nt_f, 0);

    // Randomly simulate in_data.
    // We assign ~10% probability to "masked" values 0x00 or 0xff
    
    for (int i = 0; i < nupfreq * nt_f; i++) {
	int x = randint(rng, -32, 256+32);
	x = max(x, 0);
	x = min(x, 255);
	in_data[i] = x;
    }

    // To test corner cases, assign ~10% probability to masking an entire (nupfreq,32) block.
    for (int it_c = 0; it_c < (nt_f/32); it_c++) {
	if (uniform_rand(rng) > 0.1)
	    continue;
	
	for (int ifreq = 0; ifreq < nupfreq; ifreq++)
	    for (int it = 32*it_c; it < 32*(it_c+1); it++)
		in_data[ifreq*nt_f + it] = (uniform_rand(rng) < 0.5) ? 0x00 : 0xff;
    }

    vector<float> in_scales = uniform_randvec<float> (rng, nt_c, 1.0e-3, 1.0e-2);
    vector<float> in_offsets = uniform_randvec<float> (rng, nt_c, -1.0, 1.0);

    vector<int> out1_mask = randintvec(rng, nupfreq * (nt_f/2), -10, 10);
    vector<float> out1_data = uniform_randvec<float> (rng, nupfreq * (nt_f/2), -1.0, 1.0);
    vector<float> out1_count = uniform_randvec<float> (rng, nupfreq * (nt_f/2), -100.0, 100.0);
    vector<float> out1_mean = uniform_randvec<float> (rng, nupfreq * (nt_f/2), -100.0, 100.0);

    vector<int> out2_mask = randintvec(rng, nupfreq * (nt_f/2), -10, 10);
    vector<float> out2_data = uniform_randvec<float> (rng, nupfreq * (nt_f/2), -1.0, 1.0);
    vector<float> out2_count = uniform_randvec<float> (rng, nupfreq * (nt_f/2), -100.0, 100.0);
    vector<float> out2_mean = uniform_randvec<float> (rng, nupfreq * (nt_f/2), -100.0, 100.0);

    // Slow kernel.
    ds_slow_kernel1(&out1_data[0], &out1_mask[0], &in_data[0], &in_offsets[0], &in_scales[0], &out1_count[0], &out1_mean[0], nupfreq, nt_f, nt_per_packet);
    
    // Fast kernel.
    _ds_kernel1(&out2_data[0], &out2_mask[0], &in_data[0], &in_offsets[0], &in_scales[0], &out2_count[0], &out2_mean[0], nupfreq, nt_f);

    for (int i = 0; i < nupfreq; i++) {
	for (int j = 0; j < nt_f/2; j++) {
	    unsigned int d0 = in_data[i*nt_f + 2*j];
	    unsigned int d1 = in_data[i*nt_f + 2*j+1];
	    int m_slow = out1_mask[i*(nt_f/2) + j];
	    int m_fast = out2_mask[i*(nt_f/2) + j];

	    float d_slow = out1_data[i*(nt_f/2) + j];
	    float d_fast = out2_data[i*(nt_f/2) + j];

	    if ((m_slow != m_fast) || (fabs(d_slow - d_fast) > 1.0e-3)) {
		cerr << "test_avx2_ds_kernel1 failed: (nupfreq,nt_f)=(" << nupfreq << "," << nt_f << "),"
		     << " (i,j)=(" << i << "," << j << "),"
		     << " (d0,d1)=(" << hex << d0 << "," << d1 << "), " << dec
		     << " (m_slow,m_fast)=(" << m_slow << "," << m_fast << "),"
		     << " (d_slow,d_fast)=(" << d_slow << "," << d_fast << ")\n";
		throw runtime_error("test_avx2_ds_kernel1() failed");
	    }
	}
    }

    for (int i = 0; i < (nt_c/2); i++) {
	if ((fabs(out1_count[i] - out2_count[i]) > 1.0e-3) || (fabs(out1_mean[i] - out2_mean[i]) > 1.0e-3)) {
	    cerr << "test_avx2_ds_kernel1 failed: (nupfreq,nt_f)=(" << nupfreq << "," << nt_f << "), i=" << i
		 << " (count_slow,count_fast)=(" << out1_count[i] << "," << out2_count[i] << "),"
		 << " (mean_slow,mean_fast)=(" << out1_mean[i] << "," << out2_mean[i] << ")\n";
	    throw runtime_error("test_avx2_ds_kernel1() failed");
	}
    }
}


static void test_avx2_ds_kernel1(std::mt19937 &rng)
{
    cout << "test_avx2_ds_kernel1..";

    for (int iouter = 0; iouter < 1000; iouter++) {
	if (iouter % 10 == 0)
	    cout << "." << flush;

	int nupfreq = randint(rng, 1, 31);
	int nt_f = 256 * randint(rng, 1, 9);
	test_avx2_ds_kernel1(rng, nupfreq, nt_f);
    }

    cout << "success" << endl;
}


void test_avx2_kernels(std::mt19937 &rng)
{
    test_avx2_ds_kernel1(rng);

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
