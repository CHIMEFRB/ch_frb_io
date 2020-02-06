#include <cassert>
#include <iostream>
#include "ch_frb_io_internals.hpp"
#include "chlog.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif



// Initializes a 'struct intensity_packet' from raw packet data.  The "pointer" fields of the
// struct intensity_packet are initialized to pointers into the 'src' buffer, so the caller is
// responsible for ensuring that this buffer doesn't get freed while the struct intensity_packet 
// is in scope.
//
// Does a bunch of sanity checks and returns 'true' if packet is good, 'false' if bad.
//
// Explicitly, the following checks are performed:
//   - protocol version == 1
//   - dimensions (nbeams, nfreq_coarse, nupfreq, ntsamp) are not large enough to lead to integer overflows
//   - packet and data byte counts are correct
//   - coarse_freq_ids are in range
//   - ntsamp is a power of two, and in the range (0,max_allowed_nt_per_packet].
//   - nbeams, nfreq_coarse, nupfreq, ntsamp, fpga_counts_per_sample are all > 0
//   - fpga_count is a multiple of (fpga_counts_per_sample * ntsamp)

bool intensity_packet::decode(const uint8_t *src, int src_nbytes)
{
    if (_unlikely(src_nbytes < 24)) {
        chlog("packet nbytes < 24");
	return false;
    }
    if (_unlikely(src_nbytes > constants::max_input_udp_packet_size)) {
        chlog("packet nbytes > " << constants::max_input_udp_packet_size);
	return false;
    }

    memcpy(this, src, 24);

    if (_unlikely(protocol_version != 1)) {
        chlog("packet protocol version bad");
	return false;
    }
    if (_unlikely((ntsamp > constants::max_allowed_nt_per_packet) || ((ntsamp & (ntsamp-1)) != 0))) {
        chlog("packet ntsamp bad");
	return false;
    }
    if (_unlikely(fpga_counts_per_sample == 0)) {
        chlog("packet fpga_counts_per_sample bad");
	return false;
    }

    // Note conversions to uint64_t, to prevent integer overflow
    uint64_t fpga_counts_per_packet = uint64_t(fpga_counts_per_sample) * uint64_t(ntsamp);
    if (_unlikely(fpga_count % fpga_counts_per_packet != 0)) {
        chlog("packet fpga_count bad");
	return false;
    }
	
    uint64_t n1 = uint64_t(nbeams);
    uint64_t n2 = uint64_t(nfreq_coarse);
    uint64_t n3 = uint64_t(nupfreq);
    uint64_t n4 = uint64_t(ntsamp);

    // Expected header, data size
    uint64_t nh = 24 + 2*n1 + 2*n2 + 8*n1*n2;
    uint64_t nd = n1 * n2 * n3 * n4;

    if (_unlikely(uint64_t(src_nbytes) != nh+nd)) {
        chlog("packet src_nbytes bad");
	return false;
    }
    if (_unlikely(uint64_t(data_nbytes) != nd)) {
        chlog("packet data_nbytes bad");
	return false;
    }

    this->beam_ids = (uint16_t *) (src + 24);
    this->coarse_freq_ids = (uint16_t *) (src + 24 + 2*n1);
    this->scales = (float *) (src + 24 + 2*n1 + 2*n2);
    this->offsets = (float *) (src + 24 + 2*n1 + 2*n2 + 4*n1*n2);
    this->data = (uint8_t *) (src + nh);

    for (int i = 0; i < nfreq_coarse; i++)
	if (_unlikely(coarse_freq_ids[i] >= constants::nfreq_coarse_tot)) {
            chlog("packet coarse_freq_ids bad");
	    return false;
        }

    return true;
}

int intensity_packet::set_pointers(uint8_t *dst) {
    int nb = this->nbeams;
    int nf = this->nfreq_coarse;
    int nu = this->nupfreq;
    int nt = this->ntsamp;

    memcpy(dst, this, 24);

    this->beam_ids = (uint16_t *)(dst + 24);
    this->coarse_freq_ids = (uint16_t *)(dst + 24 + 2*nb);
    this->scales  = (float *) (dst + 24 + 2*nb + 2*nf);
    this->offsets = (float *) (dst + 24 + 2*nb + 2*nf + 4*nb*nf);
    this->data = dst + 24 + 2*nb + 2*nf + 8*nb*nf;

    return 24 + 2*nb + 2*nf + 8*nb*nf + nb*nf*nu*nt;
}

// Encodes a floating-point array of intensities into raw packet data, before sending packet.
// The precise semantics aren't very intuitive, see extended comment in ch_frb_io_internals.hpp for details!

int intensity_packet::encode(uint8_t *dst, const float *intensity, int beam_istride, int freq_istride, const float *weights, int beam_wstride, int freq_wstride, float wt_cutoff)
{
    int nb = this->nbeams;
    int nf = this->nfreq_coarse;
    int nu = this->nupfreq;
    int nt = this->ntsamp;

    // similar but not exactly the same as set_pointers()...
    memcpy(dst, this, 24);
    memcpy(dst + 24, this->beam_ids, 2*nb);
    memcpy(dst + 24 + 2*nb, this->coarse_freq_ids, 2*nf);

    this->scales  = (float *) (dst + 24 + 2*nb + 2*nf);
    this->offsets = (float *) (dst + 24 + 2*nb + 2*nf + 4*nb*nf);
    this->data = dst + 24 + 2*nb + 2*nf + 8*nb*nf;

    int nbytes = 24 + 2*nb + 2*nf + 8*nb*nf + nb*nf*nu*nt;

    for (int b = 0; b < nb; b++) {
	for (int f = 0; f < nf; f++) {
	    uint8_t *sub_data = data + (b*nf+f) * (nu*nt);
	    const float *sub_int = intensity + b*beam_istride + f*nu*freq_istride;
	    const float *sub_wt  = weights   + b*beam_wstride + f*nu*freq_wstride;

	    float acc0 = 0.0;
	    float acc1 = 0.0;
	    float acc2 = 0.0;

	    for (int u = 0; u < nu; u++) {
		for (int t = 0; t < nt; t++) {
		    float x = sub_int[u*freq_istride+t];
		    float w = (sub_wt[u*freq_wstride+t] >= wt_cutoff) ? 1.0 : 0.0;

		    acc0 += w;
		    acc1 += w * x;
		    acc2 += w * x * x;
		}
	    }
		    
	    if (acc0 <= 0.0) {
		this->scales [b*nf+f] = 1.0;
		this->offsets[b*nf+f] = 0.0;
		memset(sub_data, 0, nu*nt);
		continue;
	    }

	    float mean = acc1/acc0;
	    float var = acc2/acc0 - mean*mean;
    
	    // Since we use single precision, the numerical error in 'var' is approx (1.0e-7 * mean*mean),
	    // so we need to regulate values of 'var' which are of this order or smaller.  The radiometer
	    // equation predicts that the ideal variance is (1.0e-3 * mean*mean) or more, depending on
	    // the correlator configuration.  We choose a regulator which is safely in between these values.

	    var = max(var, float(1.0e-5) * mean*mean);

	    float scale = sqrt(var) / 25.;
	    float offset = -128.*scale + mean;   // 0x80 -> mean

	    this->scales [b*nf+f] = scale;
	    this->offsets[b*nf+f] = offset;

	    for (int u = 0; u < nu; u++) {
		for (int t = 0; t < nt; t++) {
		    float x = sub_int[u*freq_istride+t];
		    float w = (sub_wt[u*freq_wstride+t] >= wt_cutoff) ? 1.0 : 0.0;

		    x = w * (x - offset) / scale;
		    x = min(x, float(255.));
		    x = max(x, float(0.));
		    sub_data[u*nt+t] = uint8_t(x+0.5);   // round to nearest integer
		}
	    }
	}
    }
    return nbytes;
}


int intensity_packet::find_coarse_freq_id(int id) const
{
    for (int i = 0; i < this->nfreq_coarse; i++)
	if (this->coarse_freq_ids[i] == id)
	    return i;
    return -1;
}


bool intensity_packet::contains_coarse_freq_id(int id) const
{
    int i = this->find_coarse_freq_id(id);
    return (i >= 0);
}


// This test is kinda silly, but checks that the byte alignment of the intensity_packet header fields
// is what I think it is.  (Just worried that the compiler might insert some padding bytes.)
void test_packet_offsets(std::mt19937 &rng)
{
    cerr << "test_packet_offsets()...";
    vector<uint8_t> buf(24, 0);

    for (int iouter = 0; iouter < 1000; iouter++) {
	intensity_packet p;

	// Randomized header fields
	uint32_t protocol_version = std::uniform_int_distribution<uint32_t>()(rng);
	int16_t data_nbytes = std::uniform_int_distribution<int16_t>()(rng);
	uint16_t fpga_counts_per_sample = std::uniform_int_distribution<uint16_t>()(rng);
	uint64_t fpga_count = std::uniform_int_distribution<uint64_t>()(rng);
	uint16_t nbeams = std::uniform_int_distribution<uint16_t>()(rng);
	uint16_t nfreq_coarse = std::uniform_int_distribution<uint16_t>()(rng);
	uint16_t nupfreq = std::uniform_int_distribution<uint16_t>()(rng);
	uint16_t ntsamp = std::uniform_int_distribution<uint16_t>()(rng);

	*((uint32_t *) &buf[0]) = protocol_version;
	*((int16_t *) &buf[4]) = data_nbytes;
	*((uint16_t *) &buf[6]) = fpga_counts_per_sample;
	*((uint64_t *) &buf[8]) = fpga_count;
	*((uint16_t *) &buf[16]) = nbeams;
	*((uint16_t *) &buf[18]) = nfreq_coarse;
	*((uint16_t *) &buf[20]) = nupfreq;
	*((uint16_t *) &buf[22]) = ntsamp;

	memcpy(&p, &buf[0], 24);

	assert(p.protocol_version == protocol_version);
	assert(p.data_nbytes == data_nbytes);
	assert(p.fpga_counts_per_sample == fpga_counts_per_sample);
	assert(p.fpga_count == fpga_count);
	assert(p.nbeams == nbeams);
	assert(p.nfreq_coarse == nfreq_coarse);
	assert(p.nupfreq == nupfreq);
	assert(p.ntsamp == ntsamp);
    }

    cerr << "success\n";
}


}   // namespace ch_frb_io
