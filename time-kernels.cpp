// TODO: add encode kernels

#include "ch_frb_io_internals.hpp"

using namespace std;
using namespace ch_frb_io;


static void time_decode(std::mt19937 &rng)
{
    // Note: 'nchunks' is a multiplier, intended to prevent the whole instance from fitting in L3 cache.
    const int beam_id = 0;
    const int nupfreq = 16;
    const int nt_per_packet = 16;
    const int fpga_counts_per_sample = 384;
    const int stride = 16384;
    const int nchunks = 2;
    const int niter = 5;

    const int nfreq_c = constants::nfreq_coarse_tot;
    const int nfreq_f = nfreq_c * nupfreq;
    const int nt_f = constants::nt_per_assembled_chunk;
    const int nt_c = nt_f / nt_per_packet;

    vector<float> intensity(nchunks * nfreq_f * stride, 0.0);
    vector<float> weights(nchunks * nfreq_f * stride, 0.0);

    vector<shared_ptr<assembled_chunk>> chunks(nchunks);

    //
    // Use fast_assembled_chunks
    //

    for (int ichunk = 0; ichunk < nchunks; ichunk++) {
	chunks[ichunk] = make_shared<fast_assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);
	chunks[ichunk]->randomize(rng);   // I don't think this matters
    }
	
    struct timeval tv0 = xgettimeofday();

    for (int iter = 0; iter < niter; iter++)
	for (int ichunk = 0; ichunk < nchunks; ichunk++)
	    chunks[ichunk]->decode(&intensity[ichunk * nfreq_f * stride],
				   &weights[ichunk * nfreq_f * stride],
				   stride);

    struct timeval tv1 = xgettimeofday();

    double cpu_time = 1.0e-6 * usec_between(tv0, tv1);
    double real_time = double(niter) * double(nchunks) * double(nt_f) * double(fpga_counts_per_sample) * constants::dt_fpga;

    cout << "decode: loadfrac/beam = " << (cpu_time / real_time) << endl;

    //
    // Use (slow) assembled_chunks
    //

    for (int ichunk = 0; ichunk < nchunks; ichunk++) {
	chunks[ichunk] = make_shared<assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);
	chunks[ichunk]->randomize(rng);   // I don't think this matters
    }
	
    tv0 = xgettimeofday();

    for (int iter = 0; iter < niter; iter++)
	for (int ichunk = 0; ichunk < nchunks; ichunk++)
	    chunks[ichunk]->decode(&intensity[ichunk * nfreq_f * stride],
				   &weights[ichunk * nfreq_f * stride],
				   stride);

    tv1 = xgettimeofday();

    cpu_time = 1.0e-6 * usec_between(tv0, tv1);
    cout << "slow decode: loadfrac/beam = " << (cpu_time / real_time) << endl;    
}


int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 rng(rd());

    time_decode(rng);
    return 0;
}
