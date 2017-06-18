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


static void time_downsample(std::mt19937 &rng)
{
    const int nt_f = constants::nt_per_assembled_chunk;

    // Note: 'nchunks' is a multiplier, intended to prevent the whole instance from fitting in L3 cache.
    const int beam_id = 0;
    const int nupfreq = 16;
    const int nt_per_packet = 16;
    const int fpga_counts_per_sample = 384;
    const int nchunks = 2;
    const int niter = 5;

    vector<shared_ptr<assembled_chunk>> dst_chunks(nchunks);
    vector<shared_ptr<assembled_chunk>> src_chunks(2 * nchunks);

    for (int ichunk = 0; ichunk < 2*nchunks; ichunk++) { 
	src_chunks[ichunk] = make_shared<fast_assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);
	src_chunks[ichunk]->randomize(rng);
    }

    //
    // Use fast_assembled_chunks
    //
    
    for (int ichunk = 0; ichunk < nchunks; ichunk++)
	dst_chunks[ichunk] = make_shared<fast_assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);

    struct timeval tv0 = xgettimeofday();

    for (int iter = 0; iter < niter; iter++)
	for (int ichunk = 0; ichunk < nchunks; ichunk++)
	    dst_chunks[ichunk]->downsample(src_chunks[2*ichunk], src_chunks[2*ichunk+1]);

    struct timeval tv1 = xgettimeofday();

    double cpu_time = 15.0 * 1.0e-6 * usec_between(tv0, tv1);   // note factor 15 here!
    double real_time = double(niter) * double(nchunks) * double(nt_f) * double(fpga_counts_per_sample) * constants::dt_fpga;

    cout << "downsample: loadfrac / (8 beams) = " << (cpu_time / real_time) << endl;

    // 
    // Use (slow) assembled_chunks
    //
    
    for (int ichunk = 0; ichunk < nchunks; ichunk++)
	dst_chunks[ichunk] = make_shared<assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);

    tv0 = xgettimeofday();

    for (int iter = 0; iter < niter; iter++)
	for (int ichunk = 0; ichunk < nchunks; ichunk++)
	    dst_chunks[ichunk]->downsample(src_chunks[2*ichunk], src_chunks[2*ichunk+1]);

    tv1 = xgettimeofday();

    cpu_time = 15.0 * 1.0e-6 * usec_between(tv0, tv1);   // note factor 15 here!
    cout << "slow downsample: loadfrac / (8 beams) = " << (cpu_time / real_time) << endl;
}


int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 rng(rd());

    time_decode(rng);
    time_downsample(rng);

    return 0;
}
