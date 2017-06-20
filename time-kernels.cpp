// TODO: add encode kernels

#include "ch_frb_io_internals.hpp"

using namespace std;
using namespace ch_frb_io;


template<typename Tchunk>
static double time_decode(std::mt19937 &rng)
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

    for (int ichunk = 0; ichunk < nchunks; ichunk++) {
	chunks[ichunk] = make_shared<Tchunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);
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

    return cpu_time / real_time;
}


static void time_assemble(std::mt19937 &rng)
{
    // Note: 'nchunks' is a multiplier, intended to prevent the whole instance from fitting in L3 cache.
    const int beam_id = 0;
    const int nupfreq = 16;
    const int nt_per_packet = 16;
    const int fpga_counts_per_sample = 384;
    const int nchunks = 2;
    const int niter = 5;

    const int nfreq_c = constants::nfreq_coarse_tot;
    const int nt_f = constants::nt_per_assembled_chunk;
    const int nt_c = nt_f / nt_per_packet;

    // List of (ifreq_c, it_c) pairs.
    vector<pair<int,int>> packet_ordering(nfreq_c/4 * nt_c);

    for (int ifreq_c = 0; ifreq_c < nfreq_c/4; ifreq_c++)
	for (int it_c = 0; it_c < nt_c; it_c++)
	    packet_ordering[ifreq_c*nt_c + it_c] = { ifreq_c, it_c };

    ssize_t npackets_per_chunk = (nfreq_c/4) * nt_c;
    ssize_t packet_size = intensity_packet::packet_size(1, 4, nupfreq, nt_per_packet);

    vector<intensity_packet> packets(nchunks * npackets_per_chunk);
    vector<uint8_t> packet_data(nchunks * npackets_per_chunk * packet_size, 0);

    for (int ichunk = 0; ichunk < nchunks; ichunk++) {
	std::shuffle(packet_ordering.begin(), packet_ordering.end(), rng);

	for (int ipacket = 0; ipacket < npackets_per_chunk; ipacket++) {
	    int ifreq_c = packet_ordering[ipacket].first;
	    int it_c = packet_ordering[ipacket].second;

	    intensity_packet &p = packets[ichunk*npackets_per_chunk + ipacket];
	    uint8_t *d = &packet_data[(ichunk*npackets_per_chunk + ipacket) * packet_size];

	    p.protocol_version = 1;
	    p.data_nbytes = 4 * nupfreq * nt_per_packet;
	    p.fpga_counts_per_sample = fpga_counts_per_sample;
	    p.fpga_count = it_c * nt_per_packet * fpga_counts_per_sample;
	    p.nbeams = 1;
	    p.nfreq_coarse = 4;
	    p.nupfreq = nupfreq;
	    p.ntsamp = nt_per_packet;

	    p.beam_ids = reinterpret_cast<uint16_t *> (d + 24);
	    p.coarse_freq_ids = reinterpret_cast<uint16_t *> (d + 26);
	    p.scales = reinterpret_cast<float *> (d + 34);
	    p.offsets = reinterpret_cast<float *> (d + 50);
	    p.data = d + 66;

	    p.beam_ids[0] = 0;

	    for (int f = 0; f < 4; f++)
		p.coarse_freq_ids[f] = ifreq_c + f * (nfreq_c/4);
	}
    }

    vector<shared_ptr<assembled_chunk>> chunks(nchunks);

    //
    // Use fast_assembled_chunks
    //

    for (int ichunk = 0; ichunk < nchunks; ichunk++)
	chunks[ichunk] = make_shared<fast_assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, 0);

    struct timeval tv0 = xgettimeofday();

    for (int iter = 0; iter < niter; iter++)
	for (int ichunk = 0; ichunk < nchunks; ichunk++)
	    for (int ipacket = 0; ipacket < npackets_per_chunk; ipacket++)
		chunks[ichunk]->add_packet(packets[ichunk*npackets_per_chunk + ipacket]);

    struct timeval tv1 = xgettimeofday();
    
    double cpu_time = 8.0e-6 * usec_between(tv0, tv1);    // note factor 8 here
    double real_time = double(niter) * double(nchunks) * double(nt_f) * double(fpga_counts_per_sample) * constants::dt_fpga;

    cout << "assemble: loadfrac / (8 beams) = " << (cpu_time / real_time) << endl;

    //
    // Use (slow) assembled_chunks
    //

    for (int ichunk = 0; ichunk < nchunks; ichunk++)
	chunks[ichunk] = make_shared<assembled_chunk> (beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, 0);

    tv0 = xgettimeofday();

    for (int iter = 0; iter < niter; iter++)
	for (int ichunk = 0; ichunk < nchunks; ichunk++)
	    for (int ipacket = 0; ipacket < npackets_per_chunk; ipacket++)
		chunks[ichunk]->add_packet(packets[ichunk*npackets_per_chunk + ipacket]);

    tv1 = xgettimeofday();
    
    cpu_time = 8.0e-6 * usec_between(tv0, tv1);    // note factor 8 here
    cout << "slow assemble: loadfrac / (8 beams) = " << (cpu_time / real_time) << endl;
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

    cout << "slow decode: loadfrac/beam = " << time_decode<assembled_chunk> (rng) << "\n";
    cout << "fast decode: loadfrac/beam = " << time_decode<fast_assembled_chunk> (rng) << "\n";

    time_assemble(rng);
    time_downsample(rng);

    return 0;
}
