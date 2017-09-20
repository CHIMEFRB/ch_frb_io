#include <cassert>
#include <algorithm>
#include "ch_frb_io_internals.hpp"

using namespace std;
using namespace ch_frb_io;


inline vector<int> vrange(int n)
{
    vector<int> ret(n);
    for (int i = 0; i < n; i++)
	ret[i] = i;
    return ret;
}


static void test_encode_decode(std::mt19937 &rng)
{
    cerr << "test_encode_decode()";

    for (int iouter = 0; iouter < 100; iouter++) {
	cerr << ".";

	// First, we make a random test instance

#ifdef __AVX2__
	bool use_fast_kernels = (iouter % 2) == 1;    // alternating rather than random
#else
	bool use_fast_kernels = false;
#endif

	int nbeams = randint(rng,1,9);
	int nfreq_coarse_per_packet = (1 << randint(rng,0,3));
	int nupfreq = use_fast_kernels ? (2 * randint(rng,1,9)) : randint(rng,1,17);
	int nt_per_packet = use_fast_kernels ? 16 : (1 << randint(rng,0,5));
	int nfreq_coarse_tot = constants::nfreq_coarse_tot;
	int nt_per_chunk = constants::nt_per_assembled_chunk;
	double wt_cutoff = uniform_rand(rng, 0.2, 0.3);
	int src_istride = randint(rng, nt_per_chunk, 2*nt_per_chunk);
	int src_wstride = randint(rng, nt_per_chunk, 2*nt_per_chunk);
	int dst_istride = randint(rng, nt_per_chunk, 2*nt_per_chunk);
	int dst_wstride = randint(rng, nt_per_chunk, 2*nt_per_chunk);
	int ichunk = randint(rng, 0, 1024);
	int fpga_counts_per_sample = randint(rng, 1, 1024);

	vector<int> send_freq_ids = vrange(constants::nfreq_coarse_tot);
	std::shuffle(send_freq_ids.begin(), send_freq_ids.end(), rng);

	vector<float> src_intensity(nbeams * nfreq_coarse_tot * nupfreq * src_istride, 0.0);
	vector<float> src_weights(nbeams * nfreq_coarse_tot * nupfreq * src_wstride, 0.0);

	for (int i = 0; i < nbeams * nfreq_coarse_tot * nupfreq; i++) {
	    for (int it = 0; it < nt_per_chunk; it++) {
		src_intensity[i*src_istride + it] = uniform_rand(rng);
		src_weights[i*src_wstride + it] = uniform_rand(rng);
	    }
	}

	intensity_network_ostream::initializer ini_params;
	ini_params.dstname = "127.0.0.1";
	ini_params.beam_ids = vrange(nbeams);
	ini_params.coarse_freq_ids = send_freq_ids;
	ini_params.nupfreq = nupfreq;
	ini_params.nt_per_chunk = nt_per_chunk;
	ini_params.nfreq_coarse_per_packet = nfreq_coarse_per_packet;
	ini_params.nt_per_packet = nt_per_packet;
	ini_params.fpga_counts_per_sample = fpga_counts_per_sample;
	ini_params.wt_cutoff = wt_cutoff;
	ini_params.target_gbps = 1.0;  // dummy

	// Note: intensity_network_ostream::make() initializes some things we don't need 
	// (e.g. spawns network thread) but that's OK.

	auto ostream = intensity_network_ostream::make(ini_params);
	auto packet_list = make_unique<udp_packet_list> (ostream->npackets_per_chunk, ostream->nbytes_per_chunk);
	vector<shared_ptr<assembled_chunk>> assembled_chunks(nbeams);

	for (int ibeam = 0; ibeam < nbeams; ibeam++) {
	    assembled_chunk::initializer ini_params;
	    ini_params.beam_id = ibeam;
	    ini_params.nupfreq = nupfreq;
	    ini_params.nt_per_packet = nt_per_packet;
	    ini_params.fpga_counts_per_sample = fpga_counts_per_sample;
	    ini_params.ichunk = ichunk;
	    ini_params.force_reference = !use_fast_kernels;
	    ini_params.force_fast = use_fast_kernels;

	    assembled_chunks[ibeam] = assembled_chunk::make(ini_params);
	}

	ostream->_encode_chunk(&src_intensity[0], src_istride, &src_weights[0], src_wstride, ichunk * nt_per_chunk * fpga_counts_per_sample, packet_list);
	
	for (int ipacket = 0; ipacket < packet_list->curr_npackets; ipacket++) {
	    uint8_t *packet_data = packet_list->get_packet_data(ipacket);
            int packet_nbytes = packet_list->get_packet_nbytes(ipacket);

	    intensity_packet packet;
	    if (!packet.decode(packet_data, packet_nbytes))
		throw runtime_error("intensity_packet::decode() failed");

	    if (packet.nbeams != nbeams)
		throw runtime_error("nbeams mismatch");
	    if (packet.nfreq_coarse != nfreq_coarse_per_packet)
		throw runtime_error("nfreq_coarse_per_packet mismatch");
	    if (packet.nupfreq != nupfreq)
		throw runtime_error("nupfreq mismatch");
	    if (packet.ntsamp != nt_per_packet)
		throw runtime_error("nt_per_packet mismatch");
	    if (packet.fpga_counts_per_sample != fpga_counts_per_sample)
		throw runtime_error("fpga_counts_per_sample mismatch");

	    for (int ibeam = 0; ibeam < nbeams; ibeam++)
		if (packet.beam_ids[ibeam] != ibeam)
		    throw runtime_error("beam_id mismatch");

	    // Danger zone
	    int new_data_nbytes = nfreq_coarse_per_packet * nupfreq * nt_per_packet;
	    packet.data_nbytes = new_data_nbytes;
	    packet.nbeams = 1;
	    
	    for (int ibeam = 0; ibeam < nbeams; ibeam++) {
		assembled_chunks[ibeam]->add_packet(packet);
		
		// Danger zone
		packet.beam_ids += 1;
		packet.scales += nfreq_coarse_per_packet;
		packet.offsets += nfreq_coarse_per_packet;
		packet.data += new_data_nbytes;
	    }
	}

	vector<float> dst_intensity(nfreq_coarse_tot * nupfreq * dst_istride, 1.0e30);
	vector<float> dst_weights(nfreq_coarse_tot * nupfreq * dst_wstride, 1.0e30);

	for (int ibeam = 0; ibeam < nbeams; ibeam++) {
	    assembled_chunks[ibeam]->decode(&dst_intensity[0], &dst_weights[0], dst_istride, dst_wstride);

	    for (int ifreq_coarse = 0; ifreq_coarse < nfreq_coarse_tot; ifreq_coarse++) {
		for (int iupfreq = 0; iupfreq < nupfreq; iupfreq++) {
		    int si = ((ibeam*nfreq_coarse_tot + ifreq_coarse) * nupfreq + iupfreq) * src_istride;
		    int sw = ((ibeam*nfreq_coarse_tot + ifreq_coarse) * nupfreq + iupfreq) * src_wstride;
		    int di = ((send_freq_ids[ifreq_coarse] * nupfreq) + iupfreq) * dst_istride;
		    int dw = ((send_freq_ids[ifreq_coarse] * nupfreq) + iupfreq) * dst_wstride;

		    for (int it = 0; it < nt_per_chunk; it++) {
			if (dst_weights[dw+it] == 0.0)
			    assert(src_weights[sw+it] <= 1.00001 * wt_cutoff);
			else if (dst_weights[dw+it] == 1.0) {
			    assert(src_weights[sw+it] >= 0.99999 * wt_cutoff);
			    assert(fabs(dst_intensity[di+it] - src_intensity[si+it]) < 0.02);
			}
			else
			    throw runtime_error("dst_weights not equal to 0 or 1?!");
		    }
		}
	    }
	}
    }

    cerr << "success\n";
}


// -------------------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 rng(rd());

    test_lexical_cast();       // defined in lexical_cast.cpp
    test_packet_offsets(rng);  // defined in intensity_packet.cpp
    test_avx2_kernels(rng);    // defined in avx2_kernels.cpp
    test_encode_decode(rng);   // defined above

    return 0;
}
