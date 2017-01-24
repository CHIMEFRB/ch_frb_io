#include <cassert>
#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;
using namespace ch_frb_io;


inline double intval(int ifreq, int ipol, int it)
{
    return 1.3183*ifreq + 1502.382*ipol - 0.328*it - 3201.91;
}

inline double wtval(int ifreq, int ipol, int it)
{
    // all signs positive
    return 1.021*ifreq + 1502.382*ipol + 0.283*it + 2038.18;
}


// -------------------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
    if (argc == 1) {
        cout << "Warning: this will create a ~100 MB file in the current directory!\n"
             << "If this is OK press return.  If not, press control-C!\n"
             << "I AM WAITING, HUMAN: "
             << flush;
        string dummy;
        getline(cin, dummy);
    }

    const char *filename = "test_intensity_hdf5_file.hdf5";
    const int nfreq = 1024;
    const int npol = 2;
    const int nt_file = 16384;
    const double freq0_MHz = 800.0;
    const double freq1_MHz = 400.0;
    const double dt_sample = 1.0e-3;
    const double initial_time = 1.328;

    std::random_device rd;
    std::mt19937 rng(rd());

    vector<string> pol = { "XX", "YY" };
    unique_ptr<intensity_hdf5_ofile> f = make_unique<intensity_hdf5_ofile> (filename, nfreq, pol, freq0_MHz, freq1_MHz, dt_sample);

    ssize_t curr_file_ipos = 137;
    vector<ssize_t> chunk_ipos_list;
    vector<ssize_t> chunk_nt_list;

    while (f->curr_nt < nt_file) {
	ssize_t chunk_ipos = curr_file_ipos + max(0, randint(rng,-100,100));
	ssize_t chunk_nt = min((ssize_t)randint(rng,1,33), nt_file - f->curr_nt);
	
	vector<float> intensity(nfreq * npol * chunk_nt);
	vector<float> weights(nfreq * npol * chunk_nt);

	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    for (int ipol = 0; ipol < npol; ipol++) {
		for (int it = 0; it < chunk_nt; it++) {
		    int i = ifreq*npol*chunk_nt + ipol*chunk_nt + it;
		    intensity[i] = intval(ifreq, ipol, chunk_ipos + it);
		    weights[i] = wtval(ifreq, ipol, chunk_ipos + it);
		}
	    }
	}

	double chunk_t0 = initial_time + dt_sample * chunk_ipos;
	f->append_chunk(chunk_nt, &intensity[0], &weights[0], chunk_ipos, chunk_t0);

	chunk_ipos_list.push_back(chunk_ipos);
	chunk_nt_list.push_back(chunk_nt);
	curr_file_ipos = chunk_ipos + chunk_nt;
    }

    // calls intensity_hdf5_ofile destructor, which flushes data to disk
    f = unique_ptr<intensity_hdf5_ofile> ();

    // check that data has really been flushed, by immediate renaming file and reading it back
    const char *filename2 = "test_intensity_hdf5_file_renamed.hdf5";

    int err = rename(filename, filename2);
    if (err)
	throw runtime_error("rename() failed");

    cerr << "renamed " << filename << " -> " << filename2 << endl;

    intensity_hdf5_file f2(filename2);
    double expected_time_lo = initial_time + dt_sample * chunk_ipos_list[0];
    double expected_time_hi = initial_time + dt_sample * curr_file_ipos;

    assert(f2.nfreq == nfreq);
    assert(f2.npol == npol);
    assert(f2.nt_file == nt_file);
    assert(!f2.frequencies_are_increasing);
    assert(fabs(f2.freq_lo_MHz - freq1_MHz) < 1.0e-4);
    assert(fabs(f2.freq_hi_MHz - freq0_MHz) < 1.0e-4);
    assert(fabs(f2.dt_sample - dt_sample) < 1.0e-10);
    assert(fabs(f2.time_lo - expected_time_lo) < 1.0e-7);
    assert(fabs(f2.time_hi - expected_time_hi) < 1.0e-7);

    ssize_t curr_it_file = 0;

    for (unsigned int ichunk = 0; ichunk < chunk_ipos_list.size(); ichunk++) {
	for (ssize_t it_chunk = 0; it_chunk < chunk_nt_list[ichunk]; it_chunk++) {
	    ssize_t it_file = curr_it_file + it_chunk;
	    ssize_t it_logical = chunk_ipos_list[ichunk] - chunk_ipos_list[0] + it_chunk;
	    ssize_t it_time = chunk_ipos_list[ichunk] + it_chunk;
	    double expected_time = initial_time + dt_sample * (chunk_ipos_list[ichunk] + it_chunk);

	    assert(f2.time_index_mapping[it_file] == it_logical);
	    assert(fabs(f2.times[it_file] - expected_time) < 1.0e-7);

	    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
		for (int ipol = 0; ipol < npol; ipol++) {
		    int i = ifreq*npol*nt_file + ipol*nt_file + it_file;
		    // cerr << "ifreq=" << ifreq << " ipol=" << ipol << " it=" << it_file << " file_int=" << f2.intensity[i] << " exp_int=" << intval(ifreq,ipol,it_logical) << endl;
		    assert(fabs(f2.intensity[i] - intval(ifreq,ipol,it_time)) < 1.0e-3);
		    assert(fabs(f2.weights[i] - wtval(ifreq,ipol,it_time)) < 1.0e-3);
		}
	    }
	}

	curr_it_file += chunk_nt_list[ichunk];
    }

    assert(curr_it_file == nt_file);

    err = unlink(filename2);
    if (err)
	throw runtime_error("couldn't unlink file");

    cerr << "deleted " << filename2 << "\n"
	 << "Test passed!\n";
    
    return 0;
}
