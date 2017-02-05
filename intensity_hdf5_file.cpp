#include <iostream>
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


// Helper function for intensity_hdf5_file constructor
static void _init_frequencies(intensity_hdf5_file &f, hdf5_group &g)
{
    vector<hsize_t> freq_shape;
    g.get_dataset_shape("freq", freq_shape);

    if (freq_shape.size() != 1)
	throw runtime_error(f.filename + ": expected index_map.freq.ndim == 1");

    f.nfreq = freq_shape[0];

    vector<double> frequency_Hz(f.nfreq);
    g.read_dataset("freq", &frequency_Hz[0], freq_shape);

    if (f.nfreq < 2)
	throw runtime_error(f.filename + ": expected nfreq >= 2");

    double nu_lo = 1.0e-6 * frequency_Hz[0];          // Hz -> MHz
    double nu_hi = 1.0e-6 * frequency_Hz[f.nfreq-1];  // Hz -> MHz
    nu_hi += (nu_hi - nu_lo) / (f.nfreq-1);           // actual endpoint of band

    for (int i = 0; i < f.nfreq; i++) {
	// Expected and actual frequency
	double nu_exp = ((f.nfreq-i)*nu_lo + i*nu_hi) / (double)f.nfreq;
	double nu_act = 1.0e-6 * frequency_Hz[i];   // Hz -> MHz

	if (fabs(nu_exp-nu_act) > 1.0e-6 * fabs(nu_hi-nu_lo))
	    throw runtime_error(f.filename + ": expected equally spaced frequencies");
    }

    if (nu_lo == nu_hi)
	throw runtime_error(f.filename + ": frequencies in file are all equal?!");

    f.frequencies_are_increasing = (nu_lo < nu_hi);
    f.freq_lo_MHz = min(nu_lo, nu_hi);
    f.freq_hi_MHz = max(nu_lo, nu_hi);    
}


// Helper function for intensity_hdf5_file constructor
static void _init_polarizations(intensity_hdf5_file &f, hdf5_group &g)
{
    vector<hsize_t> pol_shape;
    g.get_dataset_shape("pol", pol_shape);

    if (pol_shape.size() != 1)
	throw runtime_error(f.filename + ": expected index_map.pol.ndim == 1");
    if ((pol_shape[0] != 1) && (pol_shape[0] != 2))
	throw runtime_error(f.filename + ": expected npol==1 or npol==2");

    f.npol = pol_shape[0];
    g.read_string_dataset("pol", f.polarizations, pol_shape);

    for (int ipol = 0; ipol < f.npol; ipol++) {
	if ((f.polarizations[ipol] != "XX") && (f.polarizations[ipol] != "YY"))
	    throw runtime_error(f.filename + ": expected polarization string to be either 'XX' or 'YY', got '" + f.polarizations[ipol] + "'");

	for (int jpol = 0; jpol < ipol; jpol++)
	    if (f.polarizations[ipol] == f.polarizations[jpol])
		throw runtime_error(f.filename + ": duplicate polarizations in file?!");
    }
}


// Helper function for intensity_hdf5_file constructor
static void _init_times(intensity_hdf5_file &f, hdf5_group &g)
{
    vector<hsize_t> time_shape;
    g.get_dataset_shape("time", time_shape);
    
    if (time_shape.size() != 1)
	throw runtime_error(f.filename + ": expected index_map.time.ndim == 1");

    f.nt_file = time_shape[0];

    if (f.nt_file < 2)
	throw runtime_error(f.filename + ": expected nt >= 2");

    f.times.resize(f.nt_file);
    g.read_dataset("time", &f.times[0], time_shape);

    // provisional dt_sample = min_i(dt[i+1]-dt[i])
    double dt_p = f.times[1] - f.times[0];
    for (int i = 2; i < f.nt_file; i++)
	dt_p = min(dt_p, f.times[i] - f.times[i-1]);

    if (dt_p <= 0.0)
	throw runtime_error(f.filename + ": expected time array to be monotone increasing");

    f.time_index_mapping.resize(f.nt_file, 0);

    for (int i = 1; i < f.nt_file; i++) {
	double dt = f.times[i] - f.times[i-1];
	if (dt > 1.0e6 * dt_p)
	    throw runtime_error(f.filename + ": file has unreasonably large time gap");

	int nstep = (int)(dt/dt_p + 0.5);
	if ((nstep < 1) || (nstep > 1000000))
	    throw runtime_error(f.filename + ": internal error ('impossible' nstep)");

	if (fabs(dt - nstep*dt_p) > 0.01*dt_p) {
	    cerr << "debug: " << dt << " " << dt_p << " " << nstep << endl;
	    throw runtime_error(f.filename + ": hmmm, time steps in file do not appear to be integer multiples of a common step size ");
	}

	f.time_index_mapping[i] = f.time_index_mapping[i-1] + nstep;
	if (f.time_index_mapping[i] > 100000000)
	    throw runtime_error(f.filename + ": file is unreasonably gappy");
    }

    f.nt_logical = f.time_index_mapping[f.nt_file-1] + 1;
    f.time_lo = f.times[0];
    f.time_hi = f.times[f.nt_file-1] + dt_p;
    f.dt_sample = (f.time_hi - f.time_lo) / f.nt_logical;
}


intensity_hdf5_file::intensity_hdf5_file(const string &filename_, bool noisy) :
    filename(filename_)
{
    hdf5_file f(filename);
    hdf5_group g_root(f, ".");
    hdf5_group g_im(f, "index_map");

    _init_frequencies(*this, g_im);
    _init_polarizations(*this, g_im);
    _init_times(*this, g_im);

    vector<hsize_t> data_shape(3);
    data_shape[0] = nfreq;
    data_shape[1] = npol;
    data_shape[2] = nt_file;

    this->intensity.resize(nfreq * npol * nt_file);
    this->weights.resize(nfreq * npol * nt_file);

    g_root.read_dataset("intensity", &intensity[0], data_shape);
    g_root.read_dataset("weight", &weights[0], data_shape);

    // Note: these need to be double-precision
    double wsum = 0.0;
    double wmax = 0.0;

    for (int i = 0; i < nfreq*npol*nt_file; i++) {
	if (weights[i] < 0.0)
	    throw runtime_error(filename + ": negative value in weight array, this is currently considered an error");
	wsum += weights[i];
	wmax = max(wmax, (double)weights[i]);
    }

    if (wmax == 0.0)
	cerr << (filename + ": warning: weight array is all zeros\n");

    this->frac_ungapped = (double)nt_file / (double)nt_logical;
    this->frac_unmasked = (wmax > 0.0) ? (wsum / wmax / (float)(nfreq*npol*nt_file)) : 0.0;

    if (noisy)
	cerr << "read " << filename << ", frac_ungapped=" << frac_ungapped << ", frac_unmasked=" << frac_unmasked << endl;
}


void intensity_hdf5_file::get_unpolarized_intensity(float *out_int, float *out_wt, int out_t0, int out_nt, int out_stride) const
{
    if (out_stride == 0)
	out_stride = out_nt;

    if (out_int==nullptr || out_wt==nullptr)
	throw runtime_error("NULL pointer passed to intensity_hdf5_file::get_unpolarized_intensity()");
    if (out_nt <= 0)
	throw runtime_error("intensity_hdf5_file::get_unpolarized_intensity(): expected out_nt > 0");
    if ((out_t0 < 0) || (out_t0 + out_nt > this->nt_logical))
	throw runtime_error("intensity_hdf5_file::get_unpolarized_intensity(): requested time range overflows file");
    if (abs(out_stride) < out_nt)
	throw runtime_error("intensity_hdf5_file::get_unpolarized_intensity(): out_stride is too small");

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	memset(out_int + ifreq*out_stride, 0, out_nt * sizeof(*out_int));
	memset(out_wt + ifreq*out_stride, 0, out_nt * sizeof(*out_wt));
    }

    //
    // First pass: dst_int = (weight*intensity) and dst_wt = (weight).
    // FIXME this loop could be written significantly more efficiently!
    //

    for (int it_f = 0; it_f < nt_file; it_f++) {
	int it_l = time_index_mapping[it_f];

	if ((it_l < out_t0) || (it_l >= out_t0+out_nt))
	    continue;

	// 1D arrays with length=nfreq and stride=out_stride
	float *dst_int = out_int + (it_l - out_t0);
	float *dst_wt = out_wt + (it_l - out_t0);

	for (int ipol = 0; ipol < npol; ipol++) {
	    // 1D arrays with length=nfreq and stride=in_stride
	    const float *src_int = &this->intensity[0] + ipol*nt_file + it_f;
	    const float *src_wt = &this->weights[0] + ipol*nt_file + it_f;
	    const int in_stride = npol * nt_file;

	    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
		int d = ifreq * out_stride;
		int s = ifreq * in_stride;
		dst_int[d] += src_int[s] * src_wt[s];
		dst_wt[d] += src_wt[s];
	    }
	}
    }

    // Second pass: dst_int = (intensity)

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	// 1D arrays with length=out_nt
	float *dst_int = out_int + ifreq * out_stride;
	float *dst_wt = out_wt + ifreq * out_stride;

	for (int it = 0; it < out_nt; it++) {
	    if (dst_wt[it] > 0.0)
		dst_int[it] /= dst_wt[it];
	}
    }
}


void intensity_hdf5_file::run_unit_tests() const
{
    cerr << filename << ": starting unit test\n";

    std::random_device rd;
    std::mt19937 rng(rd());

    for (int it_f = 0; it_f < nt_file; it_f++) {
	int it_l = time_index_mapping[it_f];
	if ((it_l < 0) || (it_l >= nt_logical))
	    throw runtime_error(filename + ": time_index_mapping index out of range");
	
	// Expected and actual times
	double t_exp = ((nt_logical-it_l)*time_lo + it_l*time_hi) / (double)nt_logical;
	double t_act = times[it_f];

	if (fabs(t_exp-t_act) > 0.01*dt_sample)
	    throw runtime_error(filename + ": timestamp check failed");
    }    

    // The rest of this routine tests get_unpolarized_intensity()

    // First we populate "reference" buffers for the whole time range,
    // using the most boneheaded procedure possible.

    vector<float> ref_int(nfreq * nt_logical, 0.0);
    vector<float> ref_wt(nfreq * nt_logical, 0.0);

    for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	for (int ipol = 0; ipol < npol; ipol++) {
	    for (int it_f = 0; it_f < nt_file; it_f++) {
		int it_l = time_index_mapping[it_f];
		int s = ifreq*npol*nt_file + ipol*nt_file + it_f;
		int d = ifreq*nt_logical + it_l;

		ref_int[d] += weights[s] * intensity[s];
		ref_wt[d] += weights[s];
	    }
	}
    }

    for (int i = 0; i < nfreq*nt_logical; i++) {
	if (ref_wt[i] > 0.0)
	    ref_int[i] /= ref_wt[i];
    }

    // Next we make some random calls to intensity_hdf5_file::get_unpolarized_intensity()
    // and compare to the reference buffers.

    vector<float> v_int(4 * nfreq * nt_logical, 0.0);
    vector<float> v_wt(4 * nfreq * nt_logical, 0.0);

    // These pointers have been chosen so that strides in the range [-2*nt_logical,2*nt_logical] will not overflow
    float *buf_int = &v_int[0] + 2 * nfreq * nt_logical;
    float *buf_wt = &v_wt[0] + 2 * nfreq * nt_logical;

    for (int n = 0; n < 100; n++) {
	cerr << ".";

	// make random (it0, nt)
	int buf_it0 = randint(rng, 0, nt_logical-1);
	int buf_it1 = randint(rng, 0, nt_logical-1);
	if (buf_it0 > buf_it1) std::swap(buf_it0, buf_it1);
	int buf_nt = buf_it1 - buf_it0 + 1;

	// random stride, can be positive or negative
	int stride = randint(rng, buf_nt, 2*nt_logical);
	if (uniform_rand(rng) > 0.5) stride *= -1;
	
	// randomize buffer
	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    uniform_rand(rng, buf_int + ifreq*stride, buf_nt);
	    uniform_rand(rng, buf_wt + ifreq*stride, buf_nt);
	}

	// randomized call to intensity_hdf5_file::get_unpolarized_intensity()
	this->get_unpolarized_intensity(buf_int, buf_wt, buf_it0, buf_nt, stride);

	// compare with reference
	for (int ifreq = 0; ifreq < nfreq; ifreq++) {
	    for (int it = 0; it < buf_nt; it++) {
		int sr = ifreq*nt_logical + it + buf_it0;
		int sb = ifreq*stride + it;

		if (fabs(ref_int[sr] - buf_int[sb]) > 1.0e-3)
		    throw runtime_error(filename + ": intensity consistency check failed");
		if (fabs(ref_wt[sr] - buf_wt[sb]) > 1.0e-3)
		    throw runtime_error(filename + ": weights consistency check failed");
	    }
	}
    }

    cerr << "\n" << filename << ": pass\n";
}


}  // namespace ch_frb_io
