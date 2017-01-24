#include <cstdio>
#include <iostream>
#include <stdarg.h> // for va_start/va_end
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


assembled_chunk::assembled_chunk(int beam_id_, int nupfreq_, int nt_per_packet_, int fpga_counts_per_sample_, uint64_t ichunk_)
    : beam_id(beam_id_), 
      nupfreq(nupfreq_), 
      nt_per_packet(nt_per_packet_),
      fpga_counts_per_sample(fpga_counts_per_sample_), 
      nt_coarse(constants::nt_per_assembled_chunk / nt_per_packet),
      nscales(constants::nfreq_coarse_tot * nt_coarse),
      ndata(constants::nfreq_coarse_tot * nupfreq * constants::nt_per_assembled_chunk),
      ichunk(ichunk_),
      isample(ichunk * constants::nt_per_assembled_chunk)
{
    if ((beam_id < 0) || (beam_id > constants::max_allowed_beam_id))
	throw runtime_error("assembled_chunk constructor: bad beam_id argument");
    if ((nupfreq <= 0) || (nupfreq > constants::max_allowed_nupfreq))
	throw runtime_error("assembled_chunk constructor: bad nupfreq argument");
    if ((nt_per_packet <= 0) || !is_power_of_two(nt_per_packet) || (nt_per_packet > constants::nt_per_assembled_chunk))
	throw runtime_error("assembled_chunk constructor: bad nt_per_packet argument");
    if ((fpga_counts_per_sample <= 0) || (fpga_counts_per_sample > constants::max_allowed_fpga_counts_per_sample))
	throw runtime_error("assembled_chunk constructor: bad fpga_counts_per_sample argument");

    this->scales = aligned_alloc<float> (nscales);
    this->offsets = aligned_alloc<float> (nscales);
    this->data = aligned_alloc<uint8_t> (ndata);
}

assembled_chunk::~assembled_chunk()
{
    free(data);
    free(scales);
    free(offsets);
}

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
    //   (BEAM)    -> %04i beam_id
    //   (CHUNK)   -> %08i ichunk
    //   (NCHUNK)  -> %02i  size in chunks
    //   (BINNING) -> %02i  size in chunks
    //   (FPGA0)   -> %012i start FPGA-counts
    //   (FPGAN)   -> %08i  FPGA-counts size
    string s = pattern;
    s = replaceAll(s, "(BEAM)",    stringprintf("%04i",        beam_id));
    s = replaceAll(s, "(CHUNK)",   stringprintf("%08"  PRIu64, ichunk));
    s = replaceAll(s, "(NCHUNK)",  stringprintf("%02i",        binning));
    s = replaceAll(s, "(BINNING)", stringprintf("%02i",        binning));
    s = replaceAll(s, "(FPGA0)",   stringprintf("%012" PRIu64, fpgacounts_begin()));
    s = replaceAll(s, "(FPGAN)",   stringprintf("%08"  PRIu64, fpgacounts_N()));
    return s;
}

void assembled_chunk::fill_with_copy(const shared_ptr<assembled_chunk> &x)
{
    if (!x)
	throw runtime_error("assembled_chunk::fill_with_copy() called with empty pointer");
    if ((this->nupfreq != x->nupfreq) || (this->nt_per_packet != x->nt_per_packet))
	throw runtime_error("assembled_chunk::fill_with_copy() called on non-conformable chunks");

    if (x.get() == this)
	return;

    memcpy(this->data, x->data, ndata);
    memcpy(this->scales, x->scales, nscales * sizeof(float));
    memcpy(this->offsets, x->offsets, nscales * sizeof(float));
}


// Used in unit tests
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
}


// virtual member function; any changes made here should be reflected in override fast_assembled_chunk::decode().
void assembled_chunk::decode(float *intensity, float *weights, int stride) const
{
    if (!intensity || !weights)
	throw runtime_error("ch_frb_io: null pointer passed to assembled_chunk::decode()");	
    if (stride < constants::nt_per_assembled_chunk)
	throw runtime_error("ch_frb_io: bad stride passed to assembled_chunk::decode()");

    for (int if_coarse = 0; if_coarse < constants::nfreq_coarse_tot; if_coarse++) {
	const float *scales_f = this->scales + if_coarse * nt_coarse;
	const float *offsets_f = this->offsets + if_coarse * nt_coarse;
	
	for (int if_fine = if_coarse*nupfreq; if_fine < (if_coarse+1)*nupfreq; if_fine++) {
	    const uint8_t *src_f = this->data + if_fine * constants::nt_per_assembled_chunk;
	    float *int_f = intensity + if_fine * stride;
	    float *wt_f = weights + if_fine * stride;

	    for (int it_coarse = 0; it_coarse < nt_coarse; it_coarse++) {
		float scale = scales_f[it_coarse];
		float offset = offsets_f[it_coarse];
		
		for (int it_fine = it_coarse*nt_per_packet; it_fine < (it_coarse+1)*nt_per_packet; it_fine++) {
		    float x = float(src_f[it_fine]);
		    int_f[it_fine] = scale*x + offset;
		    wt_f[it_fine] = ((x==0) || (x==255)) ? 0.0 : 1.0;
		}
	    }
	}
    }
}

void assembled_chunk::decode_subset(float *intensity, float *weights,
                                    int t0, int NT, int stride) const {
    if (!intensity || !weights)
	throw runtime_error("ch_frb_io: null pointer passed to assembled_chunk::decode_subset()");
    if (stride < NT)
	throw runtime_error("ch_frb_io: bad stride passed to assembled_chunk::decode_subset()");
    if ((t0 < 0) || (NT < 0) || (t0 + NT > constants::nt_per_assembled_chunk))
	throw runtime_error("ch_frb_io: bad (t0,NT) passed to assembled_chunk::decode_subset()");

    for (int if_coarse = 0; if_coarse < constants::nfreq_coarse_tot; if_coarse++) {
	const float * scales_f = this->scales  + if_coarse * nt_coarse;
	const float *offsets_f = this->offsets + if_coarse * nt_coarse;

	for (int if_fine = if_coarse*nupfreq; if_fine < (if_coarse+1)*nupfreq; if_fine++) {
	    const uint8_t *src_f = this->data + if_fine * constants::nt_per_assembled_chunk;
	    float *int_f = intensity + if_fine * stride;
	    float * wt_f = weights   + if_fine * stride;

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

assembled_chunk* assembled_chunk::downsample(assembled_chunk* dest,
                                             const assembled_chunk* src1,
                                             const assembled_chunk* src2) {

    if (src1->beam_id != src2->beam_id)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched beam_id");
    if (src1->nupfreq != src2->nupfreq)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched nupfreq");
    if (src1->nt_coarse != src2->nt_coarse)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched nt_coarse");
    if (src1->nt_per_packet != src2->nt_per_packet)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched nt_per_packet");
    if (src1->nscales != src2->nscales)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched nscales");
    if (src1->ndata != src2->ndata)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched ndata");

    if (src1->binning != src2->binning)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched binning");

    if (src1->ichunk >= src2->ichunk)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: expected src1 to have earlier ichunk than src2");

    if (!dest) {
        unique_ptr<assembled_chunk> up = assembled_chunk::make(src1->beam_id, src1->nupfreq, src1->nt_per_packet, 2 * src1->fpga_counts_per_sample, src1->ichunk);
        dest = up.release();
    }

    if (src1->beam_id != dest->beam_id)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched dest beam_id");
    if (src1->nupfreq != dest->nupfreq)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched dest nupfreq");
    if (src1->nt_coarse != dest->nt_coarse)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched dest nt_coarse");
    if (src1->nt_per_packet != dest->nt_per_packet)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched dest nt_per_packet");
    if (src1->nscales != dest->nscales)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched dest nscales");
    if (src1->ndata != dest->ndata)
        throw runtime_error("ch_frb_io: assembled_chunk::downsample: mismatched dest ndata");

    float* dest_scales  = dest->scales;
    float* dest_offsets = dest->offsets;
    bool temp_scales = false;
    if (src1 == dest) {
        // enable in-place downsampling
        dest_scales  = (float*)malloc(src1->nscales * sizeof(float));
        dest_offsets = (float*)malloc(src1->nscales * sizeof(float));
        temp_scales = true;
    }

    // Compute destination offset + scale.
    // first half of destination offset + scale comes from src1.
    for (int i=0; i<src1->nscales/2; i++) {
        float scale_1  = src1->scales [2*i];
        float offset_1 = src1->offsets[2*i];
        float scale_2  = src1->scales [2*i+1];
        float offset_2 = src1->offsets[2*i+1];

        // Lower and upper limits of each offset, scale
        // choice... ignoring the fact that 0 and 255 are
        // marker values for masked values.
        float lo = min(offset_1, offset_2);
        float hi = max(offset_1 + scale_1 * 255., offset_2 + scale_2 * 255.);

        // Make the new range contain the old range.  (The
        // data may fill a smaller range than this...)
        dest_scales [i]  = (hi - lo) / 255.;
        dest_offsets[i] = lo;
    }
    // second half of destination offset + scale comes from src2.
    for (int i=0; i<src2->nscales/2; i++) {
        float scale_1  = src2->scales [2*i];
        float offset_1 = src2->offsets[2*i];
        float scale_2  = src2->scales [2*i+1];
        float offset_2 = src2->offsets[2*i+1];

        // Lower and upper limits of each offset, scale
        // choice... ignoring the fact that 0 and 255 are
        // marker values for masked values.
        float lo = min(offset_1, offset_2);
        float hi = max(offset_1 + scale_1 * 255., offset_2 + scale_2 * 255.);

        // Make the new range contain the old range.  (The
        // data may fill a smaller range than this...)
        dest_scales [src1->nscales/2 + i]  = (hi - lo) / 255.;
        dest_offsets[src1->nscales/2 + i] = lo;
    }

    const int nupfreq = src1->nupfreq;
    const int nt_coarse = src1->nt_coarse;
    const int nt_per_packet = src1->nt_per_packet;

    for (int if_coarse = 0; if_coarse < constants::nfreq_coarse_tot; if_coarse++) {
	const float *scales_1  = src1->scales  + if_coarse * nt_coarse;
	const float *offsets_1 = src1->offsets + if_coarse * nt_coarse;
	const float *scales_2  = src2->scales  + if_coarse * nt_coarse;
	const float *offsets_2 = src2->offsets + if_coarse * nt_coarse;
	const float *scales_d  = dest_scales   + if_coarse * nt_coarse;
	const float *offsets_d = dest_offsets  + if_coarse * nt_coarse;

	for (int if_fine = if_coarse*nupfreq; if_fine < (if_coarse+1)*nupfreq; if_fine++) {
	    const uint8_t *data_1 = src1->data + if_fine * constants::nt_per_assembled_chunk;
	    const uint8_t *data_2 = src2->data + if_fine * constants::nt_per_assembled_chunk;

	    uint8_t *data_d = dest->data + if_fine * constants::nt_per_assembled_chunk;
            // First half of data comes from src1
	    for (int it_coarse = 0; it_coarse < nt_coarse; it_coarse++) {
		float scale  = scales_1 [it_coarse];
		float offset = offsets_1[it_coarse];

                float iscale_d = 1./scales_d[it_coarse / 2];
                float offset_d =   offsets_d[it_coarse / 2];

		for (int it_fine = (it_coarse*nt_per_packet)/2; it_fine < ((it_coarse+1)*nt_per_packet)/2; it_fine++) {
		    uint8_t d1 = data_1[it_fine * 2    ];
		    uint8_t d2 = data_1[it_fine * 2 + 1];
                    int wtd = 0;
                    float xd = 0;
                    if (!(d1 == 0 || d1 == 255)) {
                        wtd++;
                        xd += offset + (float)d1 * scale;
                    }
                    if (!(d2 == 0 || d2 == 255)) {
                        wtd++;
                        xd += offset + (float)d2 * scale;
                    }
                    if (wtd == 0)
                        data_d[it_fine] = 0;
                    else {
                        xd /= (float)wtd;
                        xd = (xd - offset_d) * iscale_d;
                        // FIXME -- round?
                        data_d[it_fine] = (uint8_t)lround(xd);
                    }
		}
	    }

	    data_d += constants::nt_per_assembled_chunk / 2;
            // Second half of data comes from src2
	    for (int it_coarse = 0; it_coarse < nt_coarse; it_coarse++) {
		float scale  = scales_2 [it_coarse];
		float offset = offsets_2[it_coarse];

                float iscale_d = 1./scales_d[nt_coarse / 2 + it_coarse / 2];
                float offset_d =   offsets_d[nt_coarse / 2 + it_coarse / 2];

		for (int it_fine = (it_coarse*nt_per_packet)/2; it_fine < ((it_coarse+1)*nt_per_packet)/2; it_fine++) {
		    uint8_t d1 = data_2[it_fine * 2    ];
		    uint8_t d2 = data_2[it_fine * 2 + 1];
                    int wtd = 0;
                    float xd = 0;
                    if (!(d1 == 0 || d1 == 255)) {
                        wtd++;
                        xd += offset + (float)d1 * scale;
                    }
                    if (!(d2 == 0 || d2 == 255)) {
                        wtd++;
                        xd += offset + (float)d2 * scale;
                    }
                    if (wtd == 0)
                        data_d[it_fine] = 0;
                    else {
                        xd /= (float)wtd;
                        xd = (xd - offset_d) * iscale_d;
                        // FIXME -- round?
                        data_d[it_fine] = (uint8_t)lround(xd);
                    }
		}
	    }
	}
    }
    if (temp_scales) {
        memcpy(dest->scales,  dest_scales,  src1->nscales * sizeof(float));
        memcpy(dest->offsets, dest_offsets, src1->nscales * sizeof(float));
        free(dest_scales);
        free(dest_offsets);

        // When downsampling in place, update the sampling.
        dest->fpga_counts_per_sample = 2 * src1->fpga_counts_per_sample;
    }
    dest->binning = src1->binning * 2;

    return dest;
}


unique_ptr<assembled_chunk> assembled_chunk::make(int beam_id_, int nupfreq_, int nt_per_packet_, int fpga_counts_per_sample_, uint64_t ichunk_, bool force_reference, bool force_fast)
{
    // FIXME -- if C++14 is available, use make_unique()
    if (force_reference && force_fast)
        throw runtime_error("ch_frb_io: assembled_chunk::make(): both force_reference and force_fast were set!");

#ifdef __AVX2__
    if (force_fast ||
        ((nt_per_packet_ == 16) && (nupfreq_ % 2 == 0) && !force_reference))
	return unique_ptr<fast_assembled_chunk>(new fast_assembled_chunk(beam_id_, nupfreq_, nt_per_packet_, fpga_counts_per_sample_, ichunk_));
#else
    if (force_fast)
        throw runtime_error("ch_frb_io: assembled_chunk::make(): force_fast set on a machine without AVX2!");
#endif

    return unique_ptr<assembled_chunk>(new assembled_chunk(beam_id_, nupfreq_, nt_per_packet_, fpga_counts_per_sample_, ichunk_));
}


void assembled_chunk::write_hdf5_file(const string &filename)
{
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

void assembled_chunk::write_msgpack_file(const string &filename)
{
    FILE* f = fopen(filename.c_str(), "w+");
    if (!f)
        throw runtime_error("ch_frb_io: failed to open file " + filename + " for writing an assembled_chunk in msgpack format: " + strerror(errno));
    // msgpack buffer that will write to file "f"
    msgpack::fbuffer buffer(f);
    // Construct a shared_ptr from this, carefully
    shared_ptr<assembled_chunk> shthis(shared_ptr<assembled_chunk>(), this);
    msgpack::pack(buffer, shthis);
    if (fclose(f))
        throw runtime_error("ch_frb_io: failed to close assembled_chunk msgpack file " + filename + string(strerror(errno)));
}

shared_ptr<assembled_chunk> assembled_chunk::read_msgpack_file(const string &filename)
{
    struct stat st;
    if (stat(filename.c_str(), &st)) {
        throw runtime_error("ch_frb_io: failed to stat file " + filename + " for reading an assembled_chunk in msgpack format: " + strerror(errno));
    }
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


}  // namespace ch_frb_io
