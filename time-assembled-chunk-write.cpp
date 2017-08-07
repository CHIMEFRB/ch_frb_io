#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <dirent.h>
#include <cstring>
#include <sstream>

#include "ch_frb_io_internals.hpp"


using namespace std;
using namespace ch_frb_io;


// -------------------------------------------------------------------------------------------------
//
// Utils (all cut-and-paste from elsewhere)


static void delete_file(const string &filename)
{
    int err = unlink(filename.c_str());
    if (err < 0)
	throw runtime_error(filename + ": " + strerror(errno));
}


static bool is_empty_directory(const string &dirname)
{
    DIR *dir = opendir(dirname.c_str());
    if (!dir)
	throw runtime_error(dirname + ": opendir() failed: " + strerror(errno));

    ssize_t name_max = pathconf(dirname.c_str(), _PC_NAME_MAX);
    name_max = min(name_max, (ssize_t)4096);

    vector<char> buf(sizeof(struct dirent) + name_max + 1);
    struct dirent *entry = reinterpret_cast<struct dirent *> (&buf[0]);
    
    for (;;) {
	struct dirent *result = nullptr;

	int err = readdir_r(dir, entry, &result);	
	if (err)
	    throw runtime_error(dirname + ": readdir_r() failed");
	if (!result)
	    return true;
	if (!strcmp(entry->d_name, "."))
	    continue;
	if (!strcmp(entry->d_name, ".."))
	    continue;
	
	return false;
    }
}


static void write_file(const string &filename, const void *buf, ssize_t count, bool clobber=true)
{
    // This cast to (uint8_t *) suppresses a superfluous compiler warning below.
    uint8_t *p = (uint8_t *) buf;
    
    if (count < 0)
	throw runtime_error("write_file(): expected count >= 0");
    if (count && !p)
	throw runtime_error("write_file(): 'buf' is a null pointer");
	
    int flags = O_WRONLY | O_CREAT | O_TRUNC;
    if (!clobber)
	flags |= O_EXCL;

    // Reasonable default?  owner=rw, group=r, other=r (note that umask will also be applied)
    mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;

    int fd = open(filename.c_str(), flags, mode);
    if (fd < 0)
	throw runtime_error(filename + ": open() failed: " + strerror(errno));

    // We loop over multiple calls to write() because, according to the manpage,
    // write() is allowed to write a subset of the buffer and return the number of
    // bytes actually written.  I've never seen this actually happen though!

    while (count > 0) {
	ssize_t n = write(fd, p, count);

	if (n <= 0) {
	    close(fd);
	    const char *msg = (n < 0) ? strerror(errno) : "write() returned 0?!";
	    throw runtime_error(filename + ": write() failed: " + msg);
	}
	
	count -= n;
	p += n;
    }

    close(fd);
}


static ssize_t get_physical_memory()
{
    ssize_t pagesize = sysconf(_SC_PAGESIZE);
    if (pagesize < 0)
	throw runtime_error(string("sysconf(_SC_PAGESIZE) failed: ") + strerror(errno));

    ssize_t npages = sysconf(_SC_PHYS_PAGES);
    if (npages < 0)
	throw runtime_error(string("sysconf(_SC_PHYS_PAGES) failed: ") + strerror(errno));

    return npages * pagesize;
}


static ssize_t get_file_size(const string &filename)
{
    struct stat s;

    int err = stat(filename.c_str(), &s);
    if (err < 0)
	throw runtime_error(filename + ": " + strerror(errno));

    return s.st_size;
}


static void sync_filesystem(const string &filename)
{
    int fd = open(filename.c_str(), O_RDONLY, 0);
    if (fd < 0)
	throw runtime_error(filename + ": open() failed: " + strerror(errno));

    int err = syncfs(fd);
    close(fd);
    
    if (err < 0)
	throw runtime_error(filename + ": syncfs() failed: " + strerror(errno));
}


// -------------------------------------------------------------------------------------------------
//
// "Boneheaded" write


static void boneheaded_write(const string &filename, const unique_ptr<assembled_chunk> &chunk, uint8_t *tmp_buffer)
{
    int *enc0 = (int *) tmp_buffer;
    *enc0++ = 1;   // version number
    *enc0++ = chunk->beam_id;
    *enc0++ = chunk->nupfreq;
    *enc0++ = chunk->nt_per_packet;
    *enc0++ = chunk->fpga_counts_per_sample;
    *enc0++ = chunk->binning;

    uint64_t *enc1 = (uint64_t *) enc0;
    *enc1++ = chunk->ichunk;

    float *enc2 = (float *) enc1;
    memcpy(enc2, chunk->scales, chunk->nscales * sizeof(float));
    enc2 += chunk->nscales;

    memcpy(enc2, chunk->offsets, chunk->nscales * sizeof(float));
    enc2 += chunk->nscales;

    uint8_t *enc3 = (uint8_t *) enc2;
    memcpy(enc3, chunk->data, chunk->ndata);
    enc3 += chunk->ndata;
    
    ssize_t nbytes = (enc3 - tmp_buffer);  // OK since buffer,enc3 are both (uint8_t *)
    write_file(filename, tmp_buffer, nbytes);
}


// -------------------------------------------------------------------------------------------------


static void usage(const char *msg = nullptr)
{
    cerr << "Usage: time-assembled-chunk-write [-zbu] <target_dir> <target_gb>\n"
	 << "   -z uses zeroed chunks (default is to randomize)\n"
	 << "   -b uses boneheaded file format (default is compressed msgpack)\n"
	 << "   -u uses uncompressed msgpack file format (default is compressed msgpack)\n"
	 << "   -p pins process to core 0 (default is unpinned)\n";

    if (msg)
	cerr << "Fatal: " << msg << endl;

    exit(2);
}


static void usage(const string &msg)
{
    usage(msg.c_str());
}


int main(int argc, char **argv)
{
    // Chunk parameters expected for full CHIME.
    constexpr int nupfreq = 16;
    constexpr int nt_per_packet = 16;
    constexpr int fpga_counts_per_sample = 384;

    // "Derived" chunk parameters.
    constexpr ssize_t nfreq_c = ch_frb_io::constants::nfreq_coarse_tot;
    constexpr ssize_t nfreq_f = nfreq_c * nupfreq;
    constexpr ssize_t nt_f = ch_frb_io::constants::nt_per_assembled_chunk;
    constexpr ssize_t nt_c = nt_f / nt_per_packet;
    constexpr ssize_t header_nbytes = 2 * nfreq_c * nt_c * sizeof(float);
    constexpr ssize_t data_nbytes = nfreq_f * nt_f;
    constexpr double gb_per_chunk = (header_nbytes + data_nbytes) / pow(2.,30.);

    // Parse command line

    vector<string> args;
    bool zflag = false;
    bool bflag = false;
    bool uflag = false;
    bool pflag = false;

    for (int i = 1; i < argc; i++) {
	const char *arg = argv[i];
	
	if (arg[0] != '-') {
	    args.push_back(arg);
	    continue;
	}

	int arglen = strlen(arg);
	if (arglen <= 1)
	    usage();

	for (int j = 1; j < arglen; j++) {
	    if (arg[j] == 'z')
		zflag = true;
	    else if (arg[j] == 'b')
		bflag = true;
	    else if (arg[j] == 'u')
		uflag = true;
	    else if (arg[j] == 'p')
		uflag = true;
	    else
		usage();
	}
    }

    if (args.size() != 2) 
	usage();

    string target_dir = args[0];
    double target_gb = lexical_cast<double> (args[1]);

    if (!is_empty_directory(target_dir))
	usage(target_dir + " must be an empty directory");
    if (target_gb <= 0.0)
	usage("target_gb must be > 0");
    if (bflag && uflag)
	usage("-b with -u does not make sense");

    double gb_physical = get_physical_memory() / pow(2.,30.0);
    if (target_gb > gb_physical)
	usage("target_gb must be <= " + to_string(gb_physical) + " (total physical memory)");

    // Command-line parsing finished!

    std::random_device rd;
    std::mt19937 rng(rd());

    if (pflag) {
	int hwcon = std::thread::hardware_concurrency();
	pin_thread_to_cores({0,hwcon/2});
    }

    int nchunks = int(target_gb / gb_per_chunk) + 1;
    double actual_gb = nchunks * gb_per_chunk;

    cout << "Generating " << nchunks << " assembled_chunks (" << actual_gb << " GB), this may take some time" << endl;

    vector<unique_ptr<assembled_chunk>> chunks(nchunks);
    vector<string> filenames(nchunks);

    for (int i = 0; i < nchunks; i++) {
	assembled_chunk::initializer ini_params;
	ini_params.nupfreq = nupfreq;
	ini_params.nt_per_packet = nt_per_packet;
	ini_params.fpga_counts_per_sample = fpga_counts_per_sample;

	chunks[i] = assembled_chunk::make(ini_params);

	if (!zflag)
	    chunks[i]->randomize(rng);

	stringstream filename;
	filename << target_dir << "/deleteme_" << i << ".msgpack";
	filenames[i] = filename.str();

	if ((i % 100 == 0) && (i > 0))
	    cout << i << "/" << nchunks << " chunks generated" << endl;
    }
    
    // Temp buffer for encoding.
    // 32MB is overkill here!
    vector<uint8_t> tmp_buffer(32 * 1024 * 1024, 0);
    
    // Sync filesystem before writing.
    sync_filesystem(target_dir);

    cout << "Writing files";
    struct timeval tv0 = xgettimeofday();

    for (int i = 0; i < nchunks; i++) {
	bool compress = !bflag && !uflag;
	cout << "." << flush;

	if (bflag)
	    boneheaded_write(filenames[i], chunks[i], &tmp_buffer[0]);
	else
	    chunks[i]->write_msgpack_file(filenames[i], compress, &tmp_buffer[0]);
    }

    // Sync filesystem after writing.
    sync_filesystem(target_dir);

    struct timeval tv1 = xgettimeofday();
    double secs = 1.0e-6 * usec_between(tv0,tv1);

    cout << "done" << endl;
    cout << "Write speed: " << (actual_gb/secs) << " GB/sec" << endl;

    double compressed_gb = 0.0;
    for (int i = 0; i < nchunks; i++)
	compressed_gb += get_file_size(filenames[i]) / pow(2.,30.);

    cout << "Actual file size / Nominal size = " << (compressed_gb / actual_gb) << endl;
    cout << "Deleting files..." << flush;

    for (int i = 0; i < nchunks; i++)
	delete_file(filenames[i]);

    cout << "done" << endl;
    return 0;
}
