#include <sys/types.h>
#include <algorithm>

#include "ch_frb_io_internals.hpp"

#include <sys/stat.h>

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif

std::shared_ptr<slow_pulsar_chunk> slow_pulsar_chunk::make_slow_pulsar_chunk(std::shared_ptr<ch_chunk_initializer> ini_params)
{
	if(!ini_params->pool)
		throw std::runtime_error("slow_pulsar_chunk: memory_slab_pool required to construct new chunk");
	
	if(!ini_params->slab){
		ini_params->slab = ini_params->pool->get_slab(true);

		if (!ini_params->slab)
	    	throw std::runtime_error("**** Too much memory pressure for this poor L1 node to survive!  Blowing up now... ****");
	}

	return std::shared_ptr<slow_pulsar_chunk>(new slow_pulsar_chunk(ini_params));
}

slow_pulsar_chunk::slow_pulsar_chunk(const std::shared_ptr<ch_chunk_initializer> ini_params) : ch_chunk(*ini_params)
{
	if(!ini_params->pool)
		throw std::runtime_error("slow_pulsar_chunk constructor: 'pool' required to construct");
	this->memory_pool = ini_params->pool;

	if(ini_params->slab){
		this->memory_slab.swap(ini_params->slab);
	}
	else{
		uint8_t *p = aligned_alloc<uint8_t> (ini_params->pool->nbytes_per_slab);
    	this->memory_slab = memory_slab_t(p);
	}
    // TODO add size check logic to ensure minimum 
    // if (ini_params.pool->nbytes_per_slab < mc.slab_size) {
    //     throw std::runtime_error("slow_pulsar_chunk constructor: memory_slab_pool::nbytes_per_slab (=" 
    //             + to_string(ini_params.pool->nbytes_per_slab) 
    //             + ") is less than required slab size (="
    //             + to_string(mc.slab_size) + ")");
    // }
}

const int slow_pulsar_chunk::commit_chunk(std::shared_ptr<sp_chunk_header> header, std::shared_ptr<std::vector<uint8_t>> idat,
                      const ssize_t compressed_data_len, std::shared_ptr<std::vector<uint8_t>> mask,
                      std::shared_ptr<std::vector<float>> means, std::shared_ptr<std::vector<float>> vars)
{
	const ssize_t bytes_slab = this->memory_pool->nbytes_per_slab;
	const ssize_t size_head = header->get_header_size();
	// note that this must be explicitly provided as we can only predict
	// moments of the sample entropy

	// pull some constants from the chunk header
	const ssize_t nfreq = header->nfreq;
	const ssize_t nsamp = nfreq * header->ntime;

	const ssize_t size_i = compressed_data_len * sizeof(uint32_t); // word length
	const ssize_t size_m = (nsamp * sizeof(uint8_t)) / 8;
	const ssize_t size_freq = nfreq * sizeof(float);
	const ssize_t byte_size = size_head + size_i + size_m + 2 * size_freq + sizeof(long int);
	
	std::lock_guard<std::mutex> lg(this->slab_mutex);
	const ssize_t islab = this->islab;
	const ssize_t islab_post = islab + byte_size;

	if(byte_size > bytes_slab){
		return 2;
	}

	if(islab_post > bytes_slab){
		return 1;
	}

	// copy header
	header->copy((void*) &(this->memory_slab[islab]));
	// copy means
	std::memcpy((void*) &(this->memory_slab[islab + size_head]), (void*) &((*means)[0]), size_freq);
	// copy vars
	std::memcpy((void*) &(this->memory_slab[islab + size_head + size_freq]), 
							(void*) &((*vars)[0]), size_freq);
	// copy raw RFI mask
	std::memcpy((void*) &(this->memory_slab[islab + size_head + 2 * size_freq]), (void*) &((*mask)[0]), size_m);

	const long int comp_len = (long int) compressed_data_len;
	//copy compressed data len as long int (python int)
	std::memcpy((void*) &(this->memory_slab[islab + size_head + 2 * size_freq + size_m]), 
							(void*) &(comp_len), sizeof(long int));
	// copy encoded intensity data
	std::memcpy((void*) &(this->memory_slab[islab + size_head + 2 * size_freq + size_m + sizeof(long int)]),
							(void*) &((*idat)[0]), size_i);
	
	this->islab = islab_post;

	// This logic ensures that:
	//   - this->file_header.start is the start time of the first chunk in the file
	//   - this->file_header.end is the end time of the last chunk in the file
	//
	// (Previously in rf_pipelines::chime_slow_pulsar_writer::_process_chunk())

	const ssize_t fpga_nano = 2560;
	const double tnow = (header->fpga0 * fpga_nano + header->frame0_nano) * 1e-9;
	const double tend = tnow + (header->fpgaN * fpga_nano * 1e-9);

	if (this->file_header.start == 0.0)
	    this->file_header.start = tnow;

	this->file_header.end = tend;
	
	return 0;
}

std::shared_ptr<std::string> get_stem(const std::string& path){
	std::shared_ptr<std::string> rpath(new std::string(path));
	std::reverse(rpath->begin(), rpath->end());
	ssize_t strlen = rpath->length();
	std::shared_ptr<std::string> stem(new std::string(rpath->substr(rpath->find("/"), strlen)));
	std::reverse(stem->begin(), stem->end());
	std::cout << *stem << " " << path << std::endl;
	return stem;
}

void mkdir_recursive(std::shared_ptr<std::string> file_path, const std::string& delim,
				std::shared_ptr<std::string> dir_path = nullptr){
	if(!dir_path){
		dir_path = std::make_shared<std::string>("");
	}

	const ssize_t pos = file_path->find(delim);

	if(pos != -1){
		dir_path->append(delim);
		dir_path->append(file_path->substr(0, pos));
		file_path = std::make_shared<std::string>(file_path->substr(pos + 1, file_path->length()));
		// mkdir(dir_path->data(), 420); // permissions 644
		// TODO: fix output dir permissions issue
		// int mkres = mkdir(dir_path->data(), S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
		int mkres = mkdir(dir_path->data(), S_IRWXU | S_IRWXG | S_IRWXO);
		mkdir_recursive(file_path, delim, dir_path);
	}
}

// virtual override
void slow_pulsar_chunk::write_msgpack_file(const std::string &filename, bool compress,
                            uint8_t* buffer)
{
	// TODO: address ignored fields, consider renaming function in superclass

	std::shared_ptr<std::string> fnameptr = std::make_shared<std::string>(filename);


	// hard-code permissions (644)
	const int ofile = open(filename.data(), O_WRONLY | O_CREAT, 420);

	if(ofile == -1){
		// attempt to make the directory structure
		mkdir_recursive(fnameptr, "/");
		// std::cout << *fnameptr << std::endl;
		const int ofile2 = open(filename.data(), O_WRONLY | O_CREAT, 420);
		if(ofile2 == -1){
			throw std::runtime_error("slow_pulsar_chunk: failed to open new file to write chunk data out");
		}
	}

	std::lock_guard<std::mutex> lg(this->slab_mutex);
	const ssize_t nwrite = this->islab + 1 + this->file_header.get_header_size();

	ssize_t nwritten = this->file_header.write_to_file(ofile);
	nwritten += write(ofile, (void*) &(this->memory_slab[0]), this->islab + 1 );

	if(nwritten != nwrite){
		throw std::runtime_error("slow_pulsar_chunk: failed to write entire contents of chunk to file");
	}

	const int cfail = close(ofile);
	if(cfail == -1){
		throw std::runtime_error("slow_pulsar_chunk: failed to close output file");
	}
}

}
