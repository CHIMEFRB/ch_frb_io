#include<fcntl.h>
#include<sys/types.h>
#include<sys/stat.h>

#include "ch_frb_io_internals.hpp"

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

bool slow_pulsar_chunk::commit_chunk(sp_chunk_header& header, std::shared_ptr<std::vector<uint32_t>> idat,
                      const ssize_t compressed_data_len, std::shared_ptr<std::vector<uint8_t>> mask,
                      std::shared_ptr<std::vector<float>> means, std::shared_ptr<std::vector<float>> vars)
{
	const ssize_t bytes_slab = this->memory_pool->nbytes_per_slab;
	const ssize_t size_head = header.get_header_size();
	// note that this must be explicitly provided as we can only predict
	// moments of the "sample entropy"
	const ssize_t size_i = compressed_data_len * sizeof(uint32_t);
	const ssize_t size_m = mask->size() * sizeof(uint8_t);
	const ssize_t size_freq = means->size() * sizeof(float);
	const ssize_t byte_size = size_head + size_i + size_m + 2 * size_freq;

	std::lock_guard<std::mutex> lg(this->slab_mutex);
	const ssize_t islab = this->islab;
	const ssize_t islab_post = islab + byte_size;

	if(islab_post > bytes_slab){
		return false;
	}

	// copy header
	std::memcpy((void*) &(this->memory_slab[islab]), (void*) &header, size_head);
	// copy encoded intensity data
	std::memcpy((void*) &(this->memory_slab[islab + size_head]), (void*) &((*idat)[0]), size_i);
	// copy raw RFI mask
	std::memcpy((void*) &(this->memory_slab[islab + size_head + size_i]), (void*) &((*mask)[0]), size_m);
	// copy means
	std::memcpy((void*) &(this->memory_slab[islab + size_head + size_i + size_m]), (void*) &((*means)[0]), size_freq);
	// copy vars
	std::memcpy((void*) &(this->memory_slab[islab + size_head + size_i + size_m + size_freq]), 
							(void*) &((*vars)[0]), size_freq);
	
	this->islab = islab_post;
	return true;
}

// virtual override
void slow_pulsar_chunk::write_msgpack_file(const std::string &filename, bool compress,
                            uint8_t* buffer)
{
	// TODO: address ignored fields, consider renaming function in superclass
	
	// hard-code permissions
	const int ofile = open(filename.data(), O_WRONLY | O_CREAT, 644);
	if(ofile == -1){
		throw std::runtime_error("slow_pulsar_chunk: failed to open new file to write chunk data out");
	}

	std::lock_guard<std::mutex> lg(this->slab_mutex);
	const ssize_t nwrite = this->islab + 1 + this->file_header.get_header_size();

	ssize_t nwritten = write(ofile, (void*) &(this->file_header), this->file_header.get_header_size());
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