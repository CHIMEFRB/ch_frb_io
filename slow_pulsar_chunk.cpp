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

bool slow_pulsar_chunk::commit_chunk(sp_header& header, std::shared_ptr<std::vector<char>> dat)
{
	const ssize_t nbytes_slab = this->memory_pool->nbytes_per_slab;
	const ssize_t nhead = header.get_header_size();
	const ssize_t ndat = header.get_data_size();
	const ssize_t n = nhead + ndat;

	std::lock_guard<std::mutex> lg(this->slab_mutex);
	const ssize_t islab = this->islab;
	const ssize_t islab_post = islab + n;

	if(islab_post > nbytes_slab){
		return false;
	}

	std::memcpy((void*) &(this->memory_slab[islab]), (void*) &header, nhead);
	std::memcpy((void*) &(this->memory_slab[islab + nhead]), (void*) &((*dat)[0]), ndat);
	this->islab = islab_post;
	return true;
}

void slow_pulsar_chunk::write_msgpack_file(const std::string &filename, bool compress,
                            uint8_t* buffer)
{
	// no-op
}

ssize_t byte_ceil(const ssize_t bits, const ssize_t nbits)
{
    return bits/nbits + (bits % nbits > 0);
}

}