#include "ch_frb_io_internals.hpp"

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif

std::shared_ptr<slow_pulsar_chunk> slow_pulsar_chunk::make_slow_pulsar_chunk(ch_chunk::initializer& ini_params)
{
	std::shared_ptr<memory_slab_pool> pool;
	if(!ini_params.pool)
		throw std::runtime_error("slow_pulsar_chunk: memory_slab_pool required to construct new chunk");
	
	if(!ini_params.slab){
		ini_params.slab = ini_params.pool->get_slab(true);

		if (!ini_params.slab)
	    	throw std::runtime_error("**** Too much memory pressure for this poor L1 node to survive!  Blowing up now... ****");
	}

	return std::shared_ptr<slow_pulsar_chunk>(new slow_pulsar_chunk(ini_params));
}

void slow_pulsar_chunk::write_msgpack_file(const std::string &filename, bool compress,
                            uint8_t* buffer)
{
	//no-op
}
}