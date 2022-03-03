#include "ch_frb_io_internals.hpp"

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif

ch_chunk::ch_chunk(const ch_chunk_initializer &ini_params) :
    fpga_counts_per_sample(ini_params.fpga_counts_per_sample),
    fpga_begin(ini_params.ichunk * constants::nt_per_assembled_chunk * ini_params.fpga_counts_per_sample),
    fpga_end((ini_params.ichunk + 1) * constants::nt_per_assembled_chunk * ini_params.fpga_counts_per_sample),
    frame0_nano(ini_params.frame0_nano),
    ichunk(ini_params.ichunk),
    beam_id(ini_params.beam_id)
{
    if ((fpga_counts_per_sample <= 0) || (fpga_counts_per_sample > constants::max_allowed_fpga_counts_per_sample))
	throw std::runtime_error("ch_chunk constructor: bad 'fpga_counts_per_sample' argument");
}

ch_chunk::~ch_chunk()
{
    if (memory_pool) {
		memory_pool->put_slab(memory_slab);
		memory_pool = std::shared_ptr<memory_slab_pool>();
    }
}

}
