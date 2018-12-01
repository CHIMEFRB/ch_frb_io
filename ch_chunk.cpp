#include "ch_frb_io_internals.hpp"

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif

ch_chunk::~ch_chunk()
{
    if (memory_pool) {
		memory_pool->put_slab(memory_slab);
		memory_pool = std::shared_ptr<memory_slab_pool>();
    }
}

}