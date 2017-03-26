#include "ch_frb_io.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


void pin_thread_to_cores(const vector<int> &core_list)
{
    if (core_list.size() == 0)
	return;

#ifdef __APPLE__
    cerr << "warning: pinning threads to cores is not implemented in osx\n";
#else
    int hwcores = std::thread::hardware_concurrency();
    pthread_t thread = pthread_self();

    cpu_set_t cs;
    CPU_ZERO(&cs);

    for (int core_id: core_list) {
	if ((core_id < 0) || (core_id >= hwcores))
	    throw runtime_error("pin_thread_to_cores: core_id=" + to_string(core_id) + " is out of range (hwcores=" + to_string(hwcores) + ")");
	CPU_SET(core_id, &cs);
    }

    int err = pthread_setaffinity_np(thread, sizeof(cs), &cs);
    if (err)
        throw runtime_error("pthread_setaffinity_np() failed");
#endif
}


}  // namespace ch_frb_io
