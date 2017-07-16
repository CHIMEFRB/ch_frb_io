#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif



memory_slab_pool::memory_slab_pool(ssize_t nbytes_per_slab_, ssize_t nslabs_, const vector<int> &allocation_cores, int verbosity_) :
    nbytes_per_slab(nbytes_per_slab_),
    nslabs(nslabs_),
    verbosity(verbosity_)
{
    double gb = pow(2.,-30.) * double(nbytes_per_slab) * double(nslabs);

    if (nbytes_per_slab <= 0)
	throw runtime_error("ch_frb_io: memory_slab_pool constructor expects nbytes_per_slab > 0");
    if (nslabs <= 0)
	throw runtime_error("ch_frb_io: memory_slab_pool constructor expects nslabs > 0");
    if (gb > 64.0)
	throw runtime_error("ch_frb_io: memory_slab_pool constructor: attempt to allocate > 64 GB, this is assumed unintentional");

    if (verbosity >= 1) {
	cout << "ch_frb_io: allocating " << gb << " GB memory pool";
	if (allocation_cores.size() > 0)
	    cout << ", cores=" << vstr(allocation_cores);
	cout << ", this may take a few seconds..." << endl;
    }

    std::thread t(std::bind(&memory_slab_pool::allocate, this, allocation_cores));
    t.join();

    if (verbosity >= 1)
	cout << "ch_frb_io: " << gb << " GB memory pool allocated" << endl;
}


memory_slab_pool::~memory_slab_pool()
{
    if (verbosity >= 1) {
	cout << "ch_frb_io: memory_slab_pool size_on_exit=" << curr_size << "/" << nslabs
	     <<", low_water_mark=" << low_water_mark << "/" << nslabs << endl;
    }
}


unique_ptr<uint8_t[]> memory_slab_pool::get_slab(bool zero, bool wait)
{
    ssize_t loc_size = 0;
    unique_ptr<uint8_t[]> ret;
    unique_lock<std::mutex> ulock(this->lock);

    for (;;) {
	if (curr_size > 0) {
	    ret.swap(slabs[curr_size-1]);
	    loc_size = --curr_size;
	    low_water_mark = min(low_water_mark, loc_size);
	    break;
	}

	if (!wait) {
	    ulock.unlock();
	    if (verbosity >= 2)
		cout << "ch_frb_io::memory_slab_pool::get_slab() FAILED" << endl;
	    return ret;
	}

	this->cv.wait(ulock);
    }

    ulock.unlock();

    if (verbosity >= 2)
       cout << "ch_frb_io::memory_slab_pool::get_slab(): " << loc_size << "/" << nslabs << endl;
    if (!ret)
	throw runtime_error("ch_frb_io: internal error: unexpected null pointer 'ret' in memory_slab_pool::get_slab()");
    if (zero)
	memset(ret.get(), 0, nbytes_per_slab);

    return ret;
}


void memory_slab_pool::put_slab(unique_ptr<uint8_t[]> &p)
{
    if (!p)
	throw runtime_error("ch_frb_io: internal error: unexpected null pointer 'p' in memory_slab_pool::put_slab()");

    unique_lock<std::mutex> ulock(this->lock);

    if (curr_size >= (int)slabs.size())
	throw runtime_error("ch_frb_io: internal error: buffer is full in memory_slab_pool::put_slab()");
    if (slabs[curr_size])
	throw runtime_error("ch_frb_io: internal error: unexpected null pointer 'slabs[curr_size]' in memory_slab_pool::put_slab()");
    
    slabs[curr_size].swap(p);
    ssize_t loc_size = ++curr_size;

    cv.notify_all();
    ulock.unlock();

    if (verbosity >= 2)
	cout << "ch_frb_io::memory_slab_pool::put_slab(): " << loc_size << "/" << nslabs << endl;
}


// Called as separate thread!
void memory_slab_pool::allocate(const vector<int> &allocation_cores)
{
    pin_thread_to_cores(allocation_cores);

    this->slabs.resize(nslabs);
    this->curr_size = nslabs;
    this->low_water_mark = nslabs;

    for (ssize_t i = 0; i < nslabs; i++) {
	uint8_t *p = aligned_alloc<uint8_t> (nbytes_per_slab);
	this->slabs[i] = unique_ptr<uint8_t[]> (p);
    }
}


}  // namespace ch_frb_io
