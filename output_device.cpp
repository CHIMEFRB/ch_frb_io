#include "chlog.hpp"
#include "ch_frb_io_internals.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


// Static factory function 
shared_ptr<output_device> output_device::make(const output_device::initializer &ini_params)
{
    output_device *p = new output_device(ini_params);
    shared_ptr<output_device> ret(p);

    ret->output_thread = std::thread(std::bind(&output_device::io_thread_main, ret));
    return ret;
}


output_device::output_device(const output_device::initializer &ini_params_) :
    ini_params(ini_params_)
{
    // Note: no sanity-checking of ini_params needed here!
    // XXX 32M is overkill!  How much space should I use here?
    this->_buffer = unique_ptr<uint8_t[]> (new uint8_t[32 * 1024 * 1024]);
}


void output_device::io_thread_main()
{
    // Note: if 'io_thread_allowed_cores' is an empty list, no core-pinning will be done.
    pin_thread_to_cores(ini_params.io_thread_allowed_cores);

    // If io_thread_name is unspecified, assign a default thread name.
    string thread_name = ini_params.io_thread_name;
    if (thread_name.size() == 0)
	thread_name = "io_thread[" + ini_params.device_name + "]";

    chime_log_set_thread_name(thread_name);

    if (ini_params.verbosity >= 2)
	chlog("thread starting up");

    for (;;) {
	// Pull a write_chunk_request off the queue!
	shared_ptr<write_chunk_request> w = this->pop_write_request();
	if (!w)
	    break;

	shared_ptr<assembled_chunk> chunk = w->chunk;
	string error_message;
	string link_src;

	unique_lock<mutex> ulock(chunk->filename_mutex);

	// If this write request is a duplicate of a previous write request, return immediately.
	if (chunk->filename_set.count(w->filename) > 0) {
	    ulock.unlock();

	    if (ini_params.verbosity >= 3)
		chlog("write request '" + w->filename + "' is a duplicate, skipping...");

	    w->status_changed(true, true, "SUCCEEDED", "Duplicate filename " + w->filename + " already written");
	    continue;
	}
	
	// If this write request is a "pseudo-duplicate", i.e. the same chunk has been
	// previously written to a different filename on the same output_device, then
	// we'll make a hard link instead of a copy.

	auto p = chunk->filename_map.find(ini_params.device_name);	
	if (p != chunk->filename_map.end())
	    link_src = p->second;

	ulock.unlock();
        bool success = true;
	try {
	    if (file_exists(w->filename))
		throw runtime_error("Assembled chunk msgpack file to be written already exists: " + w->filename);
	    else if (link_src.size() == 0)
		chunk->write_msgpack_file(w->filename, false, this->_buffer.get());  // compress=false
	    else {
		if (ini_params.verbosity >= 3)
		    chlog("write request '" + w->filename + "' is a pseudo-duplicate, will hard-link from '" + link_src + "' instead of writing new copy");
		hard_link(link_src, w->filename);
	    }
	} catch (exception &e) {
	    // Note: we now include the exception text in the error_message.
            success = false;
	    error_message = "Write msgpack file '" + w->filename + "' failed: " + e.what();
	}

        if (success) {
            ulock.lock();
            chunk->filename_set.insert(w->filename);
            chunk->filename_map[ini_params.device_name] = w->filename;
            ulock.unlock();
        }

	if (error_message.size() > 0 && ini_params.verbosity >= 1)
	    chlog(error_message);
	else if (error_message.size() == 0 && ini_params.verbosity >= 3)
	    chlog("wrote " + w->filename);

	w->status_changed(true, success, success ? "SUCCEEDED" : "FAILED", error_message);
    }

    if (ini_params.verbosity >= 2)
	chlog("received NULL write_chunk_request; exiting.");
}


// Called by an "external" thread (assembler thread or RPC thread).
bool output_device::enqueue_write_request(const shared_ptr<write_chunk_request> &req)
{
    if (!req)
	throw runtime_error("ch_frb_io::output_device::enqueue_write_request(): req is null");
    if (!req->chunk)
	throw runtime_error("ch_frb_io::output_device::enqueue_write_request(): req->chunk is null");
    if (req->filename.size() == 0)
	throw runtime_error("ch_frb_io::output_device::enqueue_write_request(): req->filename is an empty string");
    if (!is_prefix(this->ini_params.device_name, req->filename))
	throw runtime_error("ch_frb_io::output_device::enqueue_write_request(): req->filename, device_name mismatch");
    
    if (req->need_rfi_mask && (req->chunk->nrfifreq <= 0)) {

	// I decided to throw an exception here, instead of returning false,
	// so that we get a "verbose" error message instead of an undiagnosed failure.
	// (The L1 server is responsible for not crashing, either by catching the exception
	// or by checking this error condition in advance.)
    
	throw runtime_error("ch_frb_io: enqueue_write_request() was called with need_rfi_mask=true,"
			    " but this server instance is not saving the RFI mask");
    }

    unique_lock<std::mutex> ulock(_lock);

    if (end_stream_called)
	return false;

    if (req->need_rfi_mask && !req->chunk->has_rfi_mask) {
        _awaiting_rfi.push_back(req);
        req->status_changed(false, true, "AWAITING_RFI", "Waiting for RFI mask to be computed");
    } else {
        _write_reqs.push(req);
        req->status_changed(false, true, "QUEUED", "Queued for writing");
        _cond.notify_all();
    }

    int n = _write_reqs.size();
    int na = _awaiting_rfi.size();
    ulock.unlock();
    
    if (ini_params.verbosity >= 3) {
	chlog("enqueued write request: filename " + req->filename
	      + ", beam " + to_string(req->chunk->beam_id) 
	      + ", chunk " + to_string(req->chunk->ichunk) 
	      + ", FPGA counts " + to_string(req->chunk->fpga_begin)
	      + ", write_queue_size=" + to_string(n)
              + ", awaiting_rfi_size=" + to_string(na));
    }

    return true;
}

int output_device::count_queued_write_requests() {
    unique_lock<std::mutex> ulock(_lock);
    int rtn = _write_reqs.size();
    ulock.unlock();
    return rtn;
}

// This gets called by the chime_mask_counter class in rf_pipelines
// to tell us that a chunk's RFI mask has been filled in.
void output_device::filled_rfi_mask(const std::shared_ptr<assembled_chunk> &chunk) {
    // FIXME? -- simplest approach -- tell i/o thread(s) waiting on the
    // queue not being empty that something may have happened!
    unique_lock<std::mutex> ulock(_lock);
    _cond.notify_all();
    ulock.unlock();
}

// Called internally by I/O thread.
shared_ptr<write_chunk_request> output_device::pop_write_request()
{
    unique_lock<std::mutex> ulock(_lock);

    for (;;) {
        // Check whether any chunks that were waiting for RFI masks to be
        // filled in have been.
        for (auto req=_awaiting_rfi.begin(); req!=_awaiting_rfi.end(); req++) {
            if ((*req)->chunk->has_rfi_mask) {
                cout << "Chunk " << (*req)->filename << " got its RFI mask!" << endl;
                _write_reqs.push((*req));
                auto toerase = req;
                req--;
                _awaiting_rfi.erase(toerase);
                (*req)->status_changed(false, true, "QUEUED", "RFI mask received; queued for writing");
            }
        }
	if (!_write_reqs.empty())
	    break;
	if (end_stream_called && _awaiting_rfi.empty())
	    return shared_ptr<write_chunk_request> ();
	_cond.wait(ulock);
    }

    shared_ptr<write_chunk_request> ret = _write_reqs.top();
    _write_reqs.pop();
    _cond.notify_all();

    int n = _write_reqs.size();
    ulock.unlock();

    if (ini_params.verbosity >= 3) {
	chlog("dequeued write request: filename " + ret->filename
	      + ", beam " + to_string(ret->chunk->beam_id) 
	      + ", chunk " + to_string(ret->chunk->ichunk) 
	      + ", FPGA counts " + to_string(ret->chunk->fpga_begin)
	      + ", write_queue_size=" + to_string(n));
    }

    return ret;
}


void output_device::end_stream(bool wait)
{
    unique_lock<std::mutex> ulock(_lock);

    // Note: multiple calls to end_stream() are OK.
    this->end_stream_called = true;

    // If 'wait' is false, cancel all pending writes and return.
    if (!wait) {
	while (!_write_reqs.empty())
	    _write_reqs.pop();
	_awaiting_rfi.clear();
	_cond.notify_all();
	return;
    }

    // If 'wait' is true, wait for pending writes to complete before returning.
    while (!_write_reqs.empty() || !_awaiting_rfi.empty())
	_cond.wait(ulock);
}


void output_device::join_thread()
{
    unique_lock<std::mutex> ulock(_lock);

    if (!end_stream_called)
	throw runtime_error("output_device::join_thread() called, with no preceding call to end_stream()");

    // Multiple calls to join_thread() are OK, but if this is not the first
    // time join_thread() has been called, we still need to check that the
    // thread has actually joined, and wait for it if not.

    if (join_thread_called) {
	while (!thread_joined)
	    _cond.wait(ulock);
	return;
    }

    // If we get here, this is the first time join_thread() has been called.
    // Set the 'join_thread_called' flag and drop the lock...
    join_thread_called = true;
    _cond.notify_all();
    ulock.unlock();

    // Now join the thread.  The sychronization logic above guarantees that no
    // other threads will call output_thread.join().  (Instead, they will wait
    // for this thread to set the 'thread_joined' flag.)
    output_thread.join();

    // Now re-acquire the lock and set the 'thread_joined' flag.
    ulock.lock();
    thread_joined = true;
    _cond.notify_all();
}


}  // namespace ch_frb_io
