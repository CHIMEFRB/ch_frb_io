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
    this->_write_reqs.reserve(32768);
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

	if (ini_params.verbosity >= 3) {
	    chlog("dequeued write request: filename " + w->filename
		  + ", beam " + to_string(w->chunk->beam_id) 
		  + ", chunk " + to_string(w->chunk->ichunk) 
		  + ", FPGA counts " + to_string(w->chunk->fpga_begin));
	}

	string error_message;

	try {
	    w->chunk->msgpack_bitshuffle = true;
	    w->chunk->write_msgpack_file(w->filename);
	} catch (exception &e) {
	    // Note: we now include the exception text in the error_message.
	    error_message = "Write msgpack file '" + w->filename + "' failed: " + e.what();
	}

	if (error_message.size() > 0 && ini_params.verbosity >= 1)
	    chlog(error_message);
	else if (ini_params.verbosity >= 3)
	    chlog("wrote " + w->filename);

	// The write_chunk_request::next pointers are used to maintain a linked list
	// of "redundant" write requests for the same (chunk, filename) pair.  Maintaining
	// this list is only necessary in order to ensure that all write_callbacks get
	// called, whether or not the request was redundant.  Here, we traverse the list 
	// and do the write_callbacks.

	for (shared_ptr<write_chunk_request> ww = w; ww; ww = ww->next)
	    ww->write_callback(error_message);
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
    if (req->next)
	throw runtime_error("ch_frb_io::output_device::enqueue_write_request(): req->next is non-null");
    if (!is_prefix(this->ini_params.device_name, req->filename))
	throw runtime_error("ch_frb_io::output_device::enqueue_write_request(): req->filename, device_name mismatch");

    unique_lock<std::mutex> ulock(_lock);

    if (end_stream_called)
	return false;

    // Two requests are considered "redundant" if they have the same (chunk, filename) pair.
    // For each pending request, we maintain ..
    // If the new write_request is redundant with an existing request, we 

    for (const auto &qreq: _write_reqs) {
	if ((qreq->chunk == req->chunk) && (qreq->filename == req->filename)) {
	    qreq->priority = max(req->priority, qreq->priority);
	    req->next = qreq->next;
	    qreq->next = req;
	    
	    ulock.unlock();

	    if (ini_params.verbosity >= 3)
		chlog("write_request for filename '" + req->filename + "' absorbed into existing write request");

	    return true;
	}
    }

    // No redundant request was found, just add the new request to the list.
    _write_reqs.push_back(req);
    _cond.notify_all();

    ulock.unlock();
    
    if (ini_params.verbosity >= 3) {
	chlog("enqueued write request: filename " + req->filename
	      + ", beam " + to_string(req->chunk->beam_id) 
	      + ", chunk " + to_string(req->chunk->ichunk) 
	      + ", FPGA counts " + to_string(req->chunk->fpga_begin));
    }

    return true;
}


// Called internally by I/O thread.
shared_ptr<write_chunk_request> output_device::pop_write_request()
{
    unique_lock<std::mutex> ulock(_lock);

    for (;;) {
	if (_write_reqs.size() > 0)
	    break;
	if (end_stream_called)
	    return shared_ptr<write_chunk_request> ();
	_cond.wait(ulock);
    }

    int n = _write_reqs.size();
    int max_pri = _write_reqs[0]->priority;
    int max_ix = 0;

    for (int i = 0; i < n; i++) {
	if (_write_reqs[i]->priority > max_pri) {
	    max_pri = _write_reqs[i]->priority;
	    max_ix = i;
	}
    }

    shared_ptr<write_chunk_request> ret = _write_reqs[max_ix];
    _write_reqs[max_ix] = _write_reqs[n-1];
    _write_reqs.resize(n-1);
    _cond.notify_all();
    return ret;
}


void output_device::end_stream(bool wait)
{
    unique_lock<std::mutex> ulock(_lock);

    // Note: multiple calls to end_stream() are OK.
    this->end_stream_called = true;

    // If 'wait' is false, cancel all pending writes and return.
    if (!wait) {
	_write_reqs.clear();
	_cond.notify_all();
	return;
    }

    // If 'wait' is true, wait for pending writes to complete before returning.
    while (_write_reqs.size() > 0)
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
