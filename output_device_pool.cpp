#include "ch_frb_io_internals.hpp"
#include "chlog.hpp"

using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif


output_device_pool::output_device_pool(const vector<shared_ptr<output_device>> &streams_) :
    streams(streams_)
{
    int nstreams = streams.size();
    this->device_names.resize(nstreams);

    for (int i = 0; i < nstreams; i++) {
	if (!streams[i])
	    throw runtime_error("ch_frb_io: empty pointer in output_device_pool constructor");
	device_names[i] = streams[i]->ini_params.device_name;
    }

    // Sanity check: if one device_name is a prefix of another, throw an exception
    for (int i = 0; i < nstreams; i++) {
	for (int j = 0; j < nstreams; j++) {
	    if ((i != j) && is_prefix(device_names[i], device_names[j])) {
		throw runtime_error("ch_frb_io: output_device_pool constructor: device name '" 
				    + device_names[i] + "' is a prefix of device name '" 
				    + device_names[j] + "'");
	    }
	}
    }
}


bool output_device_pool::enqueue_write_request(shared_ptr<write_chunk_request> req)
{
 //    if (req->need_wait && (req->chunk->nrfifreq <= 0)) {

	// // I decided to throw an exception here, instead of returning false,
	// // so that we get a "verbose" error message instead of an undiagnosed failure.
	// // (The L1 server is responsible for not crashing, either by catching the exception
	// // or by checking this error condition in advance.)
    
	// throw runtime_error("ch_frb_io: enqueue_write_request() was called with need_rfi_mask=true,"
	// 		    " but this server instance is not saving the RFI mask");
 //    }
	  
    if (req->filename.size() == 0)
	return false;

    for (unsigned int i = 0; i < device_names.size(); i++) {
	if (is_prefix(device_names[i], req->filename))
	    return streams[i]->enqueue_write_request(req);
    }

    if (device_names.size() == 0)
	throw runtime_error("ch_frb_io: enqueue_write_request() was called, but no output_devices were registered");

    chlog("Request to write to file \"" << req->filename << "\" did not match any of the output_devices registered; request dropped.");

    return false;
}


void output_device_pool::end_streams(bool wait)
{
    for (const auto &s: streams)
	s->end_stream(wait);
}


void output_device_pool::join_threads()
{
    for (const auto &s: streams)
	s->join_thread();
}


}  // namespace ch_frb_io
