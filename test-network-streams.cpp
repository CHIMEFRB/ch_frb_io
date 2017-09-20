// This is a long-running unit test which creates (input_stream, output_stream) pairs
// on the same machine, which communicate over the loopback interface 127.0.0.1.
//
// We check that the data captured and decoded by the input_stream is consistent
// with the data that was given to the output_stream to encode.  This is an end-to-end
// test of essentially all the network code, but it takes a long time to run!
//
// Note: this code isn't very well commented!

#include <cassert>
#include <algorithm>
#include "ch_frb_io_internals.hpp"

using namespace std;
using namespace ch_frb_io;


inline bool array_contains(const int *arr, int len, int x)
{
    for (int i = 0; i < len; i++)
	if (arr[i] == x)
	    return true;
    return false;
}

// An arbitrary but deterministic way of generating intensity values,
// so that the input and output streams can be compared.
inline float intval(int beam_id, int ifreq, int it)
{
    return sin(0.823*beam_id + 1.319*ifreq + 1.023*it);
}

// An arbitrary but deterministic way of generating weights.
inline float wtval(int beam_id, int ifreq, int it)
{
    return 0.6 + 0.4 * sin(1.328*beam_id + 2.382*ifreq + 0.883*it);
}


// -------------------------------------------------------------------------------------------------


struct unit_test_instance {
    static constexpr int maxbeams = 8;

    bool use_fast_kernels = false;
    int nbeams = 0;
    int nupfreq = 0;
    int nfreq_coarse_per_packet = 0;
    int nt_per_packet = 0;
    int nt_per_chunk = 0;
    int nt_tot = 0;
    int fpga_counts_per_sample = 0;
    uint64_t initial_t0 = 0;
    float wt_cutoff = 0.0;
    double target_gbps = 0.0;

    vector<int> recv_beam_ids;
    vector<int> send_beam_ids;
    vector<int> send_freq_ids;

    int send_istride = 0;
    int send_wstride = 0;
    int recv_istride = 0;
    int recv_wstride = 0;

    int nbytes_per_packet = 0;
    int npackets_per_chunk = 0;

    // not protected by lock
    pthread_t processing_threads[maxbeams];    
    shared_ptr<ch_frb_io::intensity_network_stream> istream;
    shared_ptr<ch_frb_io::intensity_network_ostream> ostream;

    pthread_mutex_t tpos_lock;
    pthread_cond_t cond_tpos_changed;
    uint64_t processing_tpos[maxbeams];

    unit_test_instance(std::mt19937 &rng, int irun, int nrun, double target_gbps);
    ~unit_test_instance();
};


unit_test_instance::unit_test_instance(std::mt19937 &rng, int irun, int nrun, double target_gbps_)
{
    const int nfreq_coarse_tot = ch_frb_io::constants::nfreq_coarse_tot;

#ifdef __AVX2__
    // In alternating iterations of the test, we choose parameters so that the "fast" kernels are used.
    this->use_fast_kernels = ((irun % 2) == 0);
#else
    this->use_fast_kernels = false;
#endif

    this->nbeams = randint(rng, 1, maxbeams+1);
    this->nupfreq = use_fast_kernels ? (2*randint(rng,1,9)) : randint(rng,1,17);
    this->nt_per_packet = use_fast_kernels ? 16 : (1 << randint(rng,0,5));

    // Now assign nfreq_coarse_per_packet, subject to packet size constraints.
    // The constants "c0" and "c1" are defined so that the packet size is c0 + c1 * nfreq_coarse_per_packet.
    int c0 = 24 + 2*nbeams;
    int c1 = 2 + 8*nbeams + nbeams*nupfreq*nt_per_packet;

    this->nfreq_coarse_per_packet = (ch_frb_io::constants::max_output_udp_packet_size - c0) / c1;
    this->nfreq_coarse_per_packet = min(nfreq_coarse_per_packet, ch_frb_io::constants::nfreq_coarse_tot);
    this->nfreq_coarse_per_packet = round_down_to_power_of_two(nfreq_coarse_per_packet);

    assert(nfreq_coarse_per_packet >= 4);
    this->nfreq_coarse_per_packet /= (1 << randint(rng,0,3));
    
    // Assign nt_per_chunk.  Each chunk should be no more than 512 samples.
    this->nt_per_chunk = nt_per_packet * randint(rng, 1, 512/nt_per_packet + 1);

    // Assign nt_tot.  We require <= 1024 chunks, and <= 1 GB total (summed over all beams).
    int packet_nbytes = intensity_packet::packet_size(nbeams, nfreq_coarse_per_packet, nupfreq, nt_per_packet);
    int chunk_nbytes = packet_nbytes * (nfreq_coarse_tot / nfreq_coarse_per_packet) * (nt_per_chunk / nt_per_packet);
    int max_nchunks = min(1024, (1<<30) / chunk_nbytes);
    this->nt_tot = nt_per_chunk * randint(rng, 1, max_nchunks+1);

    // We now require that initial_t0 is a multiple of nt_per_packet.
    this->initial_t0 = randint(rng, 0, 4097) * nt_per_packet;
    this->fpga_counts_per_sample = randint(rng, 1, 1025);
    this->wt_cutoff = uniform_rand(rng, 0.3, 0.7);
    this->target_gbps = target_gbps_;

    this->send_istride = randint(rng, nt_per_chunk, 2*nt_per_chunk+1);
    this->send_wstride = randint(rng, nt_per_chunk, 2*nt_per_chunk+1);
    this->recv_istride = randint(rng, constants::nt_per_assembled_chunk, 2 * constants::nt_per_assembled_chunk);
    this->recv_wstride = randint(rng, constants::nt_per_assembled_chunk, 2 * constants::nt_per_assembled_chunk);

#if 0
    // Sometimes it's convenient to debug a specific test case...
    this->use_fast_kernels = 0;
    this->nbeams = 6;
    this->nupfreq = 2;
    this->nfreq_coarse_per_packet = 16;
    this->nt_per_packet = 2;
    this->nt_per_chunk = 236;
    this->nt_tot = 8968;
    this->initial_t0 = 4282;
    this->fpga_counts_per_sample = 244;
    this->wt_cutoff = 0.63584;
    this->target_gbps = 0.1;
    this->send_istride = 445;
    this->send_wstride = 445;
    this->recv_istride = 1455;
    this->recv_wstride = 1455;
#endif

    // Clunky way of generating random beam_ids
    this->recv_beam_ids.resize(nbeams);
    for (int i = 0; i < nbeams; i++) {
	do {
	    recv_beam_ids[i] = randint(rng, 0, ch_frb_io::constants::max_allowed_beam_id);
	} while (array_contains(&recv_beam_ids[0], i, recv_beam_ids[i]));
    }

    // To slightly strengthen the test, we permute the receiver beams relative to sender
    this->send_beam_ids = recv_beam_ids;
    std::shuffle(send_beam_ids.begin(), send_beam_ids.end(), rng);

    // Randomly permute frequencies, just to strengthen the test
    this->send_freq_ids.resize(nfreq_coarse_tot);
    for (int i = 0; i < nfreq_coarse_tot; i++)
	send_freq_ids[i] = i;
    std::shuffle(send_freq_ids.begin(), send_freq_ids.end(), rng);

    this->nbytes_per_packet = intensity_packet::packet_size(nbeams, nfreq_coarse_per_packet, nupfreq, nt_per_packet);
    this->npackets_per_chunk = (nt_per_chunk / nt_per_packet) * (nfreq_coarse_tot / nfreq_coarse_per_packet);

    pthread_mutex_init(&this->tpos_lock, NULL);
    pthread_cond_init(&this->cond_tpos_changed, NULL);
    
    for (int ithread = 0; ithread < nbeams; ithread++)
	this->processing_tpos[ithread] = initial_t0;

    cout << "\nStarting test run " << irun << "/" << nrun << endl;

    cout << "    use_fast_kernels=" << use_fast_kernels << endl
	 << "    nbeams=" << nbeams << endl
	 << "    nupfreq=" << nupfreq << endl
	 << "    nfreq_coarse_per_packet=" << nfreq_coarse_per_packet << endl
	 << "    nt_per_packet=" << nt_per_packet << endl
	 << "    nt_per_chunk=" << nt_per_chunk << endl
	 << "    nt_tot=" << nt_tot << endl
	 << "    initial_t0=" << initial_t0 << endl
	 << "    fpga_counts_per_sample=" << fpga_counts_per_sample << endl
	 << "    wt_cutoff=" << wt_cutoff << endl
	 << "    target_gbps=" << target_gbps << endl
	 << "    send_istride=" << send_istride << endl
	 << "    send_wstride=" << send_wstride << endl
	 << "    recv_istride=" << recv_istride << endl
	 << "    recv_wstride=" << recv_wstride << endl
	 << "    nbytes_per_packet=" << nbytes_per_packet << endl
	 << "    npackets_per_chunk=" << npackets_per_chunk << endl;

#if 0
    // In a previous version of this test, we imposed the requirement that the
    // unassembled_ringbuf be large enough to contain the max allowed number of
    // in-flight packets.  Given the current level of optimization, this should
    // be overkill, so I removed the requirement but left it commented out.

    const int nt_assembler = 2 * ch_frb_io::constants::nt_per_assembled_chunk;

    // Worst-case storage requirements for unassembled ringbuf.
    int wc_nchunks = min(nt_assembler/nt_per_chunk + 1, nt_tot/nt_per_chunk);
    int wc_npackets = wc_nchunks * npackets_per_chunk;
    int wc_nbytes = wc_npackets * nbytes_per_packet;
    
    // Storage actually allocated
    int npackets_alloc = ch_frb_io::constants::unassembled_ringbuf_capacity * ch_frb_io::constants::max_unassembled_packets_per_list;
    int nbytes_alloc = ch_frb_io::constants::unassembled_ringbuf_capacity * ch_frb_io::constants::max_unassembled_nbytes_per_list;
    
    if ((npackets_alloc < wc_npackets) || (nbytes_alloc < wc_nbytes)) {
	cout << "    npackets_needed=" << wc_npackets << endl
	     << "    npackets_allocated=" << npackets_alloc << endl
	     << "    nbytes_needed=" << wc_nbytes << endl
	     << "    nbytes_alloc=" << nbytes_alloc << endl
	     << "Fatal: unassembled_packet_buf is underallocated" << endl;
	
	exit(1);
    }
#endif

    cout << endl;
}


unit_test_instance::~unit_test_instance()
{
    pthread_mutex_destroy(&tpos_lock);
    pthread_cond_destroy(&cond_tpos_changed);
}


// -------------------------------------------------------------------------------------------------


struct processing_thread_context {
    shared_ptr<unit_test_instance> tp;
    int ithread = 0;

    pthread_mutex_t lock;
    pthread_cond_t cond_running;
    bool is_running = false;

    processing_thread_context(const shared_ptr<unit_test_instance> &tp_, int ithread_)
	: tp(tp_), ithread(ithread_)
    {
	pthread_mutex_init(&lock, NULL);
	pthread_cond_init(&cond_running, NULL);
    }

    ~processing_thread_context()
    {
	pthread_mutex_destroy(&lock);
	pthread_cond_destroy(&cond_running);
    }
};


static void *processing_thread_main(void *opaque_arg)
{
    processing_thread_context *context = reinterpret_cast<processing_thread_context *> (opaque_arg);

    // Note: the processing thread startup logic works like this:
    //
    //   - parent thread puts a context struct on its stack, in spawn_processing_thread()
    //   - parent thread calls pthread_create() to spawn processing thread
    //   - parent thread blocks waiting for processing thread to set context->is_running
    //   - when parent thread unblocks, spawn_processing_thread() removes and the context struct becomes invalid
    //
    // Therefore, the processing thread is only allowed to access the context struct _before_
    // setting context->is_running to unblock the parent thread.  The first thing we do is
    // extract all members of the context struct so we don't need to access it again.

    shared_ptr<unit_test_instance> tp = context->tp;
    int ithread = context->ithread;
    
    // Now we can set context->is_running and unblock the parent thread.
    pthread_mutex_lock(&context->lock);
    context->is_running = true;
    pthread_cond_broadcast(&context->cond_running);
    pthread_mutex_unlock(&context->lock);

    int nfreq_f = ch_frb_io::constants::nfreq_coarse_tot * tp->nupfreq;
    vector<float> all_intensities(nfreq_f * tp->recv_istride, 0.0);
    vector<float> all_weights(nfreq_f * tp->recv_wstride, 0.0);

    double wt_cutoff = tp->wt_cutoff;
    int test_t0 = tp->initial_t0;
    int test_t1 = tp->initial_t0 + tp->nt_tot;
    int beam_id = tp->recv_beam_ids[ithread];

    uint64_t tpos = 0;
    bool tpos_initialized = false;

    for (;;) {
	auto chunk = tp->istream->get_assembled_chunk(ithread);

	if (!chunk)
	    break;

	chunk->decode(&all_intensities[0], &all_weights[0], tp->recv_istride, tp->recv_wstride);

	assert(chunk->nupfreq == tp->nupfreq);
	assert(chunk->nt_per_packet == tp->nt_per_packet);
	assert(chunk->fpga_counts_per_sample == tp->fpga_counts_per_sample);
	assert(chunk->beam_id == beam_id);

	if (tpos_initialized)
	    assert(chunk->isample == tpos);
	else
	    assert(chunk->isample <= tp->initial_t0);

	// tpos = expected isample in the next assembled_chunk.
	tpos = chunk->isample + ch_frb_io::constants::nt_per_assembled_chunk;
	tpos_initialized = true;

	pthread_mutex_lock(&tp->tpos_lock);
	tp->processing_tpos[ithread] = tpos;
	pthread_cond_broadcast(&tp->cond_tpos_changed);
	pthread_mutex_unlock(&tp->tpos_lock);

	int chunk_t0 = chunk->isample;

	for (int ifreq = 0; ifreq < nfreq_f; ifreq++) {
	    const float *int_row = &all_intensities[0] + ifreq * tp->recv_istride;
	    const float *wt_row = &all_weights[0] + ifreq * tp->recv_wstride;

	    for (int it = 0; it < ch_frb_io::constants::nt_per_assembled_chunk; it++) {
		// Out of range
		if ((it+chunk_t0 < test_t0) || (it+chunk_t0 >= test_t1)) {
		    assert(wt_row[it] == 0.0);
		    continue;
		}

		float ival = intval(beam_id, ifreq, it+chunk_t0);
		float wval = wtval(beam_id, ifreq, it+chunk_t0);

		// Check intensity
		if ((wt_row[it] > 0.0) && fabs(int_row[it] - ival) > 0.021) {
		    stringstream ss;
		    ss << "Test failure in weights array: beam_id=" << beam_id << ", ifreq=" << ifreq << ", it=" << (it+chunk_t0) << "\n"
		       << "   intval(...)=" << ival << ", intensity_received=" << int_row[it] << "\n";

		    cerr << ss.str();
		    exit(1);
		}

		// Check weights
		if ((wt_row[it] == 0.0) && (wval <= 1.00001 * wt_cutoff))
		    continue;
		if ((wt_row[it] == 1.0) && (wval >= 0.99999 * wt_cutoff))
		    continue;

		// If we get here, the weights check failed
		stringstream ss;
		ss << "Test failure in weights array: beam_id=" << beam_id << ", ifreq=" << ifreq << ", it=" << (it+chunk_t0) << "\n"
		   << "   wtval(...)=" << wval << ", wt_cutoff=" << wt_cutoff << ", wt_received=" << wt_row[it] << "\n";

		if ((wt_row[it] == 0.0) && (wval > wt_cutoff)) {
		    ss << "A possible reason for this failure is that some packet loss occurred.\n"
		       << "This unit test fails if there is any packet loss at all!\n"
		       << "This assumption is too extreme for the CHIME realtime environment, but it's useful in unit\n"
		       << "tests for sniffing out bugs.  However, it requires running the test at low throughput.\n"
		       << "By default, the tests run at 0.1 Gbps, but the -t command-line flag overrides the default.\n"
		       << "Suggest rerunning at lower throughput, to see whether the unit test failure is a \"true\"\n"
		       << "failure or an artifact of packet loss.\n";
		}

		cerr << ss.str();
		exit(1);
	    }
	}

	// more chunk processing will go here
    }

    assert(tpos_initialized);
    assert(tpos >= tp->initial_t0 + tp->nt_tot);

    return NULL;
}


static void spawn_processing_thread(const shared_ptr<unit_test_instance> &tp, int ithread)
{
    processing_thread_context context(tp, ithread);

    int err = pthread_create(&tp->processing_threads[ithread], NULL, processing_thread_main, &context);
    if (err)
	throw runtime_error(string("pthread_create() failed to create processing thread: ") + strerror(errno));

    pthread_mutex_lock(&context.lock);
    while (!context.is_running)
	pthread_cond_wait(&context.cond_running, &context.lock);
    pthread_mutex_unlock(&context.lock);
}


static void spawn_all_receive_threads(const shared_ptr<unit_test_instance> &tp)
{
    ch_frb_io::intensity_network_stream::initializer initializer;
    initializer.beam_ids = tp->recv_beam_ids;
    initializer.nupfreq = tp->nupfreq;
    initializer.nt_per_packet = tp->nt_per_packet;
    initializer.fpga_counts_per_sample = tp->fpga_counts_per_sample;
    initializer.force_reference_kernels = !tp->use_fast_kernels;
    initializer.force_fast_kernels = tp->use_fast_kernels;
    initializer.throw_exception_on_buffer_drop = true;
    initializer.throw_exception_on_assembler_miss = true;

    tp->istream = intensity_network_stream::make(initializer);
    
    for (int ithread = 0; ithread < tp->nbeams; ithread++)
	spawn_processing_thread(tp, ithread);
}


// -------------------------------------------------------------------------------------------------


static void send_data(const shared_ptr<unit_test_instance> &tp)
{
    intensity_network_ostream::initializer ini_params;
    ini_params.dstname = "127.0.0.1";
    ini_params.beam_ids = tp->send_beam_ids;
    ini_params.coarse_freq_ids = tp->send_freq_ids;
    ini_params.nupfreq = tp->nupfreq;
    ini_params.nt_per_chunk = tp->nt_per_chunk;
    ini_params.nfreq_coarse_per_packet = tp->nfreq_coarse_per_packet;
    ini_params.nt_per_packet = tp->nt_per_packet;
    ini_params.fpga_counts_per_sample = tp->fpga_counts_per_sample;
    ini_params.wt_cutoff = tp->wt_cutoff;
    ini_params.target_gbps = tp->target_gbps;

    // spawns network thread
    tp->ostream = intensity_network_ostream::make(ini_params);

    const int nt_assembler = 2 * ch_frb_io::constants::nt_per_assembled_chunk;
    const int nfreq_coarse_tot = ch_frb_io::constants::nfreq_coarse_tot;

    const int nbeams = tp->nbeams;
    const int nupfreq = tp->nupfreq;
    const int nt_chunk = tp->nt_per_chunk;
    const int nchunks = tp->nt_tot / tp->nt_per_chunk;
    const int istride = tp->send_istride;
    const int wstride = tp->send_wstride;
    const int si2 = nupfreq * istride;
    const int sw2 = nupfreq * wstride;
    const int si3 = nfreq_coarse_tot * si2;
    const int sw3 = nfreq_coarse_tot * sw2;

    vector<float> intensity(nbeams * si3, 0.0);
    vector<float> weights(nbeams * sw3, 0.0);

    for (int ichunk = 0; ichunk < nchunks; ichunk++) {
	int chunk_t0 = tp->initial_t0 + ichunk * nt_chunk;   // start of chunk

	for (int ibeam = 0; ibeam < nbeams; ibeam++) {
	    int beam_id = tp->send_beam_ids[ibeam];

	    for (int ifreq_coarse = 0; ifreq_coarse < nfreq_coarse_tot; ifreq_coarse++) {
		int coarse_freq_id = tp->send_freq_ids[ifreq_coarse];

		for (int iupfreq = 0; iupfreq < nupfreq; iupfreq++) {
		    int ifreq_logical = coarse_freq_id * nupfreq + iupfreq;
		    
		    float *i_row = &intensity[ibeam*si3 + ifreq_coarse*si2 + iupfreq*istride];
		    float *w_row = &weights[ibeam*sw3 + ifreq_coarse*sw2 + iupfreq*wstride];

		    for (int it = 0; it < nt_chunk; it++) {
			i_row[it] = intval(beam_id, ifreq_logical, chunk_t0 + it);
			w_row[it] = wtval(beam_id, ifreq_logical, chunk_t0 + it);
		    }
		}
	    }
	}

	// Wait for processing threads if necessary.
	// Note that for some choices of unit_test_instance parameters, this can test the timeout logic.
	pthread_mutex_lock(&tp->tpos_lock);
	for (int i = 0; i < nbeams; i++) {
	    while (tp->processing_tpos[i] + nt_assembler < uint64_t(chunk_t0))
		pthread_cond_wait(&tp->cond_tpos_changed, &tp->tpos_lock);
	}
	pthread_mutex_unlock(&tp->tpos_lock);

	uint64_t fpga_count = (tp->initial_t0 + ichunk * nt_chunk) * tp->fpga_counts_per_sample;
	tp->ostream->send_chunk(&intensity[0], istride, &weights[0], wstride, fpga_count);
	cout << "sent chunk " << ichunk << "/" << nchunks << endl;
    }

    // joins network thread
    tp->ostream->end_stream(true);
}


// -------------------------------------------------------------------------------------------------


// This unit test fails if there is any packet loss at all!
// This assumption is too extreme for the CHIME realtime environment, but it's useful in unit
// tests for sniffing out bugs.  However, it requires running the test at low throughput
// (~0.1 Gbps).  The -t command-line flag can override this default.
//
// Note that on linux, an alternative to a low throughput would be to use a local unix socket
// with SOCK_SEQPACKET.  I didn't implement this because it's not supported in osx, and I like
// to run the unit tests on my laptop.

static const double default_target_gbps = 0.1;


static void usage(const char *extra=nullptr)
{
    cerr << "usage: ./test-network-streams [-t TARGET_GBPS]\n"
	 << "   if -t is unspecified, then target_gbps defaults to " << default_target_gbps << "\n";

    if (extra != nullptr)
	cerr << extra << "\n";

    exit(2);
}


int main(int argc, char **argv)
{
    const int nrun = 100;
    double target_gbps = default_target_gbps;

    // Low-budget command-line parsing.
    int iarg = 1;
    while (iarg < argc) {
	if ((iarg < argc-1) && !strcmp(argv[iarg], "-t")) {
	    if (!lexical_cast(argv[iarg+1], target_gbps))
		usage();
	    if (target_gbps <= 0.0)
		usage("Fatal: expected target_gbps > 0");
	    iarg += 2;
	}
	else
	    usage();
    }

#if 1
    std::random_device rd;
    std::mt19937 rng(rd());
#else
    // Sometimes it's convenient to use the same seed every time, for debugging.
    std::mt19937 rng;
#endif

    cout << "Warning: this cpu-intensive, multithreaded test will take over your machine for about an hour!\n"
	 << "If this is OK press return.  If not, press control-C!\n"
	 << "I AM WAITING, HUMAN: "
	 << flush;

    string dummy;
    getline(cin, dummy);

    for (int irun = 0; irun < nrun; irun++) {
	auto tp = make_shared<unit_test_instance> (rng, irun, nrun, target_gbps);

	spawn_all_receive_threads(tp);
	tp->istream->start_stream();

	send_data(tp);	
	tp->istream->join_threads();

	for (int ibeam = 0; ibeam < tp->nbeams; ibeam++) {
	    int err = pthread_join(tp->processing_threads[ibeam], NULL);
	    if (err)
		throw runtime_error("pthread_join() failed");	    
	}

	vector<int64_t> counts = tp->istream->get_event_counts();

	typedef ch_frb_io::intensity_network_stream::event_type ev_type;
	assert(counts[ev_type::packet_bad] == 0);
	assert(counts[ev_type::packet_bad] == 0);
	assert(counts[ev_type::packet_dropped] == 0);
	assert(counts[ev_type::beam_id_mismatch] == 0);
	assert(counts[ev_type::stream_mismatch] == 0);
	assert(counts[ev_type::assembler_miss] == 0);
	assert(counts[ev_type::assembled_chunk_dropped] == 0);

	int expected_npackets = (tp->nt_tot / tp->nt_per_chunk) * tp->npackets_per_chunk;
	assert(counts[ev_type::packet_received] - counts[ev_type::packet_end_of_stream] == expected_npackets);
    }

    cout << "\n    ****  network test passed!! ****\n\n";

    return 0;
}
