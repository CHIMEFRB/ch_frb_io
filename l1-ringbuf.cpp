#include <unistd.h>
#include <iostream>
#include "ch_frb_io.hpp"
#include "ringbuf.hpp"
#include "l1-ringbuf.hpp"

using namespace ch_frb_io;
using namespace std;

namespace ch_frb_io {
#if 0
};  // pacify emacs c-mode!
#endif

std::ostream& operator<<(std::ostream& s, const assembled_chunk& ch) {
    s << "assembled_chunk(beam " << ch.beam_id << ", ichunk " << ch.ichunk << " at " << (void*)(&ch) << ")";
    return s;
}

static bool
assembled_chunk_overlaps_range(const shared_ptr<assembled_chunk> ch, 
                               uint64_t min_fpga_counts,
                               uint64_t max_fpga_counts) {
    if (min_fpga_counts == 0 && max_fpga_counts == 0)
        return true;
    uint64_t fpga0 = ch->fpgacounts_begin();
    uint64_t fpga1 = ch->fpgacounts_end();
    if ((max_fpga_counts && (fpga0 > max_fpga_counts)) ||
        (min_fpga_counts && (fpga1 < min_fpga_counts)))
        return false;
    return true;
}

AssembledChunkRingbuf::AssembledChunkRingbuf(int binlevel, L1Ringbuf* parent, int maxsize) :
    Ringbuf<assembled_chunk>(maxsize),
    _binlevel(binlevel),
    _parent(parent)
{}

AssembledChunkRingbuf::~AssembledChunkRingbuf() {}

void AssembledChunkRingbuf::dropping(shared_ptr<assembled_chunk> t) {
    _parent->dropping(_binlevel, t);
}

L1Ringbuf::L1Ringbuf(uint64_t beam_id, vector<int> ringbuf_n) :
    _beam_id(beam_id),
    _q(),
    _rb(),
    _dropped()
{
    if (ringbuf_n.size() == 0) {
        // ringbuffer sizes per binning level -- if specified,
        // overrides the number of levels and their sizes.  Otherwise,
        // take defaults
        for (int i=0; i<constants::assembled_ringbuf_nlevels; i++)
            ringbuf_n.push_back(constants::assembled_ringbuf_capacity);
    }
    _nbins = ringbuf_n.size();
    // Create the ring buffer objects for each time binning
    // (0 = native rate, 1 = binned by 2, ...)
    for (size_t i=0; i<_nbins; i++)
        _rb.push_back(shared_ptr<AssembledChunkRingbuf>
                      (new AssembledChunkRingbuf(i, this, ringbuf_n[i])));
    // Fill the "_dropped" array with empty shared_ptrs.
    for (size_t i=0; i<_nbins-1; i++)
        _dropped.push_back(shared_ptr<assembled_chunk>());
}

/*
 Tries to enqueue an assembled_chunk.  If no space can be allocated,
 returns an empty shared_ptr.  The ring buffer now assumes ownership
 of the assembled_chunk.
 */
shared_ptr<assembled_chunk> L1Ringbuf::push(assembled_chunk* ch) {
    shared_ptr<assembled_chunk> p = _rb[0]->push(ch);
    if (!p)
        return p;
    _q.push_back(p);
    return p;
}

/*
 Returns the next assembled_chunk for downstream processing.
 */
shared_ptr<assembled_chunk> L1Ringbuf::pop() {
    if (_q.empty())
        return shared_ptr<assembled_chunk>();
    shared_ptr<assembled_chunk> p = _q.front();
    _q.pop_front();
    return p;
}

/*
 Peeks at (returns) the next assembled_chunk for downstream processing,
 without popping it out of the queue.
 */
shared_ptr<assembled_chunk> L1Ringbuf::peek() {
    if (_q.empty())
        return shared_ptr<assembled_chunk>();
    shared_ptr<assembled_chunk> p = _q.front();
    return p;
}

int L1Ringbuf::n_ready() {
    return _q.size();
}

int L1Ringbuf::total_capacity() {
    int n = 0;
    for (auto it = _rb.begin(); it != _rb.end(); it++) {
        n += (*it)->maxsize();
    }
    return n;
}

int L1Ringbuf::total_size() {
    int n = 0;
    for (auto it = _rb.begin(); it != _rb.end(); it++) {
        n += (*it)->size();
    }
    return n;
}

// helper function used in the next function for visiting ringbuffer chunks
static void visit(uint64_t* fpga_min, uint64_t* fpga_max, shared_ptr<assembled_chunk> ch) {
    *fpga_min = std::min(*fpga_min, ch->fpgacounts_begin());
    *fpga_max = std::max(*fpga_max, ch->fpgacounts_end());
}

void L1Ringbuf::fpga_counts_range(uint64_t* min_fpga, uint64_t* max_fpga) {
    // set up visitor function args
    uint64_t mn, mx;
    mn = numeric_limits<uint64_t>::max();
    mx = 0;
    std::function<void(shared_ptr<assembled_chunk>)> func = std::bind(visit, &mn, &mx, std::placeholders::_1);

    for (auto it = _rb.begin(); it != _rb.end(); it++)
        (*it)->visit(func);

    if (min_fpga)
        *min_fpga = mn;
    if (max_fpga)
        *max_fpga = mx;
}

/*
 Prints a report of the assembled_chunks currently queued.
 */
void L1Ringbuf::print() {
    cout << "L1 ringbuf:" << endl;
    cout << "  downstream: [ ";
    for (auto it = _q.begin(); it != _q.end(); it++) {
        cout << (*it)->ichunk << " ";
    }
    cout << "]" << endl;
    for (size_t i=0; i<_nbins; i++) {
        vector<shared_ptr<assembled_chunk> > v = _rb[i]->snapshot(NULL);
        cout << "  binning " << i << ":  [ ";
        for (auto it = v.begin(); it != v.end(); it++) {
            cout << (*it)->ichunk << " ";
        }
        cout << "]" << endl;
        if (i < _nbins-1) {
            cout << "  dropped " << i << ":  ";
            if (_dropped[i])
                cout << _dropped[i]->ichunk << endl;
            else
                cout << "none" << endl;
        }
    }
}

void L1Ringbuf::retrieve(uint64_t min_fpga_counts, uint64_t max_fpga_counts,
                         vector<pair<shared_ptr<assembled_chunk>, uint64_t> >& chunks) {
    uint64_t where;
    // Check chunks queued for downstream processing
    where = L1RB_DOWNSTREAM;
    for (auto it = _q.begin(); it != _q.end(); it++)
        if (assembled_chunk_overlaps_range(*it, min_fpga_counts, max_fpga_counts))
            chunks.push_back(make_pair(*it, where));

    // Check ring buffers
    for (size_t i=0; i<_nbins; i++) {
        vector<shared_ptr<assembled_chunk> > ch;
        // uhhh, nice
        where = ((i == 0) ? L1RB_LEVEL1 : ((i == 1) ? L1RB_LEVEL2 : ((i == 2) ? L1RB_LEVEL3 : L1RB_LEVEL4)));
        _rb[i]->snapshot(ch, std::bind(assembled_chunk_overlaps_range, placeholders::_1, min_fpga_counts, max_fpga_counts));
        for (auto it = ch.begin(); it != ch.end(); it++)
            chunks.push_back(make_pair(*it, where));

        // Check the chunks that have been "dropped" but are waiting to be binned down.
        if ((i < _nbins-1) && (_dropped[i]) &&
            assembled_chunk_overlaps_range(_dropped[i], min_fpga_counts, max_fpga_counts)) {
            where = ((i == 0) ? L1RB_WAIT1 : ((i == 1) ? L1RB_WAIT2 : L1RB_WAIT3));
            chunks.push_back(make_pair(_dropped[i], where));
        }
    }
}

// Called from the AssembledChunkRingbuf objects when a chunk is
// about to be dropped from one binning level of the ringbuf.  If
// the chunk does not have a partner waiting (in _dropped), then
// it is saved in _dropped.  Otherwise, the two chunks are merged
// into one new chunk and added to the next binning level's
// ringbuf.
void L1Ringbuf::dropping(int binlevel, shared_ptr<assembled_chunk> ch) {
    // If it's the last binning level, just drop it
    if (binlevel >= (int)(_nbins-1))
        return;

    // Was there another chunk waiting to be binned down and merged with this one?
    if (_dropped[binlevel]) {

        // FIXME -- could call can_push() to the destination before
        // going to the effort of downsampling...

        // Allocate a new binned-down chunk.
        assembled_chunk* binned = assembled_chunk::downsample(NULL, _dropped[binlevel].get(), ch.get());
        // drop the old one we were saving
        _dropped[binlevel].reset();
        // drop the new one too
        ch.reset();
        // try to push onto _rb[level+1]
        if (!_rb[binlevel+1]->push(binned))
            // push failed; free
            delete binned;

    } else {
        // Keep this one until its partner arrives!
        _dropped[binlevel] = ch;
    }
}

} // namespace

