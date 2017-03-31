#ifndef _CH_FRB_L1_RINGBUF_HPP
#define _CH_FRB_L1_RINGBUF_HPP

#include <vector>
#include <deque>
#include <iostream>

#include "ch_frb_io.hpp"
#include "ringbuf.hpp"

namespace ch_frb_io {
#if 0
}; // pacify emacs c-mode
#endif

std::ostream& operator<<(std::ostream& s, const ch_frb_io::assembled_chunk& ch);

class L1Ringbuf;

class AssembledChunkRingbuf : public Ringbuf<ch_frb_io::assembled_chunk> {

public:
    AssembledChunkRingbuf(int binlevel, L1Ringbuf* parent, int maxsize);
    virtual ~AssembledChunkRingbuf();

protected:
    // my time-binning. level: 0 = original intensity stream; 1 =
    // binned x 2, 2 = binned x 4.
    int _binlevel;
    L1Ringbuf* _parent;

    // Called when the given frame *t* is being dropped off the buffer
    // to free up some space for a new frame.
    virtual void dropping(std::shared_ptr<ch_frb_io::assembled_chunk> t);
};

enum l1_ringbuf_level {
    L1RB_DOWNSTREAM = 1,
    L1RB_LEVEL1 = 2,
    L1RB_LEVEL2 = 4,
    L1RB_LEVEL3 = 8,
    L1RB_LEVEL4 = 0x10,
    L1RB_WAIT1 = 0x20,
    L1RB_WAIT2 = 0x40,
    L1RB_WAIT3 = 0x80,
    // queued for writing in the L1 RPC system
    L1RB_WRITEQUEUE = 0x100,
};

class L1Ringbuf {
    friend class AssembledChunkRingbuf;

public:
    L1Ringbuf(uint64_t beam_id, std::vector<int> ringbuf_n=std::vector<int>());

    /*
     Tries to enqueue an assembled_chunk.  If no space can be
     allocated, returns an empty shared_ptr.  The ring buffer now
     assumes ownership of the assembled_chunk.
     */
    std::shared_ptr<ch_frb_io::assembled_chunk> push(ch_frb_io::assembled_chunk* ch);

    /*
     Returns the next assembled_chunk for downstream processing.
     */
    std::shared_ptr<ch_frb_io::assembled_chunk> pop();

    /*
     Returns the next assembled_chunk for downstream processing,
     without removing it.
     */
    std::shared_ptr<ch_frb_io::assembled_chunk> peek();

    /*
     Returns the number of assembled_chunks queued for downstream processing.
     */
    int n_ready();

    /*
     Returns the total number (across all time binnings) of
     assembled_chunks this ring buffer is allowed to contain.
     */
    int total_capacity();

    /*
     Returns the total number (across all time binnings) of
     assembled_chunks currently being held in this ring buffer.
     */
    int total_size();

    /*
     Returns the oldest (smallest FPGA-counts) and newest (largest
     FPGA-counts) value of data in all levels of this ring buffer.
     */
    void fpga_counts_range(uint64_t* fpga_min, uint64_t* fpga_max);

    /*
     Prints a report of the assembled_chunks currently queued.
     */
    void print();

    /*
     Retrieves assembled_chunks that overlap the given range of FPGA
     counts values.  If min_fpga_counts is zero, no lower limit is
     applied; likewise for max_fpga_counts.
     */
    void retrieve(uint64_t min_fpga_counts, uint64_t max_fpga_counts,
                  std::vector<std::pair<std::shared_ptr<ch_frb_io::assembled_chunk>, uint64_t> > &chunks);

public:
    uint64_t _beam_id;
    
protected:
    // The number of binning levels
    size_t _nbins;

    // The queue for downstream
    std::deque<std::shared_ptr<ch_frb_io::assembled_chunk> > _q;

    // The ring buffers for each time-binning.  Length fixed at Nbins.
    std::vector<std::shared_ptr<AssembledChunkRingbuf> > _rb;

    // The assembled_chunks that have been dropped from the ring
    // buffers and are waiting for a pair to be time-downsampled.
    // Length fixed at Nbins-1.
    std::vector<std::shared_ptr<ch_frb_io::assembled_chunk> > _dropped;

    // Called from the AssembledChunkRingbuf objects when a chunk is
    // about to be dropped from one binning level of the ringbuf.  If
    // the chunk does not have a partner waiting (in _dropped), then
    // it is saved in _dropped.  Otherwise, the two chunks are merged
    // into one new chunk and added to the next binning level's
    // ringbuf.
    void dropping(int binlevel, std::shared_ptr<ch_frb_io::assembled_chunk> ch);
};


}  // namespace ch_frb_io

#endif
