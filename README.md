### ch_frb_io

C++/python library for CHIME-FRB file and network streams.

[![Build Status](https://travis-ci.org/CHIMEFRB/ch_frb_io.png?branch=master)](https://travis-ci.org/CHIMEFRB/ch_frb_io)
[![Coverage Status](https://coveralls.io/repos/github/CHIMEFRB/ch_frb_io/badge.svg?branch=master)](https://coveralls.io/github/CHIMEFRB/ch_frb_io?branch=master)

Currently there is no shared code between the C++ and python parts, so this is really two
independent libraries in the same git repository.  In fact the C++ and python parts use
different build systems, so you have to build and install them independently.  This is
something that should be fixed later!

CHIME FRB files have a similar format to CHIME cosmology data and can be read
using the tools in 'caput', in particular the caput.tod module
(https://github.com/radiocosmology/caput).
This package contains additional tools, such as simple file data streamers and a C++
interface.  

Network streams (read/write) are C++-only, but in the single-beam case, python wrappers
are available in the `rf_pipelines` repository.  The networking code is fairly mature
and well optimized/tested, but there are some missing features and loose ends noted
later in this README.

There isn't much documentation for the networking code right now, but the code is
pretty well commented, so reading the source code shouldn't be too painful.  The following
may also help:
  - There are some pdf slides in docs/
  - In the ch_frb_l1 github repo, there is a toy program which sends a network stream (ch-simulate-l0)
    and a toy program which receives a network stream (ch-frb-l1)
  - In the rf_pipelines repo, a python interface to the networking code is defined (but not
    as complete as the C++ interface).  In examples/example5_network, there is a pair of toy
    python pipelines which send and receive data over the network, making waterfall plots on
    both sides so that the input and output can be visually compared.
This code is mostly used as a library in other places (e.g. https://github.com/kmsmith137/ch_vdif_assembler
or https://github.com/kmsmith137/rf_pipelines), but there are a few standalone command-line utilities here:
```
ch-show-intensity-file   print summary statistics for an hdf5 intensity file
ch-plot-intensity-file   make a quick waterfall plot from an hdf5 intensity file (requires rf_pipelines)
decompress-chfrb-data    bitshuffle-decompress an hdf5 intensity file
```


### DEPENDENCIES

  1. libhdf5 (https://www.hdfgroup.org/HDF5/release/obtainsrc518.html)

     Note that this is a link to HDF5 v1.8.  I imagine v1.10 also works but haven't tested it yet.
     Assuming you want to use bitshuffle (see below), you'll need to install a very recent hdf5,
     so you'll need to compile one by hand instead of using 'yum'.  The following worked for me 
     (assuming non-root privs):
     ```
     cd hdf5-1.8.18    # after downloading the source code .tar.gz from the url above
     ./configure --prefix=$HOME
     make
     make install

     # The --prefix=$HOME flag installs libhdf5 in $HOME/lib, so I suggest adding this to .bashrc
     export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH
     ```

  2. The lz4 compression library.  In CentOS this is a one-liner: `sudo yum install lz4-devel`.

     In osx, this is also a one-liner: `brew install lz4`.

  3. The msgpack library.  In CentOS this is a one-liner: `sudo yum install msgpack-devel.x86_64`.
     
     In osx, this is also a one-liner: `brew install msgpack`.

     Note that you need a fairly new version (v?.??) -- the version in
     Ubuntu Trusty (14.04), v0.5.4, for example, is too old and fails
     because it is missing "msgpack/fbuffer.hpp".  Fortunately, you
     can download the source package from, eg,
     https://github.com/msgpack/msgpack-c/releases/download/cpp-2.1.0/msgpack-2.1.0.tar.gz
     and then extract it and add the "msgpack-2.1.0/include" into the
     include path -- no need to build it.

  4. Optional but recommended: bitshuffle (https://github.com/kiyo-masui/bitshuffle)
     You'll need this if you want to use bitshuffle-compressed files (note that CHIME pathfinder
     data is generally bitshuffle-compresed).

     A hint for installing bitshuffle:
     ```
     git clone https://github.com/kiyo-masui/bitshuffle.git
     cd bitshuffle/

     # The HDF5 library can dynamically load the bitshuffle plugin, i.e. you don't need
     # to link the bitshuffle library when you compile ch_frb_io, but you need to set this
     # environment variable to tell libhdf5 where to look.  Suggest adding this to .bashrc!

     export HDF5_PLUGIN_PATH=$HOME/lib/hdf5_plugins

     # If you have root privs and want to install "system-wide", omit the --user flag
     # The --h5plugin* flags will build/install the plugin needed to use bitshuffle from C++

     python setup.py install --user --h5plugin --h5plugin-dir=$HOME/lib/hdf5_plugins
     ```


### INSTALLATION (C++)

  - Create a file ./Makefile.local containing compiler flags, library locations, etc.
    The details are described in the Makefile.  There are some examples in the site/
    directory.  (You may be able to just symlink one of these examples to ./Makefile.local)

  - Compile and install with
    ```
    make all 
    make install
    ```

  - Here are some unit tests which you may or may not want to run:
    ```
    ./test-misc                       # no problem
    ./test-assembled-chunk            # no problem
    ./test-intensity-hdf5-file        # warning: creates 100MB temp file in current directory
    ./test-network-streams            # warning: takes ~1 hour to run, CPU-intensive
    ```

  - Note that running test-intensity-hdf5-file has the side effect of testing your
    bitshuffle installation.  If bitshuffle has _not_ been installed correctly, then
    you'll see the warning "couldn't load bitshuffle plugin, data will be written uncompressed".

  - You probably don't want to run test-network-streams unless you're actively
    working on the network code!

### INSTALLATION (PYTHON)

  - To build and install, do: `python setup.py install --user`
    (If you have root privs and want to install "system-wide", omit the --user flag)


### LOOSE ENDS IN NETWORKING CODE

  - We plan to implement RPC's which operate on the assembled_ringbufs, and support
    operations such as flushes to disk, and retrieval of ring-buffered data over the network.
    The first step here is doing some research to decide which RPC framework to use (maybe
    google RPC?)

  - We don't have thread-safe logging, so diagnostic messages from the various threads
    sometimes interleave each other and are unreadable.  (Implementing thread-safe
    logging is a general todo item for the whole CHIMEFRB backend, not just this repository!)

  - There is a proposal to bitshuffle-compress the packets, which may reduce bandwidth
    by ~20%, but this is currently unimplemented.

  - Open-ended item: there are lots of things that can go wrong in a realtime system,
    such as temporary network failures, and threads running slow so that ring buffers
    overfill.  We need to think carefully about different failure modes and figure out
    how best to handle them.

  - When the network stream is running, it maintains "event counts" for many types of
    events, such as packet drops, assembler hits/misses etc.  Right now we don't really
    do anything with this information but it's intended to be a starting point for some
    sort of RPC-driven dashboard which can give a visual summary of how the backend is
    performing.

    Related: we probably want to generalize the event counts (currently cumulative) to
    keep track of the event rate for some choice of timescale (say 10 sec).

    Another possible improvement: another way of making events more granular is to
    bin them by source IP address.

    One more!  It would be nice to keep track of the fraction of intensity samples which
    are masked (either 0x00 or 0xff bytes), or missing.  A natural place to put this counting 
    is assembled_chunk::add_packet(), but this is a little nontrivial since it involves
    modifying assembly language kernels, so I haven't done it yet.

  - Eventually it would be nice to do an end-to-end test of the FRB backend, by having
    a "simulator" node generate timestreams containing noise + FRB's, and sending them
    over the network.  However, the current simulation code is too slow to do this!

    The most important thing is to multithread the part of the code which simulates 
    the timestream (e.g. we could have one rf_pipeline per beam, with outputs combined
    into a ring buffer which feeds the intensity_network_ostream).  It may also help
    a little to write an assembly language packet encoding kernel (the correlator
    developers will probably do this eventually, so maybe we can just steal theirs).

  - There are two ring buffer data structures in the network code, a udp_packet_ringbuf
    and an assembled_chunk_ringbuf.  The udp_packet_ringbuf has been designed so that a 
    fixed pool of buffers is recycled throughout the lifetime of the ring buffer, whereas
    the assembled_chunk_ringbuf continually frees and allocates buffers.  It would be
    better to change the assembled_chunk_ringbuf to use a fixed buffer pool, in order to avoid
    the page-faulting cost of Linux malloc.

  - There are Linux-specific system calls sendmmsg(), recvmmsg() which send/receive
    multiple UDP packets, avoiding the overhead of one system call per packet.  This
    may help speed things up, or help reduce packet drops.

  - In the full CHIME backend, do we want to pin the network and/or assembler threads
    to specific cores?  Do we want to increase the scheduling priority of these threads?

  - It might be possible to further optimize the assembly-language kernels for packet assembly
    and decoding, using streaming writes.  See Chapter 7 ("optimizing cache usage") of the
    intel optimization manual.  We also might be able to improve performance a little using
    aligned loads/stores, but this would impose pointer alignment requirements on callers of
    assembled_chunk::decode().

  - The assembler should handle packets which arrive in an arbitrary order, but our
    unit test doesn't fully test this.  (It does permute the coarse frequencies, since
    this was easy to implement, but this isn't as general as an arbitrary packet permutation.)
    It would be great to strengthen the unit tests, by putting a flag in the network_ostream
    which randomly reorders packets prior to sending.

  - It would be great to switch from pthreads to C++11 threads, which are much nicer!
    However, there are some things to check.  First, we want to make sure that the C++11
    API supports low-level things like setting the scheduling affinity and priority of
    a thread.  Second, I'm not sure if C++11 threads and pthreads can interoperate.  If
    not, then we have to make the switch in many libraries at once!

  - Nuisance issue: if a chime_network_stream is constructed from python, then it doesn't
    respond to control-C (not sure if this is a 'ch_frb_io' loose end, or an 'rf_pipelines' 
    loose end).

  - Minor: implement an optimization to intensity_network_ostream which doesn't send
    a packet which is entirely masked (i.e. data array is all zeros)

  - It would be natural to include event counting / logging in the packet output stream,
    along the lines of what has already been implemented for the input stream.  For example,
    we could count
       - dropped packets
       - intensity samples which are assigned weight zero in the input
       - intensity samples which are assigned weight zero because they're below the wt_cutoff
       - intensity samples which are assigned weight zero because they're 5 sigma outliers

    This is low priority since we currently only use the output code for testing!
