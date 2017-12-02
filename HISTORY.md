-v9:
  - 'beam_ids' argument to intensity_network_stream::stream_files()

- v8:
  - from Dustin: per-sender packet count history, retrievable by RPC

- v7:
  - implement "sleep hack", a temporary kludge that will go away soon!

- v6:

  - this large commit is part of the "2017 Mega Merge" affecting many parts of the CHIMEFRB pipeline
  - minor API changes, now that we have independent strides for the intensity/weights ring buffers.

- v5: lots of under-the-hood updates.

  - fast assembly language downsampling kernels (needed to make telescoping ring buffer run in real time).
  - simplified telescoping ring buffer logic.
  - use assembled_chunk free pool, instead of dynamic malloc/free.
  - more work on write_requests, assembled_chunk write path.

- v4: large update from Dustin.

  - Backward incompatible: lz4 and msgpack are now dependencies.
  - Import bitshuffle code (which adds the "lz4" library as a dependency).
  - Add assembled_chunk.downsample() method.
  - Add hdf5 output for assembled_chunks.
  - Add msgpack i/o for assembled_chunks.
  - Add multi-level, memory-footprint-limited ring buffer of assembled_chunks.
  - Add calls to retrieve statistics and assembled_chunks to support RPC calls.
  - Add debugging calls to inject assembled_chunks into the assembler.

- v3: large update containing an intial implementation of chimefrb networking code.

  - Backward incompatible: the C++ source files have been moved from cpp/ to the
    top-level directory, for consistency with other chimefrb repos.  This is a 
    "backward incompatible" change only in the sense that your Makefile.local
    should now be in the toplevel directory, not in cpp/.

  - Backwards incompatible: compiler flags must now include -pthread

  - Networking code has been written.  The most important classes are:

      - intensity_network_stream, which receives a packet stream over the 
        network, decodes it and assembles a sequence of regular arrays for 
        each beam.

      - intensity_network_ostream, which receives a sequence of regular arrays,
        encodes and packetizes the data, and sends it over the network.

- v2: implement intensity hdf5 file writing (backwards-compatible update)
  - writing intensity hdf5 files from C++
  - proper unit test of the C++ part
  - new utility ch-plot-intensity-file which makes a quick waterfall plot from an intensity hdf5 file.

- v1: initial version, contains classes for reading intensity hdf5 files and not much else.
