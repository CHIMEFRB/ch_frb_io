#include <iostream>
#include "ch_frb_io.hpp"
#include "ch_frb_io_internals.hpp"

using namespace ch_frb_io;
using namespace std;

int main() {
    int nupfreq = 4;
    int nt_per_packet = 16;
    int fpgacounts = 400;
    
    output_device::initializer out_ini;
    //memset(&out_ini, 0, sizeof(output_device::initializer));
    out_ini.device_name = "/tmp";
    out_ini.verbosity = 3;

    int nbytes_per_memory_slab = assembled_chunk::get_memory_slab_size(nupfreq, nt_per_packet);
    shared_ptr<memory_slab_pool> pool = make_shared<memory_slab_pool>(nbytes_per_memory_slab, 10, vector<int>(), 1);

    intensity_network_stream::initializer ini_params;
    ini_params.output_devices.push_back(output_device::make(out_ini));
    ini_params.memory_pool = pool;
    ini_params.nupfreq = nupfreq;
    ini_params.nt_per_packet = nt_per_packet;
    ini_params.fpga_counts_per_sample = fpgacounts;

    shared_ptr<assembled_chunk_ringbuf> assembler = make_shared<assembled_chunk_ringbuf>(ini_params, 1, 2);
    assembler->stream_to_files("/tmp/(CHUNK).msgpack", -1000);

    assembled_chunk::initializer chunk_ini;
    chunk_ini.pool = pool;
    for (int i=0; i<1; i++) {
      chunk_ini.beam_id = 1;
      chunk_ini.nupfreq = nupfreq;
      chunk_ini.nt_per_packet = nt_per_packet;
      chunk_ini.fpga_counts_per_sample = fpgacounts;
      chunk_ini.ichunk = i;
      chunk_ini.slab = pool->get_slab(true);
      unique_ptr<assembled_chunk> uch = assembled_chunk::make(chunk_ini);
      assembled_chunk* ch = uch.release();
      assembler->inject_assembled_chunk(ch);
    }
    cout << "Deleting assembler" << endl;
    assembler.reset();
    cout << "Deleted assembler" << endl;
    ini_params.memory_pool.reset();
    ini_params.output_devices.clear();
    chunk_ini.pool.reset();
    chunk_ini.slab.reset();
    pool.reset();

    usleep(3000000);
}
