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
    out_ini.device_name = "/test";
    out_ini.verbosity = 3;

    intensity_network_stream::initializer ini_params;
    ini_params.output_devices.push_back(output_device::make(out_ini));
    ini_params.nupfreq = nupfreq;
    ini_params.nt_per_packet = nt_per_packet;
    ini_params.fpga_counts_per_sample = fpgacounts;
    
    shared_ptr<assembled_chunk_ringbuf> assembler = make_shared<assembled_chunk_ringbuf>(ini_params, 1, 2);
    assembler->stream_to_files("/test/x.msgpack", 0);

    assembled_chunk::initializer chunk_ini;
    chunk_ini.beam_id = 1;
    chunk_ini.nupfreq = nupfreq;
    chunk_ini.nt_per_packet = nt_per_packet;
    chunk_ini.fpga_counts_per_sample = fpgacounts;
    chunk_ini.ichunk = 1;
    unique_ptr<assembled_chunk> uch = assembled_chunk::make(chunk_ini);
    assembled_chunk* ch = uch.release();
    assembler->inject_assembled_chunk(ch);
    cout << "Deleting assembler" << endl;
    assembler.reset();
    cout << "Deleted assembler" << endl;

    usleep(3000000);
}
