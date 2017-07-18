#include <cassert>
#include <iostream>
#include "ch_frb_io.hpp"

#include "assembled_chunk_msgpack.hpp"

#include "ch_frb_io_internals.hpp"


using namespace std;
using namespace ch_frb_io;

int main(int argc, char **argv)
{
    /*
     cout << "Warning: this will create a ~100 MB file in the current directory!\n"
     << "If this is OK press return.  If not, press control-C!\n"
     << "I AM WAITING, HUMAN: "
     << flush;
    string dummy;
    getline(cin, dummy);
     */

    uint64_t beam_id = 42;
    int nupfreq = 16;
    int nt_per_packet = 16;
    int fpga_counts_per_sample = 400;
    uint64_t ichunk = 1000;

    std::random_device rd;
    std::mt19937 rng(rd());

    unique_ptr<assembled_chunk> chunk = assembled_chunk::make(beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk);
    chunk->randomize(rng);

    /*
     const char *filename = "test_assembled_chunk.hdf5";
     chunk->write_hdf5_file(string(filename));
     */

    chunk->msgpack_bitshuffle = true;
    string fn = "test_assembled_chunk.msgpack";
    chunk->write_msgpack_file(fn);
    cout << "Wrote to " << fn << endl;

    shared_ptr<assembled_chunk> inchunk = assembled_chunk::read_msgpack_file(fn);
    //cout << "Read " << inchunk << endl;
    if (memcmp(inchunk->data, chunk->data, chunk->ndata)) {
        cout << "MISMATCH in data" << endl;
    }

    // Now fill with constant values and see that it compresses
    memset(chunk->data, 42, chunk->ndata);

    fn = "test_assembled_chunk_2.msgpack";
    chunk->write_msgpack_file(fn);
    cout << "Wrote to " << fn << endl;

    shared_ptr<assembled_chunk> inchunk2 = assembled_chunk::read_msgpack_file(fn);
    //cout << "Read " << inchunk2 << endl;
    if (memcmp(inchunk2->data, chunk->data, chunk->ndata)) {
        cout << "MISMATCH in data 2" << endl;
    }

    // Set the shared compression buffer and write compressed data.
    chunk->compression_buffer = shared_ptr<uint8_t>((uint8_t*)malloc(chunk->max_compressed_size()));
    chunk->msgpack_bitshuffle = true;
    fn = "test_assembled_chunk_3.msgpack";
    chunk->write_msgpack_file(fn);
    cout << "Wrote to " << fn << endl;
    
    shared_ptr<assembled_chunk> inchunk3 = assembled_chunk::read_msgpack_file(fn);
    if (memcmp(inchunk3->data, chunk->data, chunk->ndata)) {
        cout << "MISMATCH in data 3" << endl;
    }

    // Just write again, testing the buffer...?
    fn = "test_assembled_chunk_4.msgpack";
    chunk->write_msgpack_file(fn);
    cout << "Wrote to " << fn << endl;

    unique_ptr<assembled_chunk> uchunk2 = assembled_chunk::make(beam_id, nupfreq, nt_per_packet, fpga_counts_per_sample, ichunk+1);

    assembled_chunk* chunk2 = uchunk2.get();
    assembled_chunk* chunk1 = chunk.get();

    // Set up both chunk1 and chunk2 to contain ramps.

    for (int i=0; i<chunk1->nscales; i++) {
        chunk1->offsets[i] = uniform_rand(rng, 0, 10);
        chunk1->scales [i] = uniform_rand(rng, 1, 2);
    }
    for (int i=0; i<chunk2->nscales; i++) {
        chunk2->offsets[i] = uniform_rand(rng, 0, 10);
        chunk2->scales [i] = uniform_rand(rng, 1, 2);
    }

    for (int i=0; i<constants::nfreq_coarse_tot * chunk1->nupfreq; i++) {
        for (int j=0; j<constants::nt_per_assembled_chunk; j++) {
            int k = ((i / chunk1->nupfreq) * chunk1->nt_coarse +
                     (j / nt_per_packet));
            float offset = chunk1->offsets[k];
            float scale  = chunk1->scales [k];

            float fi = (float)i / (float)(constants::nfreq_coarse_tot * chunk1->nupfreq);
            float fj = (float)j / (float)constants::nt_per_assembled_chunk;
            float val = 20 + fi * 25 + fj * 50;

            chunk1->data[i * constants::nt_per_assembled_chunk + j] = (val - offset) / scale;

            offset = chunk2->offsets[k];
            scale  = chunk2->scales [k];

            val = 20 + fi * 25 + (fj+1) * 50;
            chunk2->data[i * constants::nt_per_assembled_chunk + j] = (val - offset) / scale;
        }
    }

    chunk1->write_msgpack_file("test-chunk1.msgpack");
    chunk2->write_msgpack_file("test-chunk2.msgpack");

    shared_ptr<assembled_chunk> chunk3 = shared_ptr<assembled_chunk>(assembled_chunk::downsample(NULL, chunk1, chunk2));
    chunk3->write_msgpack_file("test-chunk3.msgpack");

    float* intensity = (float*)malloc(chunk1->ndata * sizeof(float));
    float* weight    = (float*)malloc(chunk1->ndata * sizeof(float));

    chunk1->decode(intensity, weight, constants::nt_per_assembled_chunk);

    FILE* f = fopen("chunk1-decoded.raw", "w");
    fwrite(intensity, 1, chunk1->ndata * sizeof(float), f);
    fclose(f);

    chunk2->decode(intensity, weight, constants::nt_per_assembled_chunk);

    f = fopen("chunk2-decoded.raw", "w");
    fwrite(intensity, 1, chunk1->ndata * sizeof(float), f);
    fclose(f);

    chunk3->decode(intensity, weight, constants::nt_per_assembled_chunk);
    f = fopen("chunk3-decoded.raw", "w");
    fwrite(intensity, 1, chunk1->ndata * sizeof(float), f);
    fclose(f);

    free(intensity);
    free(weight);
    
    return 0;
}
