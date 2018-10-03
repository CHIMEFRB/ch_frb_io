#ifndef _ASSEMBLED_CHUNK_MSGPACK_HPP
#define _ASSEMBLED_CHUNK_MSGPACK_HPP

#include <vector>
#include <iostream>

#include <msgpack.hpp>

extern "C" {
  // UGH: c99
#define __STDC_VERSION__ 199901L
#include <bitshuffle.h>
}

#include <ch_frb_io.hpp>

/** Code for packing objects into msgpack mesages, and vice versa. **/

enum compression_type {
    comp_none = 0,
    comp_bitshuffle = 1
};

// Our own function for packing an assembled_chunk into a msgpack stream,
// including an optional buffer for compression.
template <typename Stream>
void pack_assembled_chunk(msgpack::packer<Stream>& o,
                          std::shared_ptr<ch_frb_io::assembled_chunk> const& ch,
                          bool compress=false,
                          uint8_t* buffer=NULL) {
    // pack member variables as an array.
    //std::cout << "Pack shared_ptr<assembled-chunk> into msgpack object..." << std::endl;
    uint8_t version = 1;
    // We are going to pack 18 items as a msgpack array (with mixed types)
    o.pack_array(18);
    // Item 0: header string
    o.pack("assembled_chunk in msgpack format");
    // Item 1: version number
    o.pack(version);

    uint8_t compression = (uint8_t)comp_none;
    int data_size = ch->ndata;

    // Create a shared pointer to the block of data to be written
    // (which defaults to this assembled_chunk's data, which is not to be deleted)
    std::shared_ptr<uint8_t> chdata(std::shared_ptr<uint8_t>(), ch->data);
    std::shared_ptr<uint8_t> data = chdata;
    if (compress) {
        compression = (uint8_t)comp_bitshuffle;
        if (buffer) {
            // We can use this buffer for compression
            data = std::shared_ptr<uint8_t>(std::shared_ptr<uint8_t>(), buffer);
        } else {
            // Try to allocate a temp buffer for the compressed data.
            // How big can the compressed data become?
            size_t maxsize = ch->max_compressed_size();
            std::cout << "bitshuffle: uncompressed size " << ch->ndata << ", max compressed size " << maxsize << std::endl;
            data = std::shared_ptr<uint8_t>((uint8_t*)malloc(maxsize));
            // unlikely...
            if (!data) {
                std::cout << "Failed to allocate a buffer to compress an assembled_chunk; writing uncompressed" << std::endl;
                compression = (uint8_t)comp_none;
                data = chdata;
                compress = false;
            }
        }
    }
    if (compress) {
        // Try compressing.  If the compressed size is not smaller than the original, write uncompressed instead.
        int64_t n = bshuf_compress_lz4(ch->data, data.get(), ch->ndata, 1, 0);
        if ((n < 0) || (n >= ch->ndata)) {
            if (n < 0)
                std::cout << "bitshuffle compression failed; writing uncompressed" << std::endl;
            else
                std::cout << "bitshuffle compression did not actually compress the data (" + std::to_string(n) + " vs orig " + std::to_string(ch->ndata) + "); writing uncompressed" << std::endl;
            data = chdata;
            compression = (uint8_t)comp_none;
            compress = false;
        } else {
            data_size = n;
            std::cout << "Bitshuffle compressed to " << n << std::endl;
        }
    }
    o.pack(compression);
    o.pack(data_size);

    o.pack(ch->beam_id);
    o.pack(ch->nupfreq);
    o.pack(ch->nt_per_packet);
    o.pack(ch->fpga_counts_per_sample);
    o.pack(ch->nt_coarse);
    o.pack(ch->nscales);
    o.pack(ch->ndata);
    o.pack(ch->fpga_begin);
    o.pack(ch->fpga_end - ch->fpga_begin);
    o.pack(ch->binning);
    // PACK FLOATS AS BINARY
    int nscalebytes = ch->nscales * sizeof(float);
    o.pack_bin(nscalebytes);
    o.pack_bin_body(reinterpret_cast<const char*>(ch->scales),
                    nscalebytes);
    o.pack_bin(nscalebytes);
    o.pack_bin_body(reinterpret_cast<const char*>(ch->offsets),
                    nscalebytes);
    o.pack_bin(data_size);
    o.pack_bin_body(reinterpret_cast<const char*>(data.get()), data_size);
    o.pack(ch->frame0_nano);
}

namespace msgpack {
MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS) {
namespace adaptor {

  // Unpack a msgpack object into an assembled_chunk.
template<>
struct convert<std::shared_ptr<ch_frb_io::assembled_chunk> > {
    msgpack::object const& operator()(msgpack::object const& o,
                                      std::shared_ptr<ch_frb_io::assembled_chunk>& ch) const {
        if (o.type != msgpack::type::ARRAY) throw msgpack::type_error();
        //std::cout << "convert msgpack object to shared_ptr<assembled_chunk>..." << std::endl;
        if (!((o.via.array.size == 17) ||
              (o.via.array.size == 18))) throw msgpack::type_error();
        msgpack::object* arr = o.via.array.ptr;

        std::string header         = arr[0].as<std::string>();
        uint8_t version            = arr[1].as<uint8_t>();
        if (version != 1)
            throw std::runtime_error("ch_frb_io: assembled_chunk msgpack version " + std::to_string(version) + ", expected 1");
        enum compression_type comp = (enum compression_type)arr[2].as<uint8_t>();
        if (!((comp == comp_none) || (comp == comp_bitshuffle)))
            throw std::runtime_error("ch_frb_io: assembled_chunk msgpack compression " + std::to_string(comp) + ", expected 0 or 1");
        int compressed_size        = arr[3].as<int>();
        int beam_id                = arr[4].as<int>();
        int nupfreq                = arr[5].as<int>();
        int nt_per_packet          = arr[6].as<int>();
        int fpga_counts_per_sample = arr[7].as<int>();
        int nt_coarse              = arr[8].as<int>();
        int nscales                = arr[9].as<int>();
        int ndata                  = arr[10].as<int>();
        uint64_t fpga0             = arr[11].as<uint64_t>();
        uint64_t fpgaN             = arr[12].as<uint64_t>();
        int binning                = arr[13].as<int>();
        int iarr = 14;

        uint64_t isample = fpga0 / (uint64_t)fpga_counts_per_sample;
        uint64_t ichunk = isample / ch_frb_io::constants::nt_per_assembled_chunk;

        uint64_t frame0_nano = 0;
        if (o.via.array.size == 18)
            frame0_nano = arr[17].as<uint64_t>();

	ch_frb_io::assembled_chunk::initializer ini_params;
	ini_params.beam_id = beam_id;
	ini_params.nupfreq = nupfreq;
	ini_params.nt_per_packet = nt_per_packet;
	ini_params.fpga_counts_per_sample = fpga_counts_per_sample;
	ini_params.binning = binning;
	ini_params.ichunk = ichunk;
        ini_params.frame0_nano = frame0_nano;

        ch = ch_frb_io::assembled_chunk::make(ini_params);

        if (ch->nt_coarse != nt_coarse)
            throw std::runtime_error("ch_frb_io: assembled_chunk msgpack nt_coarse mismatch");
	if (ch->nscales != nscales)
            throw std::runtime_error("ch_frb_io: assembled_chunk msgpack nscales mismatch");
	if (ch->ndata != ndata)
            throw std::runtime_error("ch_frb_io: assembled_chunk msgpack nscales mismatch");
        if (ch->fpga_begin != fpga0)
            throw std::runtime_error("ch_frb_io: assembled_chunk msgpack fpga_begin mismatch");
        if (ch->fpga_end != fpga0 + fpgaN)
            throw std::runtime_error("ch_frb_io: assembled_chunk msgpack fpga_end mismatch");

        if (arr[iarr + 0].type != msgpack::type::BIN) throw msgpack::type_error();
        if (arr[iarr + 1].type != msgpack::type::BIN) throw msgpack::type_error();
        if (arr[iarr + 2].type != msgpack::type::BIN) throw msgpack::type_error();

        uint nsdata = nscales * sizeof(float);
        if (arr[iarr + 0].via.bin.size != nsdata) throw msgpack::type_error();
        if (arr[iarr + 1].via.bin.size != nsdata) throw msgpack::type_error();

        memcpy(ch->scales,  arr[iarr + 0].via.bin.ptr, nsdata);
        memcpy(ch->offsets, arr[iarr + 1].via.bin.ptr, nsdata);

        if (comp == comp_none) {
            if (arr[iarr + 2].via.bin.size != (uint)ndata) throw msgpack::type_error();
            memcpy(ch->data, arr[iarr + 2].via.bin.ptr, ndata);
        } else if (comp == comp_bitshuffle) {
            if (arr[iarr + 2].via.bin.size != (uint)compressed_size) throw msgpack::type_error();
            //std::cout << "Bitshuffle: decompressing " << compressed_size << " to " << ch->ndata << std::endl;
            int64_t n = bshuf_decompress_lz4(reinterpret_cast<const void*>(arr[iarr + 2].via.bin.ptr), ch->data, ch->ndata, 1, 0);
            if (n != compressed_size)
                throw std::runtime_error("ch_frb_io: assembled_chunk msgpack bitshuffle decompression failure, code " + std::to_string(n));
        }

        return o;
    }
};

// Pack an assembled_chunk object into a msgpack stream.
template<>
struct pack<std::shared_ptr<ch_frb_io::assembled_chunk> > {
    template <typename Stream>
    packer<Stream>& operator()(msgpack::packer<Stream>& o, std::shared_ptr<ch_frb_io::assembled_chunk>  const& ch) const {
        pack_assembled_chunk(o, ch);
        return o;
    }
};

    /* Apparently not needed yet?
template <>
struct object_with_zone<std::shared_ptr<ch_frb_io::assembled_chunk> > {
    void operator()(msgpack::object::with_zone& o, std::shared_ptr<ch_frb_io::assembled_chunk>  const& v) const {
        o.type = type::ARRAY;
        std::cout << "Convert shared_ptr<assembled_chunk> into msgpack object_with_zone" << std::endl;
...
     */

} // namespace adaptor
} // MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS)
} // namespace msgpack

#endif // _ASSEMBLED_CHUNK_MSGPACK_HPP
