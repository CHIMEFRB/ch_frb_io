#ifndef _MSGPACK_BINARY_VECTOR_HPP
#define _MSGPACK_BINARY_VECTOR_HPP

#include <vector>
#include <iostream>

#include <msgpack.hpp>

namespace ch_frb_io {

template<typename T>
class msgpack_binary_vector : public std::vector<T>
{};

}

namespace msgpack {
MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS) {
namespace adaptor {

template<typename T>
struct convert<ch_frb_io::msgpack_binary_vector<T> > {
    msgpack::object const& operator()(msgpack::object const& o,
                                      ch_frb_io::msgpack_binary_vector<T>& v) const {
        //std::cout << "msgpack_binary_vector: type " << o.type << std::endl;
        if (o.type != msgpack::type::ARRAY)
            throw std::runtime_error("msgpack_binary_vector: expected type ARRAY");
        // Make sure array is big enough to check version
        //std::cout << "msgpack_binary_vector: array size " << o.via.array.size << std::endl;
        if (o.via.array.size != 3)
            throw std::runtime_error("msgpack_binary_vector: expected array size 3");
        msgpack::object* arr = o.via.array.ptr;
        uint8_t version = arr[0].as<uint8_t>();
        //std::cout << "version " << version << std::endl;
        if (version != 1)
            throw std::runtime_error("msgpack_binary_vector: expected version=1");
        size_t n = arr[1].as<size_t>();
        //std::cout << "msgpack_binary_vector: vector size " << n << std::endl; //", type " << arr[2].type << std::endl;
        v.resize(n);
        if (arr[2].type != msgpack::type::BIN)
            throw msgpack::type_error();
        //std::cout << "binary size " << arr[2].via.bin.size << " vs " << n << " x " <<
        //sizeof(T) << " = " << (n * sizeof(T)) << std::endl;
        if (arr[2].via.bin.size != n * sizeof(T))
            throw msgpack::type_error();
        memcpy(reinterpret_cast<void*>(v.data()), arr[2].via.bin.ptr, n * sizeof(T));
        //std::cout << "msgpack_binary_vector: returned vector size " << v.size() << std::endl;
        return o;
    }
};

template<typename T>
struct pack<ch_frb_io::msgpack_binary_vector<T> > {
    template <typename Stream>
    packer<Stream>& operator()(msgpack::packer<Stream>& o, ch_frb_io::msgpack_binary_vector<T> const& v) const {
        uint8_t version = 1;
        o.pack_array(3);
        o.pack(version);
        o.pack(v.size());
        o.pack_bin(v.size() * sizeof(T));
        o.pack_bin_body(reinterpret_cast<const char*>(v.data()));
        return o;
    }
};

} // namespace adaptor
} // MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS)
} // namespace msgpack

#endif



