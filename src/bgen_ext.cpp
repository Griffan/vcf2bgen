/*The MIT License (MIT)
Copyright (c) 2017 Fan Zhang, Hyun Min Kang
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

#include "genfile/bgen/bgen_ext.hpp"
#include "fstream"
#ifdef CHAR_BIT
#if (CHAR_BIT != 8)
#error CHAR_BIT "Sorry, this implementation assumes 8-bit bytes. It won't work on your platform"
#endif
#endif

namespace genfile {
    namespace bgen {

        void compress_probability_data(
                Context const &context,
                const std::vector<byte_t> & uncompressed_data,
                std::vector<byte_t> * compressed_data
        ) {
            byte_t const *begin = &uncompressed_data[0];
            byte_t const *const end = &uncompressed_data[0] + uncompressed_data.size();
            std::vector<byte_t> * buffer = compressed_data;

            // compressed_data contains the (compressed or uncompressed) probability data.
            uint32_t const compressionType = (context.flags & bgen::e_CompressedSNPBlocks);
            if (compressionType != e_NoCompression) {//need compression
                if (compressionType == e_ZlibCompression) {
                    zlib_compress(begin, end, buffer, 8);
                } else if (compressionType == e_ZstdCompression) {
                    zstd_compress(begin, end, buffer, 8);
                }
            } else {
                // copy the data between buffers.
                compressed_data->resize(uncompressed_data.size()+4);
                std::copy(uncompressed_data.begin(), uncompressed_data.end(), compressed_data->begin()+4);
            }

            uint32_t payload_size = 0 ;
            int data_size = 0;
            if( (context.flags & e_Layout) == e_Layout2 || ((context.flags & e_CompressedSNPBlocks) != e_NoCompression ) ) {
                payload_size = 4 + compressed_data->size();
                write_little_endian_integer(&(buffer->at(0)), &(buffer->at(4)), payload_size);
                data_size = uncompressed_data.size();
                write_little_endian_integer(&(buffer->at(4)), &(buffer->at(8)), data_size);
            } else {
                payload_size = 6 * context.number_of_samples ;
                assert(payload_size == uncompressed_data.size());
                write_little_endian_integer(&(buffer->at(0)), &(buffer->at(4)), payload_size);
            }

        }
    }
}
