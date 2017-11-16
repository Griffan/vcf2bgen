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

#ifndef BGEN_EXT_REFERENCE_IMPLEMENTATION_HPP
#define BGEN_EXT_REFERENCE_IMPLEMENTATION_HPP

#include "bgen.hpp"
#include <fstream>

/*
* This file contains a extended reference implementation of the BGEN file format
* specification described at:
* http://www.well.ox.ac.uk/~gav/bgen_format/bgen_format.html
*
*/


///////////////////////////////////////////////////////////////////////////////////////////
// INTERFACE
///////////////////////////////////////////////////////////////////////////////////////////

namespace genfile {
	namespace bgen {

	}
}	

///////////////////////////////////////////////////////////////////////////////////////////
// INTERFACE
///////////////////////////////////////////////////////////////////////////////////////////

namespace genfile {
	namespace bgen {
        namespace v12{
            namespace v12_HDS{
                template< typename Setter >
                void encode_probability_data(uint16_t numOfAlleles,
                                             const std::vector<std::vector<double> > &probs,
                                             Context const &context,
                                             Setter &setter
                ){
                    setter.initialise(context.number_of_samples, uint32_t(numOfAlleles));

                    int ploidy = 2;
                    bool phased = true;
                    bool missing = false;
                    for (uint32_t i = 0; i < context.number_of_samples; ++i) {
                        if (probs[i].empty())
                            missing = true;
                        else
                            missing = false;


                        uint32_t const valueCount//theoretically you need this much
                                = phased
                                  ? (ploidy * numOfAlleles)
                                  : genfile::bgen::impl::n_choose_k(uint32_t(ploidy + numOfAlleles - 1),
                                                                    uint32_t(numOfAlleles - 1));

                        uint32_t const storedValueCount =
                                valueCount - (phased ? ploidy : 1);//actually you store this much

                        if (setter.set_sample(i)) {
                            setter.set_number_of_entries(
                                    ploidy,
                                    valueCount,
                                    phased ? ePerPhasedHaplotypePerAllele : ePerUnorderedGenotype,
                                    eProbability
                            );
                            if (missing) {
                                // Consume dummy zero values, emit missing values.
                                for (uint32_t h = 0; h < valueCount; ++h) {
                                    setter.set_value(h, genfile::MissingValue());
                                }
                            } else {
                                // Consume values and interpret them.
                                double sum = 0.0;
                                uint32_t reportedValueCount = 0;
                                for (uint32_t h = 0; h < storedValueCount; ++h) {
                                    double const value = probs[i][h];
                                    setter.set_value(reportedValueCount++, value);
                                    sum += value;
                                    if ((phased &&
                                         ((h + 1) % (numOfAlleles - 1)) == 0)//last item doesn't need store or read
                                        || ((phased) && (h + 1) == storedValueCount)) {
                                        assert(sum <= 1.00000001);
                                        setter.set_value(reportedValueCount++, 1.0 - sum);
                                        sum = 0.0;
                                    }
                                }
                            }
                        }
                    }
                    call_finalise(setter);
                }
            }
        }




        void compress_probability_data(
                Context const &context,
                const std::vector<byte_t> & uncompressed_data,
                std::vector<byte_t> * compressed_data
        );
        template< typename Setter >
        void write_genotype_data_block(
                std::fstream &aStream,
                Setter &setter
        ) {
            uint32_t payload_size = setter.repr().second-setter.repr().first;
            aStream.write(const_cast<char*>(reinterpret_cast< const char *>(setter.repr().first)) , payload_size);
        }

        namespace impl {

            struct AlleleGetter {
                AlleleGetter( const std::vector<std::string>& alleles ):
                        m_alleles( alleles )
                {
                }

                std::string operator()( uint16_t i ) {
                    return m_alleles[i];
                }
                std::vector<std::string> m_alleles;
            } ;
        }

	}
}

#endif
