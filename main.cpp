#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <memory>
#include "libVcf/VcfFile.h"
#include "libVcf/libStatGenLite/Error.h"
#include "src/genfile/bgen/bgen_ext.hpp"

/****typedef****/
typedef std::unordered_map<std::string, int> HashStringInt;
typedef std::unordered_set<std::string> SetString;
typedef std::set<int> SetInt;

/****typedef****/



// BgenEncoder is a thin wrapper around the core functions in genfile/bgen/bgen.hpp.
// This class tracks file state and handles passing the right callbacks.
struct BgenEncoder {

    BgenEncoder( std::string const& filename ):
            m_filename( filename ),
            m_state( e_NotOpen ),
            m_have_sample_ids( false ),
            m_stream(filename)
    {
        if(!m_stream.is_open())
        {
            throw std::invalid_argument( filename ) ;
        }
        m_state = e_Open ;
        //Set context information.TODO: initialize m_context
        m_context.flags|=genfile::bgen::e_SampleIdentifiers;
        m_context.flags|=genfile::bgen::e_Layout2;
        m_context.flags|=genfile::bgen::e_ZstdCompression;
        m_context.magic="bgen";

        // Write the offset, header, and sample IDs if present.
        m_offset = 4;
        m_offset += m_context.header_size();
        genfile::bgen::write_offset( m_stream, m_offset );//m_offset is the size of header block
        genfile::bgen::write_header_block( m_stream, m_context ) ;

    }

    std::ostream& summarise( std::ostream& o ) const {
        o << "BgenEncoder: bgen file ("
          << ( m_context.flags & genfile::bgen::e_Layout2 ? "v1.2 layout" : "v1.1 layout" )
          << ", "
          << ( m_context.flags & genfile::bgen::e_CompressedSNPBlocks ? "compressed" : "uncompressed" ) << ")"
          << " with "
          << m_context.number_of_samples << " " << ( m_have_sample_ids ? "named" : "anonymous" ) << " samples and "
          << m_context.number_of_variants << " variants.\n" ;
        return o ;
    }

    void update_header()
    {
        m_stream.seekg(4);
        genfile::bgen::write_header_block( m_stream, m_context );
    }

    std::size_t number_of_samples() const {
        return m_context.number_of_samples ;
    }

    // Fill in the sample IDs 
    // (If there are no sample IDs in the file, we report a dummy identifier).
    void set_sample_ids( const std::vector<std::string> &  sampleVec, const SetInt & retainedSample/*index in original vcf*/) {
        if( m_have_sample_ids ) {
            for( std::size_t i = 0; i < sampleVec.size();++i ) {
                if(retainedSample.find(i)!= retainedSample.end())
                    m_sample_ids[i]=sampleVec[i];//assign name to relative index in new bgen file
            }
        } else {
            for( std::size_t i = 0; i < m_context.number_of_samples; ++i ) {
                m_sample_ids[i]= "(unknown_sample_" + std::to_string( i+1 ) + ")";
            }
        }
    }

    bool write_sample_identifier(const std::vector<std::string> &  sampleVec, const SetInt & retainedSample)
    {
        m_context.number_of_samples=retainedSample.size();
        m_have_sample_ids = true;
        m_sample_ids=std::vector<std::string>(retainedSample.size(),"retainedSample");
        set_sample_ids(sampleVec, retainedSample);
        genfile::bgen::write_sample_identifier_block(m_stream, m_context,m_sample_ids) ;
        // We keep track of state (though it's not really needed for this implementation.)
        m_state = e_ReadyForVariant ;
        return true;
    }
    // Attempt to read identifying information about a variant from the bgen file, returning
    // it in the given fields.
    // If this method returns true, data was successfully read, and it should be safe to call read_probs()
    // or ignore_probs().
    // If this method returns false, data was not successfully read indicating the end of the file.
    bool write_variant(
            const std::string& chromosome,
            const uint32_t position,
            const std::string& rsid,
            const std::vector< std::string >& alleles,
            const std::string& SNPID
    ) {
        assert( m_state == e_ReadyForVariant ) ;
        genfile::bgen::impl::AlleleGetter alleleGetter(alleles);
        std::vector< genfile::byte_t> snp_identifier_buffer;
        if( genfile::bgen::write_snp_identifying_data(
                &snp_identifier_buffer, m_context,
                SNPID, rsid, chromosome, position,
                alleles.size(),
                alleleGetter))
        {
            m_stream.write(reinterpret_cast<char*>(&snp_identifier_buffer[0]), sizeof(genfile::byte_t)*snp_identifier_buffer.size());
            m_context.number_of_variants++;
            m_state = e_ReadyForProbs ;
            return true ;
        } else {
            return false ;
        }
    }
    // Write genotype probability data for the SNP
    void write_probs( std::vector< std::vector< double > > probs, uint16_t numOfAlleles ) {
        assert( m_state == e_ReadyForProbs ) ;

        genfile::bgen::GenotypeDataBlockWriter dataBlockWriter(&m_buffer1,&m_buffer2,m_context, 16);
        genfile::bgen::v12::v12_HDS::encode_probability_data(numOfAlleles, probs, m_context, dataBlockWriter);
        //genfile::bgen::compress_probability_data(m_context,m_buffer1, &m_buffer2);
        genfile::bgen::write_genotype_data_block(m_stream,dataBlockWriter);
        m_state = e_ReadyForVariant ;
    }

    void release_resource()
    {
        m_stream.close();
    }

private:
    std::string const m_filename ;
    std::fstream m_stream ;

    // bgen::Context object holds information from the header block,
    // including bgen flags
    genfile::bgen::Context m_context ;

    // offset byte from top of bgen file.
    uint32_t m_offset ;

    // We keep track of our state in the file.
    // Not strictly necessary for this implentation but makes it clear that
    // calls must be read_variant() followed by read_probs() (or ignore_probs())
    // repeatedly.
    enum State { e_NotOpen = 0, e_Open = 1, e_ReadyForVariant = 2, e_ReadyForProbs = 3, eComplete = 4 } ;
    State m_state ;

    // If the BGEN file contains samples ids, they will be read here.
    bool m_have_sample_ids ;
    std::vector< std::string > m_sample_ids ;

    // Buffers, these are used as working space by bgen implementation.
    std::vector< genfile::byte_t > m_buffer1, m_buffer2 ;
} ;



using namespace libVcf;

std::unordered_set<std::string> sampleWhiteList;//if none, retain all the input samples
std::unordered_set<std::string> markerWhiteList;

int main(int argc, char ** argv) {
    std::string inputFilePath(argv[1]);//quick and dirty
    std::string outputFilePath(argv[2]);//quick and dirty

    printf("Load VCF\n\n");
    /****write header*****/
    BgenEncoder encoder(outputFilePath);
    /****write header done****/
    try {
        VcfFile *pVcf = new VcfFile;
        pVcf->bSiteOnly = false;
        pVcf->bParseGenotypes = false;
        pVcf->bParseDosages = false;
        pVcf->bParseValues = true;
        pVcf->openForRead(inputFilePath.c_str());

        bool keepAllSample=true;
        bool keepAllMarker=true;
        if(!sampleWhiteList.empty()) keepAllSample=false;
        if(!markerWhiteList.empty()) keepAllMarker=false;

        // check the sanity of data
        if (pVcf->getSampleCount() == 0) {
            throw VcfFileException("No individual genotype information exist in the input VCF file %s",
                                   inputFilePath.c_str());
        }

        int nSamples = pVcf->getSampleCount();
        HashStringInt sampleOrderInVCF; // key: sampleID, value: original order (0 based);
        std::vector<std::string> sampleVec;
        int person = 0;
        for (int i = 0; i < nSamples; i++) {
            sampleOrderInVCF.insert(std::make_pair(pVcf->vpVcfInds[i]->sIndID, person));
            sampleVec.push_back(pVcf->vpVcfInds[i]->sIndID.c_str());
            person++;
        }

        SetInt retainedSampleIndex;

        if(not keepAllSample)//only keep those in while list
        {
            for (const auto &k:sampleWhiteList) {//first N samples in ped, not necessary the unphased sample
                int index = sampleOrderInVCF[k.c_str()];
                if (index != -1) {//find index of this sample in current vcf
                    retainedSampleIndex.insert(index);//map index in current vcf to index in ped file
                }
            }
        }
        else
            for (const auto &s : sampleOrderInVCF)
                retainedSampleIndex.insert(s.second);

        /***write samples****/
        encoder.write_sample_identifier(sampleVec,retainedSampleIndex);
        /****write samples done***/
        int markerNumber = 0;
        VcfMarker *pMarker = nullptr;//new VcfMarker;
        String markerName;

        while (pVcf->iterateMarker()) {//for each marker

            pMarker = pVcf->getLastMarker();
            markerName.printf("%s:%d:%s", pMarker->sChrom.c_str(), pMarker->nPos,pMarker->asAlts[0].c_str());//assume tri-alleles was divided into two lines

            if (not keepAllMarker and markerWhiteList.find(markerName.c_str())==markerWhiteList.end())
                continue;
            //int AFidx = pMarker->asInfoKeys.Find("AF");

            int GTidx = pMarker->asFormatKeys.Find("GT");
            int DSidx = pMarker->asFormatKeys.Find("DS");
            int HDSidx = pMarker->asFormatKeys.Find("HDS");
            int GPidx = pMarker->asFormatKeys.Find("GP");


            if (DSidx < 0 || GPidx < 0 || GTidx <0 ||HDSidx <0) {
                throw VcfFileException("Cannot recognize GT, DS, GP or HDS key in FORMAT field");
            }

            int formatLength = pMarker->asFormatKeys.Length();
//            int idx11 = 0, idx12 = 1, idx22 = 2;

            StringArray tmpStrArray;
//            int genoindex = markerindex * 3;

//            long allele1(-1), allele2(-1);
//            double ds(0.);
//            double hap1(0.), hap2(0.);
//            double gp11(-1), gp12(-1), gp22(-1);

            // Output variants
            std::string chromosome = pMarker->sChrom.c_str();
            uint32_t position = pMarker->nPos;
            std::string rsid = pMarker->sID.c_str();
            std::vector< std::string > alleles;
            alleles.push_back(pMarker->sRef.c_str());
            for (int k = 0; k < pMarker->asAlts.Length(); ++k) {
                alleles.push_back(pMarker->asAlts[k].c_str());
            }

            /*******write variant identifier*****/
            encoder.write_variant(chromosome,position,rsid,alleles,
                                  std::string(chromosome+"_"+std::to_string(position)+"_"+alleles[0]));
            /*******write variant identifier done*****/

            std::vector< std::vector< double > > probs ;
            std::vector<double> gp;
            for (int i = 0; i < nSamples; i++)//for each individual in current vcf
            {
                if (keepAllSample or retainedSampleIndex.find(i) != retainedSampleIndex.end()) {//not in index mapping relations, which is impossible

//                    ds=pMarker->asSampleValues[DSidx + i * formatLength].AsDouble();

//                    tmpStrArray.ReplaceTokens(pMarker->asSampleValues[GPidx + i * formatLength], ",");
//
//                    gp11= tmpStrArray[idx11].AsDouble();
//                    gp12= tmpStrArray[idx12].AsDouble();
//                    gp22= tmpStrArray[idx22].AsDouble();


//                    tmpStrArray.ReplaceTokens(pMarker->asSampleValues[GTidx + i * formatLength], "|/");
//                    allele1=tmpStrArray[0].AsInteger();
//                    allele2=tmpStrArray[1].AsInteger();

                    tmpStrArray.ReplaceTokens(pMarker->asSampleValues[HDSidx + i * formatLength], ",");
                    if(tmpStrArray[0]=="." or tmpStrArray[1] == ".")
                        probs.push_back(std::vector<double>());
                    else {
                        gp.push_back(tmpStrArray[0].AsDouble());
                        gp.push_back(tmpStrArray[1].AsDouble());
                        probs.push_back(gp);
                    }
                    gp.clear();

                }
            }
            encoder.write_probs(probs,2);
            markerNumber++;
        }
        delete pVcf;
    }
    catch (VcfFileException e) {
        error(e.what());
    }
    encoder.summarise(std::cerr);
    encoder.update_header();
    encoder.release_resource();
    std::cout << "Hello, World!" << std::endl;
    return 0;
}