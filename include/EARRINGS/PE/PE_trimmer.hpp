#pragma once
#include <EARRINGS/common.hpp>
#include <EARRINGS/PE/task_processor.hpp>
#include <Biovoltron/format/fastq.hpp>
#include <Biovoltron/format/fasta_peat.hpp>
#include <Nucleona/parallel/asio_pool.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <tuple>
#include <iostream>

using namespace EARRINGS;
using namespace biovoltron::format;
using namespace nucleona::parallel;
using BitStr = biovoltron::vector<
               biovoltron::char_type, 
               simdpp::aligned_allocator<simdpp::uint8<16>, 16>
               >;

namespace EARRINGS
{
void PE_trim()
{
    using BIO_filtering_istream = boost::iostreams::filtering_istream;
    using BIO_filtering_ostream = boost::iostreams::filtering_ostream;

    constexpr size_t const_thread_num = 2;
    std::tuple<float, float, float> trimmer_param = std::make_tuple(match_rate, seq_cmp_rate, adapter_cmp_rate);
    std::vector<std::string> default_adapter{DEFAULT_ADAPTER1, DEFAULT_ADAPTER2};
    
    #define INIT_TASK_PROCESSOR task_processor(thread_num\
                            , block_size\
                            , record_line\
                            , ifs_name\
                            , ofs_name\
                            , trimmer_param\
                            , loc_tail\
                            , default_adapter\
                            , DETECT_N_READS\
                            , is_sensitive)

    if (is_fastq)
    {
        if (is_gz_input)
        {
            if (is_gz_output)
            {
                TaskProcessor<
                    FASTQ, BitStr, 
                    BIO_filtering_istream, BIO_filtering_ostream
                > INIT_TASK_PROCESSOR;
                task_processor.process();
            }
            else
            {
                TaskProcessor<
                    FASTQ, BitStr, 
                    BIO_filtering_istream, std::ofstream
                > INIT_TASK_PROCESSOR;
                task_processor.process();
            }
        }
        else
        {
            if (is_gz_output)
            {
                TaskProcessor<
                    FASTQ, BitStr, 
                    std::ifstream, BIO_filtering_ostream
                > INIT_TASK_PROCESSOR;
                task_processor.process();
            }
            else
            {
                TaskProcessor<
                    FASTQ, BitStr, 
                    std::ifstream, std::ofstream
                > INIT_TASK_PROCESSOR;
                task_processor.process();
            }
        }
    }
    else
    {
        if (is_gz_input)
        {
            if (is_gz_output)
            {
                TaskProcessor<
                    FASTA_PE, BitStr, 
                    BIO_filtering_istream, BIO_filtering_ostream
                > INIT_TASK_PROCESSOR;
                task_processor.process();
            }
            else
            {
                TaskProcessor<
                    FASTA_PE, BitStr, 
                    BIO_filtering_istream, std::ofstream
                > INIT_TASK_PROCESSOR;
                task_processor.process();
            }
        }
        else
        {
            if (is_gz_output)
            {
                TaskProcessor<
                    FASTA_PE, BitStr, 
                    std::ifstream, BIO_filtering_ostream
                > INIT_TASK_PROCESSOR;
                task_processor.process();
            }
            else
            {
                TaskProcessor<
                    FASTA_PE, BitStr, 
                    std::ifstream, std::ofstream
                > INIT_TASK_PROCESSOR;
                task_processor.process();
            }
        }
    }
}

}