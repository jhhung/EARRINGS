#pragma once
#include <EARRINGS/common.hpp>
#include <EARRINGS/PE/task_processor.hpp>
#include <Biovoltron/format/fastq.hpp>
#include <Biovoltron/format/fasta_peat.hpp>
#include <Nucleona/parallel/asio_pool.hpp>
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
    constexpr size_t const_thread_num = 2;
    std::tuple<float, float, float> trimmer_param = std::make_tuple(match_rate, seq_cmp_rate, adapter_cmp_rate);
    std::vector<std::string> default_adapter{DEFAULT_ADAPTER1, DEFAULT_ADAPTER2};

    if (is_fastq)
    {
        TaskProcessor<FASTQ, BitStr> task_processor(thread_num
                                                , block_size
                                                , record_line
                                                , ifs_name
                                                , ofs_name
                                                , trimmer_param
                                                , loc_tail
                                                , default_adapter
                                                , DETECT_N_READS
                                                , is_sensitive);
        task_processor.process();
    }
    else
    {
        TaskProcessor<FASTA_PE, BitStr> task_processor(thread_num
                                                , block_size
                                                , record_line
                                                , ifs_name
                                                , ofs_name
                                                , trimmer_param
                                                , loc_tail
                                                , default_adapter
                                                , DETECT_N_READS
                                                , is_sensitive);
        task_processor.process();
    }
    
}

}