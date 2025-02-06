#pragma once

#include <EARRINGS/common.hpp>
#include <EARRINGS/PE/task_processor.hpp>
#include <EARRINGS/PE/format/fasta.hpp>
#include <EARRINGS/PE/format/fastq.hpp>
#include <Nucleona/parallel/asio_pool.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <tuple>
#include <iostream>

using namespace EARRINGS;
using namespace nucleona::parallel;

using BitStr = EARRINGS::vector<
        EARRINGS::char_type,
        simdpp::aligned_allocator<simdpp::uint8 < 16>, 16>
>;

namespace EARRINGS {
    void PE_trim() {
        using BIO_filtering_istream = boost::iostreams::filtering_istream;
        using BIO_filtering_ostream = boost::iostreams::filtering_ostream;

        constexpr
        size_t const_thread_num = 2;
        std::tuple<float, float, float> trimmer_param = std::make_tuple(match_rate, seq_cmp_rate, adapter_cmp_rate);
        std::vector <std::string> default_adapter{DEFAULT_ADAPTER1, DEFAULT_ADAPTER2};

        #define INIT_TASK_PROCESSOR task_processor(thread_num,\
                                                   block_size,\
                                                   record_line,\
                                                   ifs_name,\
                                                   ofs_name,\
                                                   trimmer_param,\
                                                   loc_tail,\
                                                   default_adapter,\
                                                   DETECT_N_READS,\
                                                   is_sensitive)

        if (is_fastq) {
            record_line = 4;
            if (is_gz_input) {
                if (is_gz_output) {
                    EARRINGS::TaskProcessor<EARRINGS::format::Fastq, BitStr, BIO_filtering_istream, BIO_filtering_ostream> INIT_TASK_PROCESSOR;
                    task_processor.process();
                } else {
                    EARRINGS::TaskProcessor<EARRINGS::format::Fastq, BitStr, BIO_filtering_istream, std::ofstream> INIT_TASK_PROCESSOR;
                    task_processor.process();
                }
            } else {
                if (is_gz_output) {
                    EARRINGS::TaskProcessor<EARRINGS::format::Fastq, BitStr, std::ifstream, BIO_filtering_ostream> INIT_TASK_PROCESSOR;
                    task_processor.process();
                } else {
                    EARRINGS::TaskProcessor<EARRINGS::format::Fastq, BitStr, std::ifstream, std::ofstream> INIT_TASK_PROCESSOR;
                    task_processor.process();
                }
            }
        } else {
            record_line = 2;
            if (is_gz_input) {
                if (is_gz_output) {
                    EARRINGS::TaskProcessor<EARRINGS::format::Fasta, BitStr, BIO_filtering_istream, BIO_filtering_ostream> INIT_TASK_PROCESSOR;
                    task_processor.process();
                } else {
                    EARRINGS::TaskProcessor<EARRINGS::format::Fasta, BitStr, BIO_filtering_istream, std::ofstream> INIT_TASK_PROCESSOR;
                    task_processor.process();
                }
            } else {
                if (is_gz_output) {
                    EARRINGS::TaskProcessor<EARRINGS::format::Fasta, BitStr, std::ifstream, BIO_filtering_ostream> INIT_TASK_PROCESSOR;
                    task_processor.process();
                } else {
                    EARRINGS::TaskProcessor<EARRINGS::format::Fasta, BitStr, std::ifstream, std::ofstream> INIT_TASK_PROCESSOR;
                    task_processor.process();
                }
            }
        }
    }
}