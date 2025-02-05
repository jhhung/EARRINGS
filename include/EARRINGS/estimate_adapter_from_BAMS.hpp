#pragma once

#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include <EARRINGS/common.hpp>
#include <EARRINGS/fasta.hpp>
#include <htslib/sam.h>
#include <biovoltron/utility/istring.hpp>

using namespace EARRINGS;

namespace EARRINGS {
/* 1. Read name, sequence and quality score from bam file.
   2. Store the reads.
*/
    class Process_uBAMs {
    private:
        inline static void write_to_fa(bam1_t* bam, std::ofstream &ofs, bool is_rev) {
            int seq_len = bam->core.l_qseq;
            std::string sequence;
            sequence.resize(seq_len);
            uint8_t *seq_ptr = bam_get_seq(bam);
            for (int i = 0; i < seq_len; ++i) {
                sequence[i] = seq_nt16_str[bam_seqi(seq_ptr, i)];
            }
            if (is_rev) {
                sequence = biovoltron::Codec::rev_comp(sequence);
            }
            ofs << ">" << bam_get_qname(bam) << "\n" << sequence << "\n";
        }

    public:
        // singled-end
        template<typename Str1, typename Str2>
        inline static size_t extract_reads_from_uBAMs(Str1 &&bam_fname, Str2 &&output_fname) {
            samFile *bam_file = sam_open(bam_fname.c_str(), "rb");
            if (!bam_file) {
                throw std::runtime_error("Can't open BAM files\n");
            }

            bam_hdr_t *bam_header = sam_hdr_read(bam_file);
            if (!bam_header) {
                sam_close(bam_file);
                throw std::runtime_error("Can't read BAM header\n");
            }

            std::ofstream ofs(output_fname);
            if (!ofs) {
                bam_hdr_destroy(bam_header);
                sam_close(bam_file);
                throw std::runtime_error("Can't open temp file which stores reads extract from bam/ubam\n");
            }

            bam1_t *bam = bam_init1();
            size_t counter = 0;
            while (sam_read1(bam_file, bam_header, bam) >= 0) {
                bool is_rev = (bam->core.flag & BAM_FREVERSE) != 0;    // == 16
                write_to_fa(bam, ofs, is_rev);
                ++counter;
            }

            bam_destroy1(bam);
            bam_hdr_destroy(bam_header);
            sam_close(bam_file);

            return counter;
        }

        // paired-end
        template<typename Str1, typename Str2, typename Str3>
        inline static size_t extract_reads_from_uBAMs(Str1&& bam_fname, Str2&& output_fname1, Str3&& output_fname2)
        {
            samFile *bam_file = sam_open(bam_fname.c_str(), "rb");
            if (!bam_file) {
                throw std::runtime_error("Can't open BAM files\n");
            }

            bam_hdr_t *bam_header = sam_hdr_read(bam_file);
            if (!bam_header) {
                sam_close(bam_file);
                throw std::runtime_error("Can't read BAM header\n");
            }

            std::ofstream ofs1(output_fname1);
            std::ofstream ofs2(output_fname2);
            if (!ofs1 || !ofs2) {
                bam_hdr_destroy(bam_header);
                sam_close(bam_file);
                throw std::runtime_error("Can't open input fastq/fasta files\n");
            }

            bam1_t *bam1 = bam_init1();
            bam1_t *bam2 = bam_init1();
            size_t counter = 0;
            while (sam_read1(bam_file, bam_header, bam1) >= 0) {
                if (sam_read1(bam_file, bam_header, bam2) < 0) {
                    break;
                }
                bool is_rev1 = (bam1->core.flag & BAM_FREVERSE) != 0;
                bool is_rev2 = (bam2->core.flag & BAM_FREVERSE) != 0;
                write_to_fa(bam1, ofs1, is_rev1);
                write_to_fa(bam2, ofs2, is_rev2);
                ++counter;
            }

            bam_destroy1(bam1);
            bam_destroy1(bam2);
            bam_hdr_destroy(bam_header);
            sam_close(bam_file);

            return counter;
        }
    };
}
