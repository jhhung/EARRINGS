#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <string_view>
#include <experimental/filesystem>

// for both SE and PE
bool is_fastq = true;
bool is_sensitive = false;
size_t record_line = 4;
constexpr size_t DETECT_N_READS = 10000;
std::string bam_fname;

// for SE
std::string index_prefix; 
size_t seed_len(18);
size_t min_multi(0);
float prune_factor(0.1);
bool enable_mismatch(true);
bool estimate_umi_len(false);

// for PE
size_t thread_num(1);
size_t block_size(8192);
size_t min_length(0);
std::vector<std::string> ifs_name(2);
std::vector<std::string> ofs_name(2);

bool loc_tail(true);  // adapter locates at tail/head
// size_t umi_loc(0);  // 1/2/3: umi seq locates at read1/read2/both
// size_t umi_len(0);
float match_rate = 0.7, seq_cmp_rate = 0.9, adapter_cmp_rate = 0.8;

// default adapters
std::string DEFAULT_ADAPTER1("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC");
std::string DEFAULT_ADAPTER2("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA");
