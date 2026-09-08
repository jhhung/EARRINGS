#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <string_view>
#include <utility>
#include <filesystem>
#include <EARRINGS/version.h>

// first whitespace-delimited token of a read header, used as the join key
// between tags.tsv and the trimmed reads / original input.
inline std::string read_qname(const std::string& name) {
    return name.substr(0, name.find(' '));
}

#define GET_STR(arg)			#arg
#define GET_VERSION(ver)		GET_STR(ver)
#define GET_EARRINGS_VERSION	GET_VERSION(EARRINGS_VERSION)

// for both SE and PE
bool is_fastq = true;
bool is_sensitive = false;
bool is_gz_input(false), is_gz_output(false);
bool is_bam(false);
size_t record_line = 4;
constexpr size_t DETECT_N_READS = 10000;

// for SE
constexpr int SA_INTV = 64;
std::string index_prefix;
size_t seed_len(18);
size_t min_multi(0);
float prune_factor(0.1);
bool no_mismatch(false);
size_t skipped_5prime_len(7);
size_t init_kmer_size(10);
size_t kmer_step(5);
std::vector<std::pair<std::string, size_t>> tag_structure5;
size_t tags5_total_len(0);
std::vector<std::pair<std::string, size_t>> tag_structure3;
size_t tags3_total_len(0);
size_t estimated_tags3_len(0);
// How the declared 3' tag region is split between two sources (set during
// adapter detection, consumed by trim_tags3):
//   tags3_kept_on_read  - trailing bases of skewer's output that are genuine
//                         tag bases and must be peeled from each read.
//   tags3_absorbed_prefix - the innermost, ~constant tag bases that got
//                         assembled into adapter3 and are therefore no longer
//                         on the read; the same value for every read, taken
//                         from adapter3's prefix (insert->adapter order).
// tags3_kept_on_read + tags3_absorbed_prefix.size() == tags3_total_len, except
// when extraction is disabled (unreliable estimate) where both are empty/0.
size_t tags3_kept_on_read(0);
std::string tags3_absorbed_prefix;

// for PE
size_t thread_num(1);
size_t block_size(8192);
size_t min_length(0);
std::vector<std::string> ifs_name(2);
std::vector<std::string> ofs_name(2);

// for smallRNA
size_t min_seed_len(21);
size_t max_seed_len(25);

bool loc_tail(true);  // adapter locates at tail/head
// size_t umi_loc(0);  // 1/2/3: umi seq locates at read1/read2/both
// size_t umi_len(0);
float match_rate = 0.7, seq_cmp_rate = 0.9, adapter_cmp_rate = 0.8;

// default adapters
std::string DEFAULT_ADAPTER1("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC");
std::string DEFAULT_ADAPTER2("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA");
