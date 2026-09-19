#pragma once
#include <cstdint>
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
// FMD index over F + revcomp(F) doubles the text length, so a whole-genome
// reference (e.g. GRCh38, ~3.1 Gbp -> ~6.2 Gbp text) can overflow uint32_t
// (max ~4.29 Gbp). Above GENOME_SCALE_BP, `build` picks the 64-bit
// instantiation instead; below it, 32-bit halves the index's memory
// footprint. The flag is persisted in the .table file (see EARRINGS.cpp /
// SE_auto_detect.hpp) so `single`/`smallRNA` know which to load() into.
using IndexSorter32 = biovoltron::PsaisSorter<std::uint32_t>;
using IndexSorter64 = biovoltron::PsaisSorter<std::uint64_t>;
using DenseIndex32   = biovoltron::BidirectionalIndex<true,  std::uint32_t, IndexSorter32>;
using SampledIndex32 = biovoltron::BidirectionalIndex<false, std::uint32_t, IndexSorter32>;
using DenseIndex64   = biovoltron::BidirectionalIndex<true,  std::uint64_t, IndexSorter64>;
using SampledIndex64 = biovoltron::BidirectionalIndex<false, std::uint64_t, IndexSorter64>;

// Shared "genome-scale" cutover: below it SA_INTV=256, above it SA_INTV=1024
// (see pick_sa_intv_bucket), and it also doubles as the 32-bit/64-bit index
// size_type boundary (see pick_use_64bit_index) -- one number to reason
// about instead of two independently-tuned ones.
constexpr std::size_t GENOME_SCALE_BP = 2'000'000'000;

// SA_INTV bucket, chosen at `build` time from the reference's total bp and
// persisted in the .table file (see EARRINGS.cpp / SE_auto_detect.hpp).
// Each value is an empirical per-scale winner, not interpolated; other
// candidates (32/64/128/512) never won at any tested scale, hence absent.
// Machine-dependent (NUMA, cache size, thread count) -- see sa_intv_bench/.
struct SaIntvBucket { bool dense; int sa_intv; };
inline SaIntvBucket pick_sa_intv_bucket(std::size_t total_bp) {
    if (total_bp <= 100'000)         return {true, 1};
    if (total_bp <= 80'000'000)      return {false, 16};
    if (total_bp <= GENOME_SCALE_BP) return {false, 256};
    return {false, 1024};
}

inline bool pick_use_64bit_index(std::size_t total_bp) {
    return total_bp > GENOME_SCALE_BP;
}

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
