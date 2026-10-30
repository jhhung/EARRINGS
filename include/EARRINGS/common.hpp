#pragma once
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <optional>
#include <string>
#include <vector>
#include <string_view>
#include <utility>
#include <filesystem>
#include <EARRINGS/version.h>

// first whitespace-delimited token of a read header, used as the join key
// between <prefix>.tags.tsv and the trimmed reads / original input.
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
size_t fix_ext_len(0);
size_t init_kmer_size(10);
size_t kmer_step(5);
// Set (bare --sam) means write <output>.sam; the value is the mapped-read
// cap for that pass, with 0 meaning unlimited (whole file).
std::optional<std::size_t> sam_cap;
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

// -o prefix without the .fastq/.fasta suffix init_single appends to ofs_name[0];
// side outputs are named <prefix>.<kind> (e.g. .sam, .tags.tsv).
inline std::string output_prefix() {
    std::string base = ofs_name[0];
    for (std::string_view ext : {".fastq", ".fasta"})
        if (base.size() > ext.size() && std::string_view(base).substr(base.size() - ext.size()) == ext)
            base.resize(base.size() - ext.size());
    return base;
}

// Default -o for the single-end modes: the input's file name without its
// directory and read-file extensions (e.g. data/S15_1.fastq.gz -> S15_1).
inline std::string default_output_prefix(const std::string& input) {
    auto name = std::filesystem::path(input).filename().string();
    auto strip = [&](std::string_view ext) {
        if (name.size() > ext.size() && std::string_view(name).substr(name.size() - ext.size()) == ext) {
            name.resize(name.size() - ext.size());
            return true;
        }
        return false;
    };
    strip(".gz");
    for (std::string_view ext : {".fastq", ".fq", ".fasta", ".fa", ".ubam", ".bam"})
        if (strip(ext)) break;
    return name;
}

// Default -o for paired mode: the R1 name as above, minus a trailing mate
// marker, plus "_trimmed" so <prefix>_1.fastq cannot be the input itself
// (data/S15_R1.fastq.gz -> S15_trimmed).
inline std::string default_paired_output_prefix(const std::string& input1) {
    auto name = default_output_prefix(input1);
    for (std::string_view mate : {"_R1_001", "_R1", ".R1", "_1", ".1"})
        if (name.size() > mate.size() && std::string_view(name).substr(name.size() - mate.size()) == mate) {
            name.resize(name.size() - mate.size());
            break;
        }
    return name + "_trimmed";
}

// Private temporary directory for this run (under $TMPDIR, else /tmp),
// created on first use and removed at normal exit, so concurrent runs never
// share temporary files.
inline const std::filesystem::path& scratch_dir() {
    // Never destroyed: a static path would be destroyed before the atexit
    // handler (registered during its own initialization) runs.
    static const auto* dir = [] {
        auto tmpl = (std::filesystem::temp_directory_path() / "EARRINGS.XXXXXX").string();
        if (!mkdtemp(tmpl.data()))
            throw std::runtime_error("Cannot create a temporary directory from " + tmpl);
        std::atexit([] {
            std::error_code ec;
            std::filesystem::remove_all(scratch_dir(), ec);
        });
        return new std::filesystem::path(tmpl);
    }();
    return *dir;
}

// skewer's trimmed reads, as named from the -o it is given (ofs_name[0]).
inline std::string trimmed_reads_path(bool fastq) {
    return ofs_name[0] + (fastq ? "-trimmed.fastq" : "-trimmed.fasta");
}

// Exit before any work if an output is an input file (also via ./, symlinks or
// hard links): opening it for writing would truncate the input first.
inline void ensure_outputs_are_not_inputs(const std::vector<std::string>& outputs,
                                          const std::vector<std::string>& inputs) {
    for (const auto& out : outputs)
        for (const auto& in : inputs) {
            std::error_code ec;
            if (std::filesystem::equivalent(out, in, ec)) {
                std::cerr << "Error: output " << out << " is the input file " << in << "; choose another -o\n";
                std::exit(1);
            }
        }
}

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
