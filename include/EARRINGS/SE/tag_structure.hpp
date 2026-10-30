#include <string>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <charconv>
#include <boost/range/join.hpp>
#include <boost/algorithm/string/join.hpp>
#include <range/v3/all.hpp>
#include <omp.h>
#include <EARRINGS/common.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/file_io/fastq.hpp>

namespace EARRINGS {

// Tag names must be unique, except that "UMI" may appear once on each end: a
// UMI split over both ends (e.g. 5N + 5N) is reported as one value.
void parse_tag_structures(const std::string& structure_str5, const std::string& structure_str3) {
    std::unordered_set<std::string> seen_tags;
    auto parse_tag_structure = [&](auto& tag_structure, const std::string& structure_str) {
        std::unordered_set<std::string> seen_here;
        std::istringstream ss(structure_str);
        std::string tag;
        while (std::getline(ss, tag, ',')) {
            const auto colon_pos = tag.find(':');
            if (colon_pos == std::string::npos) {
                throw std::invalid_argument("Missing ':' separator when parsing the tag structure string.");
            }
            auto length = size_t{};
            const auto [ptr, ec] = std::from_chars(tag.data() + colon_pos + 1, tag.data() + tag.size(), length);
            if (ec != std::errc() || length == 0) {
                throw std::invalid_argument("Invalid length value in tag structure.");
            }
            const auto name = tag.substr(0, colon_pos);
            const bool repeated_here = !seen_here.insert(name).second;
            if (repeated_here || (!seen_tags.insert(name).second && name != "UMI")) {
                throw std::invalid_argument("Duplicate tag name in tag structures: " + name);
            }
            tag_structure.emplace_back(name, length);
        }
    };
    parse_tag_structure(tag_structure5, structure_str5);
    for (const auto& [name, length] : tag_structure5) tags5_total_len += length;
    if (!tag_structure5.empty()) skipped_5prime_len += tags5_total_len;
    parse_tag_structure(tag_structure3, structure_str3);
    for (const auto& [name, length] : tag_structure3) tags3_total_len += length;
}

void trim_tags3(auto& tags3, auto& names3) {
    const std::string filename = trimmed_reads_path(is_fastq);
    std::ifstream ifs(filename);
    if (!ifs) throw std::runtime_error("Cannot open " + filename + " to extract 3' tags");
    const std::string tmpfile = filename + ".tmp";
    std::ofstream ofs(tmpfile);
    if (!ofs) throw std::runtime_error("Cannot open " + tmpfile + " for writing");

    const auto num_tags3 = tag_structure3.size();
    const auto umi_idx = umi_index(tag_structure3);
    // A UMI declared on the 5' end but not put into the read IDs (no 5' region
    // was detected) is padded with N, so every read carries a UMI of one length.
    const auto umi5_padding = umi_index(tag_structure5) != tag_structure5.size() && !umi5_in_read_id
        ? std::string(tag_structure5[umi_index(tag_structure5)].second, 'N') : std::string{};

    // Reads are processed in batches so only one batch of records is in memory.
    constexpr auto batch_size = size_t{1} << 20;

    auto extract_tags3 = [&]<typename Record>() {
        // The declared tag region is split by adapter detection into:
        //   - tags3_kept_on_read  trailing bases still on each skewer-trimmed
        //                         read (peeled per-read here), and
        //   - tags3_absorbed_prefix  the innermost ~constant bases, folded into
        //                         adapter3 and hence already gone from the read;
        //                         the same value for every read.
        // Their concatenation, in insert->adapter order, is the full tag region.
        // When extraction was disabled (unreliable estimate) tags3_kept_on_read
        // is 0 and tags3_absorbed_prefix empty, so every tag comes out empty and
        // reads are left untouched.
        const auto& absorbed = tags3_absorbed_prefix;
        const auto keep_on_read = tags3_kept_on_read;

        std::vector<Record> records;
        records.reserve(batch_size);
        Record next;
        for (bool more = true; more;) {
            records.clear();
            while (records.size() < batch_size && (more = static_cast<bool>(ifs >> next))) {
                // Capture the name as skewer emitted it (i.e. exactly what
                // trim_heads_to_stream wrote), before trim_tags3's own UMI append.
                names3.push_back(read_qname(next.name));
                records.push_back(std::move(next));
            }
            const auto num_records = records.size();
            const auto first = tags3.size();

            // Missing tags3 (read too short after adapter trimming) are left as
            // empty strings, an unambiguous "no data" marker distinct from any
            // real tag value (declared tag lengths are always > 0).
            tags3.resize(first + num_records, std::vector<std::string>(num_tags3));

            #pragma omp parallel for schedule(static)
            for (auto idx = size_t{}; idx < num_records; ++idx) {
                auto& rec = records[idx];
                auto& rec_tags = tags3[first + idx];

                const auto seq_len = rec.seq.size();
                if (seq_len >= keep_on_read) {
                    const auto insert_len = seq_len - keep_on_read;
                    const auto full_tag = rec.seq.substr(insert_len, keep_on_read) + absorbed;
                    for (auto i = size_t{}, pos = size_t{}; i < num_tags3; ++i) {
                        const auto len = tag_structure3[i].second;
                        if (pos + len <= full_tag.size()) rec_tags[i] = full_tag.substr(pos, len);
                        pos += len;
                    }
                    rec.seq = rec.seq.substr(0, insert_len);
                    if constexpr (std::is_same_v<Record, biovoltron::FastqRecord<>>) rec.qual = rec.qual.substr(0, insert_len);
                }
                if (umi_idx != num_tags3) {
                    // A missing 3' UMI (read too short, or extraction disabled) is all N.
                    auto umi = rec_tags[umi_idx];
                    umi.resize(tag_structure3[umi_idx].second, 'N');
                    add_umi_to_read_id(rec.name, umi5_padding + umi, !umi5_in_read_id);
                }
            }

            for (auto idx = size_t{}; idx < num_records; ++idx) {
                ofs << records[idx] << '\n';
            }
        }
    };

    if (is_fastq) {
        extract_tags3.template operator()<biovoltron::FastqRecord<>>();
    } else {
        extract_tags3.template operator()<biovoltron::FastaRecord<>>();
    }

    std::filesystem::rename(tmpfile, filename);
}

void export_tags_tsv(const auto& names5, const auto& tags5, const auto& names3, const auto& tags3) {
    const auto tsv_path = output_prefix() + ".tags.tsv";
    std::ofstream ofs(tsv_path);
    if (!ofs) throw std::runtime_error("Cannot open " + tsv_path + " for writing");

    // A UMI on both ends gets one column per end (UMI_5p, UMI_3p) plus a UMI
    // column with the combined value as it appears in the read IDs.
    const auto umi5_idx = umi_index(tag_structure5);
    const auto umi3_idx = umi_index(tag_structure3);
    const bool split_umi = umi5_idx != tag_structure5.size() && umi3_idx != tag_structure3.size();

    auto write_row = [&](const std::string& name, const auto& row5, const auto& row3) {
        ofs << name;
        for (const auto& v : row5) ofs << '\t' << v;
        for (const auto& v : row3) ofs << '\t' << v;
        if (split_umi) {
            auto umi5 = row5[umi5_idx], umi3 = row3[umi3_idx];
            umi5.resize(tag_structure5[umi5_idx].second, 'N');
            umi3.resize(tag_structure3[umi3_idx].second, 'N');
            ofs << '\t' << umi5 << umi3;
        }
        ofs << '\n';
    };

    ofs << "name";
    for (const auto& p : tag_structure5) ofs << '\t' << p.first << (split_umi && p.first == "UMI" ? "_5p" : "");
    for (const auto& p : tag_structure3) ofs << '\t' << p.first << (split_umi && p.first == "UMI" ? "_3p" : "");
    if (split_umi) ofs << "\tUMI";
    ofs << '\n';

    const std::vector<std::string> empty5(tag_structure5.size());
    const std::vector<std::string> empty3(tag_structure3.size());

    if (!names5.empty()) {
        // Drive off the 5' side (one entry per input read). skewer only drops
        // reads, never reorders, so a forward cursor over names3 stays aligned.
        // names3 come from skewer's output, whose IDs carry the 5' UMI.
        auto j = size_t{};
        for (auto i = size_t{}; i < names5.size(); ++i) {
            const auto written = umi5_in_read_id ? names5[i] + ":" + tags5[i][umi5_idx] : names5[i];
            const bool has3 = j < names3.size() && names3[j] == written;
            write_row(names5[i], tags5[i].empty() ? empty5 : tags5[i], has3 ? tags3[j] : empty3);
            if (has3) ++j;
        }
    } else {
        // No head trimming happened; only 3' tags are available.
        for (auto j = size_t{}; j < names3.size(); ++j) {
            write_row(names3[j], empty5, tags3[j]);
        }
    }
}

}