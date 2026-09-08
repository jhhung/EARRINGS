#include <string>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <exception>
#include <charconv>
#include <boost/range/join.hpp>
#include <boost/algorithm/string/join.hpp>
#include <range/v3/all.hpp>
#include <omp.h>
#include <EARRINGS/common.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/file_io/fastq.hpp>

namespace EARRINGS {

void parse_tag_structures(const std::string& structure_str5, const std::string& structure_str3) {
    std::unordered_set<std::string> seen_tags;
    auto parse_tag_structure = [&](auto& tag_structure, const std::string& structure_str) {
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
            if (!seen_tags.insert(name).second) {
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
    const std::string filename = std::string("trimmed_se.") + (is_fastq ? "fastq-trimmed.fastq" : "fasta-trimmed.fasta");
    std::ifstream ifs(filename);
    if (!ifs) return;
    const std::string tmpfile = filename + ".tmp";
    std::ofstream ofs(tmpfile);
    if (!ofs) return;

    const auto num_tags3 = tag_structure3.size();
    const size_t umi_idx = ranges::find_if(tag_structure3, [](const auto& p) { return p.first == "UMI"; }) - tag_structure3.begin();

    auto extract_tags3 = [&]<typename Record>() {
        std::vector<Record> records;
        Record rec;
        while (ifs >> rec) {
            // Capture the name as skewer emitted it (i.e. exactly what
            // trim_heads_to_stream wrote), before trim_tags3's own UMI append.
            names3.push_back(read_qname(rec.name));
            records.push_back(std::move(rec));
        }
        const auto num_records = records.size();

        // Missing tags3 (read too short after adapter trimming) are left as
        // empty strings, an unambiguous "no data" marker distinct from any
        // real tag value (declared tag lengths are always > 0).
        tags3.assign(num_records, std::vector<std::string>(num_tags3));

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

        #pragma omp parallel for schedule(static)
        for (auto idx = size_t{}; idx < num_records; ++idx) {
            auto& rec = records[idx];
            auto& rec_tags = tags3[idx];

            const auto seq_len = rec.seq.size();
            if (seq_len < keep_on_read) continue;

            const auto insert_len = seq_len - keep_on_read;
            const auto full_tag = rec.seq.substr(insert_len, keep_on_read) + absorbed;
            for (auto i = size_t{}, pos = size_t{}; i < num_tags3; ++i) {
                const auto len = tag_structure3[i].second;
                if (pos + len <= full_tag.size()) rec_tags[i] = full_tag.substr(pos, len);
                pos += len;
            }
            if (umi_idx != num_tags3 && !rec_tags[umi_idx].empty()) rec.name += ":" + rec_tags[umi_idx];

            rec.seq = rec.seq.substr(0, insert_len);
            if constexpr (std::is_same_v<Record, biovoltron::FastqRecord<>>) rec.qual = rec.qual.substr(0, insert_len);
        }

        for (auto idx = size_t{}; idx < num_records; ++idx) {
            ofs << records[idx] << '\n';
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
    std::ofstream ofs("tags.tsv");
    if (!ofs) return;

    auto write_row = [&](const std::string& name, auto&& tag_view) {
        ofs << name;
        for (auto&& v : tag_view) ofs << '\t' << v;
        ofs << '\n';
    };

    ofs << "name";
    for (const auto& p : tag_structure5) ofs << '\t' << p.first;
    for (const auto& p : tag_structure3) ofs << '\t' << p.first;
    ofs << '\n';

    const std::vector<std::string> empty5(tag_structure5.size());
    const std::vector<std::string> empty3(tag_structure3.size());

    if (!names5.empty()) {
        // Drive off the 5' side (one entry per input read). skewer only drops
        // reads, never reorders, so a forward cursor over names3 stays aligned.
        auto j = size_t{};
        for (auto i = size_t{}; i < names5.size(); ++i) {
            const bool has3 = j < names3.size() && names3[j] == names5[i];
            write_row(names5[i], ranges::views::concat(tags5[i], has3 ? tags3[j] : empty3));
            if (has3) ++j;
        }
    } else {
        // No head trimming happened; only 3' tags are available.
        for (auto j = size_t{}; j < names3.size(); ++j) {
            write_row(names3[j], ranges::views::concat(empty5, tags3[j]));
        }
    }
}

}