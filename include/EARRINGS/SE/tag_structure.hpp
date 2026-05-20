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
    auto parse_tag_structure = [&](auto& tag_structure, const std::string& structure_str, size_t& tags_total_len) {
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
            tags_total_len += length;
            tag_structure.emplace_back(name, length);
        }
    };
    parse_tag_structure(tag_structure5, structure_str5, tags5_total_len);
    if (!tag_structure5.empty()) skipped_5prime_len += tags5_total_len;
    parse_tag_structure(tag_structure3, structure_str3, tags3_total_len);
}

void trim_tags3(auto& tags3, const std::string& adapter3) {
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
            records.push_back(std::move(rec));
        }
        const auto num_records = records.size();

        std::vector<std::string> template_tags(num_tags3);
        for (auto i = size_t{}, adapter3_i = size_t{}, acc_len = size_t{}; i < num_tags3; ++i) {
            const auto len = tag_structure3[i].second;
            template_tags[i].resize(len);
            for (auto j = size_t{}; j < len; ++j, ++acc_len) {
                template_tags[i][j] = acc_len >= estimated_tags3_len ? adapter3[adapter3_i++] : ' ';
            }
        }
        tags3.assign(num_records, template_tags);

        #pragma omp parallel for schedule(static)
        for (auto idx = size_t{}; idx < num_records; ++idx) {
            auto& rec = records[idx];
            auto& rec_tags = tags3[idx];

            const auto seq_len = rec.seq.size();
            if (seq_len < estimated_tags3_len) continue;

            rec_tags.resize(num_tags3);
            auto remaining = estimated_tags3_len;
            for (auto i = size_t{}, seq_i = seq_len - estimated_tags3_len; i < num_tags3 && remaining > 0; ++i) {
                const auto copy_len = std::min(tag_structure3[i].second, remaining);
                rec_tags[i].resize(copy_len);
                for (auto j = size_t{}; j < copy_len; ++j) {
                    rec_tags[i][j] = rec.seq[seq_i++];
                }
                remaining -= copy_len;
            }
            if (umi_idx != num_tags3) rec.name += ":" + rec_tags[umi_idx];

            const auto insert_len = seq_len - estimated_tags3_len;
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

void export_tags_tsv(const auto& tags5, const auto& tags3) {
    std::ofstream ofs("tags.tsv");
    if (!ofs) return;

    auto write_row = [&](const auto&& view) {
        auto it = view.begin(), end = view.end();
        if (it == end) return;
        ofs << *it++;
        for (; it != end; ++it) {
            ofs << '\t' << *it;
        }
        ofs << '\n';
    };

    write_row(ranges::views::concat(
        tag_structure5 | ranges::views::transform(&std::pair<std::string, size_t>::first),
        tag_structure3 | ranges::views::transform(&std::pair<std::string, size_t>::first)
    ));
    for (const auto& [rec_tags5, rec_tags3] : ranges::views::zip(tags5, tags3)) {
        write_row(ranges::views::concat(rec_tags5, rec_tags3));
    }
}

}