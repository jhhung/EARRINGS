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
#include <cstdio>
#include <cstring>
#include <memory>
#include <string_view>
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

// General version: parses with the Biovoltron stream operators, so it accepts any FASTA/FASTQ
// layout. trim_tags3 below only falls back to it for files that are not skewer's own output.
void trim_tags3_generic(auto& tags3, auto& names3) {
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

// Fast version for the file skewer just wrote, whose layout is fixed: 4-line FASTQ records
// ("@id\nseq\n+\nqual\n") or 2-line FASTA records (">id\nseq\n"), the ID line ending in '\n'
// (no comments of its own are lost: the ID line is kept whole). The file is read in large
// chunks, split into records, transformed on several threads and written back in order.
// Gives exactly what trim_tags3_generic gives; returns false (nothing written, tags3 and
// names3 as before) if any record does not have that layout, so the caller can fall back.
bool trim_tags3_fast(auto& tags3, auto& names3, const std::string& filename, const std::string& tmpfile,
                     const size_t chunk_bytes = size_t{32} << 20) {
    const bool fastq = is_fastq;
    const size_t lines_per_record = fastq ? 4 : 2;
    const auto num_tags3 = tag_structure3.size();
    const auto umi_idx = umi_index(tag_structure3);
    const bool has_umi = umi_idx != num_tags3;
    const auto umi5_padding = umi_index(tag_structure5) != tag_structure5.size() && !umi5_in_read_id
        ? std::string(tag_structure5[umi_index(tag_structure5)].second, 'N') : std::string{};
    const bool umi_new_field = !umi5_in_read_id;
    const auto& absorbed = tags3_absorbed_prefix;
    const auto keep_on_read = tags3_kept_on_read;
    const size_t umi_len = has_umi ? tag_structure3[umi_idx].second : 0;

    struct FileCloser { void operator()(FILE* f) const { if (f) fclose(f); } };
    std::unique_ptr<FILE, FileCloser> in(fopen(filename.c_str(), "rb"));
    if (!in) throw std::runtime_error("Cannot open " + filename + " to extract 3' tags");
    std::unique_ptr<FILE, FileCloser> out(fopen(tmpfile.c_str(), "wb"));
    if (!out) throw std::runtime_error("Cannot open " + tmpfile + " for writing");

    // One record: where its lines are in the chunk (offsets, so they survive buffer growth).
    struct Rec { size_t name, seq, qual; uint32_t name_len, seq_len; };

    const size_t num_threads = std::max<size_t>(1, thread_num);
    std::vector<char> buf;
    std::vector<Rec> recs;
    std::vector<std::string> outs(num_threads);
    const auto first_records = tags3.size();
    size_t carry = 0;       // bytes of an unfinished record kept from the previous chunk
    bool eof = false;
    bool ok = true;

    while (ok && !eof) {
        buf.resize(carry + chunk_bytes);
        const auto got = fread(buf.data() + carry, 1, chunk_bytes, in.get());
        if (got < chunk_bytes) {
            if (ferror(in.get())) throw std::runtime_error("Cannot read " + filename);
            eof = true;
        }
        const auto total = carry + got;

        // Split the complete records of this chunk.
        recs.clear();
        size_t pos = 0, used = 0;
        for (;;) {
            size_t line_start[4], line_len[4];
            auto p = pos;
            bool complete = true;
            for (size_t l = 0; l < lines_per_record; ++l) {
                const auto* nl = p < total ? static_cast<const char*>(memchr(buf.data() + p, '\n', total - p)) : nullptr;
                if (!nl) { complete = false; break; }
                line_start[l] = p;
                line_len[l] = static_cast<size_t>(nl - (buf.data() + p));
                p = static_cast<size_t>(nl - buf.data()) + 1;
            }
            if (!complete) break;
            const char* base = buf.data();
            if (fastq) {
                if (line_len[0] == 0 || base[line_start[0]] != '@' || line_len[2] != 1 || base[line_start[2]] != '+' ||
                    line_len[1] != line_len[3]) { ok = false; break; }
                recs.push_back({line_start[0] + 1, line_start[1], line_start[3], static_cast<uint32_t>(line_len[0] - 1), static_cast<uint32_t>(line_len[1])});
            } else {
                if (line_len[0] == 0 || base[line_start[0]] != '>' || (line_len[1] > 0 && base[line_start[1]] == '>')) { ok = false; break; }
                recs.push_back({line_start[0] + 1, line_start[1], 0, static_cast<uint32_t>(line_len[0] - 1), static_cast<uint32_t>(line_len[1])});
            }
            pos = used = p;
        }
        if (!ok) break;
        if (eof && used != total) { ok = false; break; }   // trailing partial record

        const auto n = recs.size();
        const auto first = tags3.size();
        names3.resize(first + n);
        tags3.resize(first + n);

        for (auto& o : outs) o.clear();   // a team smaller than num_threads leaves some unused
        #pragma omp parallel num_threads(num_threads)
        {
            const auto tid = static_cast<size_t>(omp_get_thread_num());
            const auto nt = static_cast<size_t>(omp_get_num_threads());
            const auto begin = n * tid / nt, end = n * (tid + 1) / nt;
            auto& o = outs[tid];
            o.clear();
            o.reserve((end - begin) * 120);
            std::string full_tag, umi;
            const char* base = buf.data();
            for (auto idx = begin; idx < end; ++idx) {
                const auto& r = recs[idx];
                const std::string_view name(base + r.name, r.name_len);
                const std::string_view seq(base + r.seq, r.seq_len);
                const auto sp = name.find(' ');
                names3[first + idx].assign(name.substr(0, sp == std::string_view::npos ? name.size() : sp));

                auto& rec_tags = tags3[first + idx];
                rec_tags.assign(num_tags3, std::string{});
                size_t out_len = seq.size();
                if (seq.size() >= keep_on_read) {
                    out_len = seq.size() - keep_on_read;
                    full_tag.assign(seq.substr(out_len, keep_on_read));
                    full_tag += absorbed;
                    for (auto i = size_t{}, tpos = size_t{}; i < num_tags3; ++i) {
                        const auto len = tag_structure3[i].second;
                        if (tpos + len <= full_tag.size()) rec_tags[i].assign(full_tag, tpos, len);
                        tpos += len;
                    }
                }

                o += fastq ? '@' : '>';
                if (has_umi) {
                    // A missing 3' UMI (read too short, or extraction disabled) is all N.
                    umi = rec_tags[umi_idx];
                    umi.resize(umi_len, 'N');
                    const auto id_end = sp == std::string_view::npos ? name.size() : sp;
                    o.append(name.substr(0, id_end));
                    if (umi_new_field) o += ':';
                    o += umi5_padding;
                    o += umi;
                    o.append(name.substr(id_end));
                } else {
                    o.append(name);
                }
                o += '\n';
                o.append(seq.substr(0, out_len));
                if (fastq) {
                    o += "\n+\n";
                    o.append(base + r.qual, out_len);
                }
                o += '\n';
            }
        }

        for (const auto& o : outs) {
            if (!o.empty() && fwrite(o.data(), 1, o.size(), out.get()) != o.size()) {
                throw std::runtime_error("Cannot write " + tmpfile);
            }
        }

        // Keep the unfinished tail for the next chunk.
        const auto rest = total - used;
        if (rest > 0) memmove(buf.data(), buf.data() + used, rest);
        carry = rest;
    }

    if (ok && fflush(out.get()) != 0) throw std::runtime_error("Cannot write " + tmpfile);
    if (!ok) {
        out.reset();
        std::filesystem::remove(tmpfile);
        tags3.resize(first_records);
        names3.resize(first_records);
    }
    return ok;
}

void trim_tags3(auto& tags3, auto& names3) {
    const std::string filename = trimmed_reads_path(is_fastq);
    const std::string tmpfile = filename + ".tmp";
    if (trim_tags3_fast(tags3, names3, filename, tmpfile)) {
        std::filesystem::rename(tmpfile, filename);
    } else {
        trim_tags3_generic(tags3, names3);
    }
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