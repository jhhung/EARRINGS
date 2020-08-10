#pragma once
#include <cinttypes>
#include <fstream>
#include <vector>
#include <memory>
#include <Nucleona/language.hpp>
#include <Nucleona/range.hpp>
#include <mutex>
#include "suffix.hpp"
#include <range/v3/istream_range.hpp>
#include <iostream>
#include <boost/filesystem.hpp>
#include <Nucleona/container/vector.hpp>
#include <Nucleona/proftool/gprofiler.hpp>
#include <Nucleona/proftool/timer.hpp>
#include "group_file.hpp"
namespace nucleona::algo::string_sort {
template<class Str>
auto write_helper(std::ostream& out, const Str& str) {
    // std::cout << str.idx() << '\n';
    out.write((char*)&str.idx(), sizeof(str.idx()));
};
template<class Str>
auto read_helper(std::istream& in, Str& str) {
    in.read((char*)&str.idx(), sizeof(str.idx()));
}
struct HashSlot {
    template<class T>
    inline void push_strs(const std::vector<Suffix<T>>& strs) {
        if(!group_file) {
            group_file.reset(new std::ofstream(group_filename, std::ios::binary));
        }
        for(auto&& s : strs) {
            // std::cout << group_filename << '\t';
            nucleona::algo::string_sort::write_helper(*group_file, s);
        }
        num += strs.size();
    }
    inline void close() {
        if(group_file) {
            group_file->close();
        }
    }
    inline void flush() {
        if(group_file) {
            group_file->flush();
        }
    }

    std::size_t                     num             {0};
    std::string                     group_filename  {""};
    std::unique_ptr<std::ofstream>  group_file      {nullptr};
    std::unique_ptr<std::mutex>     mux             {new std::mutex()};

    std::unique_ptr<
        std::vector<HashSlot>
    >                               sub_table       {new std::vector<HashSlot>()};
};
template<class Hasher, class IdxInt = std::size_t>
struct HashMSDRadixFileSplit {
    HashMSDRadixFileSplit(
        const Hasher& hasher, 
        int recursive_max_level = 2, 
        std::size_t group_size_max = 1024 * 1024 * 1024 / sizeof(IdxInt) 
    )
    : hasher_               (hasher)
    , recursive_max_level_  (recursive_max_level)
    , group_size_max_       (group_size_max)
    {}
    template<class StringRng>
    auto core_sync(
        StringRng&& seg, 
        std::vector<HashSlot>& table, int level = 0
    ) const {
        using String = nucleona::range::ValueT<StringRng>;
        auto size = ranges::distance(seg);
        auto timer = nucleona::proftool::make_timer([&size](auto du){

            std::cout << "classify size: " << size << "    " << std::chrono::duration_cast<
                std::chrono::milliseconds
            >(du).count() << "ms\n";
        });
        std::vector<std::vector<String>> buffer(hasher_.table_size());
        for(auto&& ch : seg | ranges::view::chunk(256 * 1024)) {
            for(auto&& str : ch) {
                auto table_idx = hasher_(str, level);
                buffer[table_idx].push_back(str);
            }
            std::size_t buffer_elements = 0;
            for(auto&& str_ids : buffer) {
                buffer_elements += str_ids.size();
            }
            while(buffer_elements > 0) {
                for(auto&& [table_idx, str_ids] : buffer | nucleona::range::indexed()) {
                    if(str_ids.size() == 0) continue;
                    auto& slot = table[table_idx];
                    if(!slot.mux->try_lock()) continue;
                    slot.push_strs(str_ids);
                    slot.mux->unlock();
                    buffer_elements -= str_ids.size();
                    str_ids.clear();
                }
            }
        }
    }
    template<class Executor, class StringRng>
    auto core(Executor& executor, StringRng&& str_rng, std::vector<HashSlot>& table, int level = 0) const {
        using namespace nucleona::container;
        std::vector<typename Executor::template Future<void>> futures;
        futures.reserve(executor.get_parallization() + 1);
        for(auto&& seg : str_rng 
            | nucleona::range::segment(executor.get_parallization() + 1)
        ) {
            futures.push_back(executor.submit([seg, this, &table, &level](){
                core_sync(seg, table, level);
            }));
        }
        for(auto&& fut : futures) {
            fut.sync();
        }
    }
    template<class Executor, class StringRng>
    auto run_with_create_groups(
        Executor& executor, StringRng&& str_rng, 
        int level = 0, const std::string& group_name_prefix = "group", bool close_file = true
    ) const {
        std::vector<HashSlot> table(hasher_.table_size());
        for(auto&& [idx, slot] : table | nucleona::range::indexed(0)) {
            slot.group_filename = group_name_prefix + "-" + std::to_string(idx);
        }
        core(executor, str_rng, table, level);
        if(close_file) {
            for(auto&& slot : table ) {
                slot.close();
            }
        } else {
            for(auto&& slot : table ) {
                slot.flush();
            }
        }
        return table;
    }
    template<class Executor, class StrRng>
    auto run_with_file_split(Executor& executor, HashSlot& slot, StrRng&& strs, int level = 0) const {
        // read a chunk of index from file and split the chunk into next level files
        // this doing chunk instead of file range is becuase the "split step is parallel"
        const std::size_t chunk_size = 1024 * 1024 * 128; // 128M
        std::vector<HashSlot> res;
        std::ifstream fin(slot.group_filename, std::ios::binary);
        for(auto&& ch : nucleona::range::irange_0(slot.num)
            | ranges::view::chunk(chunk_size) 
        ) {
            std::vector<IdxInt> idx_vec;
            idx_vec.reserve(chunk_size);
            for(auto&& i : ch) {
                IdxInt idx;
                fin.read((char*)&idx, sizeof(IdxInt));
                idx_vec.push_back(idx);
            }
            auto inner_strs = ranges::view::transform(idx_vec, [&strs](auto&& idx){
                return strs[idx];
            });
            if(res.empty()) {
                res = run_with_create_groups(executor, inner_strs, 
                    level, slot.group_filename, false);
            } else {
                core(executor, strs, res, level );
            }
        }
        for(auto&& slot : res) {
            slot.group_file->close();
        }
        boost::filesystem::remove(slot.group_filename);
        return res;
    }
    template<class Executor, class StrRng>
    auto run_with_strs_split(Executor& executor, StrRng&& strs, int level = 0) const {
        auto groups = run_with_create_groups(executor, strs, level, "group", true);
        return groups;
    }
    template<class Executor, class StrRng>
    std::vector<HashSlot> recursive_file_split(
        Executor& executor, HashSlot& slot, StrRng&& strs, int level
    ) const {
        if(level > recursive_max_level_) return {}; // TODO: Threshold parameterize
        if(slot.num < group_size_max_ || slot.num < 1) return {};
        auto tab = run_with_file_split(executor, slot, strs, level);
        for(auto&& sub_slot : tab) {
            sub_slot.sub_table.reset(new std::vector<HashSlot>(
                recursive_file_split(
                    executor, sub_slot, strs, level + 1
                )
            ));
        }
        return tab;
    }
    template<class Executor, class StrRng>
    auto group_strings(Executor& executor, StrRng&& strs) const {
        auto groups = run_with_create_groups(executor, strs, 0, "group", true);
        for(auto&& slot : groups) {
            slot.sub_table.reset(new std::vector<HashSlot>(
                recursive_file_split(executor, slot, strs, 1)
            ));
        }
        return groups;
    }
    template<class Executor, class String>
    auto group_suffixes(Executor& executor, String&& str) const {
        auto strs = nucleona::range::irange_0(
                IdxInt(str.size())
            )
            | nucleona::range::transformed([&str](auto&& idx){
                return make_suffix(str, idx);
            })
        ;
        return group_strings(executor, strs);
    }
    void file_flatten_impl(
        const std::vector<HashSlot>&    tab, 
        std::vector<GroupFile>&         file_list,
        std::size_t                     level
    ) const {
        for(auto&& slot : tab) {
            if(slot.sub_table->size() > 0) {
                file_flatten_impl(*slot.sub_table, file_list, level + 1);
            } else {
                if(slot.num > 0) {
                    file_list.emplace_back(GroupFile{slot.group_filename, slot.num, level});
                } else {
                    boost::filesystem::remove(slot.group_filename);
                }
            }
        }
    }
    auto file_flatten(const std::vector<HashSlot>& tab) const {
        std::vector<GroupFile> res;
        file_flatten_impl(tab, res, 0);
        return res;
    }
    template<class Executor, class Arg>
    auto operator()(Executor& executor, Arg&& arg) const {
        if constexpr(std::is_same_v<nucleona::range::ValueT<std::decay_t<Arg>>, char>) {
            return file_flatten(group_suffixes(executor, FWD(arg)));
        } else {
            return file_flatten(group_strings(executor, FWD(arg)));
        }
    }
    
private:
    Hasher      hasher_;
    int         recursive_max_level_;
    std::size_t group_size_max_;
};
}
