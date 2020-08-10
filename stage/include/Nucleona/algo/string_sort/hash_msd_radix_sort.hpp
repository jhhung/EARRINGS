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
#include <Nucleona/pointer.hpp>
#include <Nucleona/container/array.hpp>

namespace nucleona::algo::string_sort {
template<class Hasher, class IdxInt = std::size_t>
struct HashMSDRadixSort {
    HashMSDRadixSort(
        const Hasher&   hasher,
        std::size_t     max_level
    )
    : hasher_(hasher)
    , max_level_(max_level)
    {}
    template<class Seg, class StrRng, class Buks, class BukMuxs>
    void core_task(
        Seg&&           seg,
        const StrRng&   str_rng,
        Buks*           buckets,
        BukMuxs*        bucket_muxs,
        int             level
    ) const {
        using namespace nucleona::container;
        // auto timer = nucleona::proftool::make_timer([&](auto du){
        //     std::cout << "task size: " << ranges::distance(seg) << '\t' 
        //         << std::chrono::duration_cast<std::chrono::milliseconds>(du).count() << "ms" << std::endl;
        // });
        static thread_local Array<std::vector<IdxInt>, Hasher::table_size_()> thread_buckets;

        for(auto&& i : seg) {
            auto b_idx = hasher_(str_rng[i], level);
            thread_buckets[b_idx].push_back(i);
        }
        auto thread_buckets_size = ranges::distance(seg);
        while(thread_buckets_size > 0) {
            for(std::size_t b_i = 0; b_i < hasher_.table_size(); b_i ++) {
                auto& b = thread_buckets[b_i];
                if(b.size() == 0) continue;
                auto& bucket = buckets[b_i];
                auto& bucket_mux = bucket_muxs[b_i];
                if(!bucket_mux.try_lock()) continue;
                bucket.insert(bucket.end(), b.begin(), b.end());
                bucket_mux.unlock();
                thread_buckets_size -= b.size();
                b.clear();
            }
        }
        for(std::size_t i = 0; i < Hasher::table_size_(); i ++ ) {
            auto& tb = thread_buckets[i];
            tb.clear();
        }
    }
    template<class StrRng, class Itr>
    void two_str_opt(
        StrRng&& str_rng,
        Itr&& beg, Itr&& end, 
        int level
    ) const {
        auto& i_a = *beg;
        auto& i_b = *(beg + 1);
        for(int l = level; l < max_level_; l ++) {
            auto v_a = hasher_(str_rng[i_a], l);
            auto v_b = hasher_(str_rng[i_b], l);
            if(v_a > v_b) {
                std::swap(i_a, i_b);
                return ;
            } else if( v_a == v_b ) {
                continue;
            } else {
                return ;
            }
        }
    }
    template<class Executor, class Itr, class StrRng>
    bool core(
        Executor& executor, 
        Itr ids_beg, Itr ids_end, 
        const StrRng& str_rng, 
        int level
    ) const {
        using namespace nucleona::container;
        if(level >= max_level_) return true;
        auto ids_num = std::distance(ids_beg, ids_end);
        if(ids_num <= 1) return true;
        if(ids_num == 2) { 
            two_str_opt(str_rng, ids_beg, ids_end, level); 
            return true;
        }
        ranges::iterator_range<std::decay_t<Itr>> ids_rng(ids_beg, ids_end);


        std::vector<typename Executor::template Future<bool>> local_futs;
        local_futs.reserve(4);

        static thread_local auto p_futures = nucleona::mkuniq(
            new std::vector<typename Executor::template Future<bool>>()
        );
        auto& futures = *p_futures;
        if(futures.capacity() == 0) {
            futures.reserve(executor.get_parallization());
        }
        static thread_local auto buckets = nucleona::mkuniq_arr(
            new std::vector<IdxInt>[hasher_.table_size()]
        );
        static thread_local auto bucket_muxs = nucleona::mkuniq_arr(
            new std::mutex[hasher_.table_size()]
        );
        static thread_local bool bucket_using = false;
        if(bucket_using) return false;
        bucket_using = true;
        for(std::size_t i = 0; i < hasher_.table_size(); i ++ ) {
            buckets[i].clear();
        }

        static constexpr std::size_t task_size = Hasher::table_size_() * 1024 * 128;

        bool no_task_split = ids_num < task_size;
        if(no_task_split) {
            core_task(ids_rng, str_rng, buckets.get(), bucket_muxs.get(), level);
        } else {
            for(auto&& seg : ids_rng 
                | nucleona::range::segment(executor.get_parallization() + 1)
            ) {
                futures.push_back(executor.submit([
                    seg, this, &str_rng, level, 
                    __buckets = buckets.get(), __bucket_muxs = bucket_muxs.get()
                ](){
                    core_task(seg, str_rng, __buckets, __bucket_muxs, level);
                    return true;
                }));
            }
            for(auto&& fut : futures) {
                fut.sync();
            }
            futures.clear();
        }
        std::size_t bounds[hasher_.table_size() + 1];
        bounds[0] = 0;
        if (no_task_split) {
            auto write_itr = ids_beg;
            for(std::size_t i = 0; i < hasher_.table_size(); i ++ ){
                auto& b = buckets[i];
                for(auto&& n : b) {
                    // std::cout << "buk[" << i << "]" << "-->" << n << std::endl; 
                    *write_itr = n;
                    write_itr ++;
                }
                bounds[i + 1] = bounds[i] + b.size();
            }
        } else {
            for(std::size_t i = 0; i < hasher_.table_size(); i ++ ) {
                auto* pb = &(buckets[i]);
                bounds[i + 1] = bounds[i] + pb->size();
                futures.push_back(executor.submit([&bounds, i, &ids_beg, &ids_end, pb]() {
                    auto beg = ids_beg + bounds[i];
                    for(auto&& n : *pb) {
                        *beg = n;
                        beg ++;
                    }
                    return true;
                }));
            }
            for(auto&& fut : futures) {
                fut.sync();
            }
            futures.clear();
        }
        bucket_using = false;

        // for(std::size_t i = 0; i < bounds[i]; i ++ ) {
        //     std::cout << "bounds[" << i << "]: " << bounds[i] << std::endl;
        // }

        if(no_task_split){
            for(std::size_t i = 0; i < hasher_.table_size(); i ++ )
            {
                auto split_beg = ids_beg + bounds[i];
                auto split_end = ids_beg + bounds[i + 1];
                core(executor, split_beg, split_end, str_rng, level + 1);
            }
        } else {
            
            for(std::size_t i = 0; i < hasher_.table_size(); i ++ )
            {
                auto split_beg = ids_beg + bounds[i];
                auto split_end = ids_beg + bounds[i + 1];
                // if(split_beg < ids_beg || split_end > ids_end)
                //     throw std::runtime_error("BUG !!!");
                local_futs.push_back(executor.submit([
                    &executor, split_beg, split_end, 
                    &str_rng, &level, this
                ](){
                    return core(executor, split_beg, split_end, str_rng, level + 1);
                }));
            }
            std::size_t not_done = local_futs.size();
            while(not_done > 0) {
                for(auto&& [i, fut] : local_futs | nucleona::range::indexed()) {
                    if(!fut.valid()) continue;
                    if(!fut.sync()) {
                        // std::cout << "a fail job detected" << std::endl;
                        auto split_beg = ids_beg + bounds[i];
                        auto split_end = ids_beg + bounds[i + 1];
                        fut = executor.submit([
                            &executor, split_beg, split_end, 
                            &str_rng, &level, this
                        ](){
                            return core(
                                executor, split_beg, split_end, str_rng, level + 1
                            );
                        });
                    } else {
                        not_done --;
                    }
                }
            }
            local_futs.clear();
        }
        return true;
    }
    template<class Executor, class StrRng>
    auto operator()(
        Executor&            executor, 
        std::vector<IdxInt>& ids, 
        StrRng&&             str_rng,
        int                  level = 0
    ) const {
        // std::cout << "sort range: " 
        //     << "[0," << ids.size() << ")\n";
        core(executor, ids.begin(), ids.end(), str_rng, level);
    }
    template<class Executor, class StrRng>
    auto operator()(
        Executor&           executor, 
        GroupFile&          group_file,
        StrRng&&            str_rng
    ) const {
        std::ifstream fin(group_file.filename, std::ios::binary);
        std::vector<IdxInt> ids(group_file.num);
        // std::cout << "read file: " << group_file.filename << std::endl;
        for(auto&& i : nucleona::range::irange_0(group_file.num)) {
            fin.read((char*)&ids[i], sizeof(IdxInt));
            // std::cout << "--->" << std::dec << ids[i] << std::endl;
        }
        operator()(executor, ids, str_rng, group_file.level + 1); // TODO: verify this
        fin.close();
        boost::filesystem::remove(group_file.filename);
        return ids;
    }
    template<class Executor, class Str>
    auto operator()(
        Executor&           executor, 
        Str&&               str
    ) const {
        auto strs = nucleona::range::irange_0(
                IdxInt(str.size())
            )
            | nucleona::range::transformed([&str](auto&& idx){
                return make_suffix(str, idx);
            })
        ;
        auto ids = nucleona::range::irange_0(
                IdxInt(str.size())
            )
            | ranges::to_vector
        ;
        operator()(executor, ids, strs);
        return ids;
    }
private:
    Hasher      hasher_;
    std::size_t max_level_;
    // static thread_local std::unique_ptr<
    //    std::vector<IdxInt>[]
    // > buckets;
    // static thread_local std::unique_ptr<
    //     std::mutex[]
    // > bucket_muxs;
};

// template<class Hasher, class IdxInt>
// thread_local std::unique_ptr<
//     std::vector<IdxInt>[]
// > HashMSDRadixSort<Hasher, IdxInt>::buckets;
// 
// template<class Hasher, class IdxInt>
// thread_local std::unique_ptr<
//     std::mutex[]
// > HashMSDRadixSort<Hasher, IdxInt>::bucket_muxs;

}
