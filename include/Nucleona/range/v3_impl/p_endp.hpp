#pragma once
#include "common.h"
#include <Nucleona/parallel/thread_pool.hpp>
#include <Nucleona/range/utils.hpp>
#include "endp.hpp"
#include "sp_endp.hpp"

namespace nucleona::range {

struct PEndPFn
{
    template<
        class Rng, 
        class Executor,
        CONCEPT_REQUIRES_(ranges::InputRange<Rng>())
    >
    void operator()(
        Rng&& rng, 
        Executor& thread_pool,
        std::size_t pack_size
    ) const {
        std::vector<typename Executor::template Future<void>> fut_vec;
        for(auto&& pack : rng | ranges::view::chunk(pack_size)) {
            fut_vec.push_back(thread_pool.submit([
                _pack = nucleona::functional::mms(std::move(pack))
            ](){
                _pack.storage | endp;
            }));
        }
        for(auto&& fut : fut_vec) {
            fut.sync();
        }
    }
    template<class Executor>
    auto operator()(
        Executor& thread_pool,
        std::size_t pack_size
    ) const {
        return ranges::make_pipeable([
            &thread_pool,
            pack_size,
            ___this = this
        ](auto&& rng){
            return (*___this)(
                FWD(rng),
                thread_pool,
                pack_size
            );
        });
    }
    template<class Executor>
    auto operator()(
        Executor& thread_pool
    ) const {
        return sp_endp(thread_pool);
    }
};
RANGES_INLINE_VARIABLE(PEndPFn, p_endp)
}