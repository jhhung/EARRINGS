#pragma once
#include "common.h"
#include <Nucleona/parallel/thread_pool.hpp>
#include <Nucleona/range/utils.hpp>
#include "segment.hpp"
#include "endp.hpp"
#include <Nucleona/functional/mutable.hpp>

namespace nucleona::range {

struct SPEndPFn
{
    template<
        class Rng, 
        class Executor,
        CONCEPT_REQUIRES_(ranges::InputRange<Rng>())
    >
    void operator()(
        Rng&& rng, 
        Executor& thread_pool
    ) const {
        std::vector<typename Executor::template Future<void>> fut_vec;
        fut_vec.reserve(thread_pool.get_parallization());
        for(auto&& pack : rng | nucleona::range::segment(
            thread_pool.get_parallization() + (thread_pool.support_owner_as_worker() ? 1 : 0)
        )) {
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
        Executor& thread_pool
    ) const {
        return ranges::make_pipeable([
            &thread_pool,
            ___this = this
        ](auto&& rng){
            return (*___this)(
                FWD(rng),
                thread_pool
            );
        });
    }
};
RANGES_INLINE_VARIABLE(SPEndPFn, sp_endp)
}