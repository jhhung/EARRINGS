#pragma once
#include "common.h"
#include <Nucleona/tuple.hpp>
#include <initializer_list>
#include <Nucleona/functional/mutable.hpp>
#include <Nucleona/range/utils.hpp>
#include <Nucleona/util/auto_ref_wrapper.hpp>

namespace nucleona::range {
template<class DatRng, class IdxRng>
struct AtView 
: ranges::view_facade<
    AtView<DatRng, IdxRng>,
    ranges::range_cardinality<IdxRng>::value
>
{
    using This = AtView<DatRng, IdxRng>;
    AtView() = default;
    AtView(DatRng&& d, IdxRng&& i) 
    : dat_rng_ (std::forward<DatRng>(d))
    , idx_rng_ (std::forward<IdxRng>(i))
    {}
private:
    friend struct ranges::range_access;
    using IdxConstIter = typename std::decay_t<IdxRng>::const_iterator;
    template<class AV>
    struct Cursor {
        // using Reference         = typename std::decay_t<DatRng>::reference;
        // using ConstReference    = typename std::decay_t<DatRng>::const_reference;
        Cursor() = default;
        Cursor( AV& rng, IdxConstIter it )
        : rng_      (&rng)
        , idx_itr   (it)
        {}
        DEFAULT_COPY(Cursor)
        DEFAULT_MOVE(Cursor)
        template<class Rng, class Idx>
        static decltype(auto) access(Rng&& rng, Idx idx) {
            return rng[idx];
        }
        template<class K, class V, class Idx>
        static decltype(auto) access(const std::map<K,V>& rng, Idx idx) {
            return rng.at(idx);
        }
        decltype(auto) read() const{
            auto& dat_rng = rng_->dat_rng_.get();
            auto idx = *idx_itr; 
            return access(dat_rng, idx);
        }
        decltype(auto) read() {
            auto& dat_rng = rng_->dat_rng_.get();
            auto idx = *idx_itr; 
            return access(dat_rng, idx);
        }
        bool equal( const Cursor& c) const {
            return idx_itr == c.idx_itr;
        }

        void next() {
            ++idx_itr;
        }

        void prev() {
            --idx_itr;
        }

        std::ptrdiff_t distance_to( const Cursor& c ) const {
            return c.idx_itr - idx_itr;
        }

        void advance(std::ptrdiff_t n ) {
            idx_itr += n;
        }

        AV*             rng_        ;
        IdxConstIter    idx_itr     ;
    };
    Cursor<This> begin_cursor() {
        return {*this, idx_rng_.begin()};
    }
    Cursor<This> end_cursor() {
        return {*this, idx_rng_.end()};
    }
    Cursor<const This> begin_cursor() const {
        return {*this, idx_rng_.begin()};
    }
    Cursor<const This> end_cursor() const {
        return {*this, idx_rng_.end()};
    }
    nucleona::util::AutoRefWrapper<DatRng> dat_rng_;
    IdxRng idx_rng_;
};
struct AtFn
{
    template<class DatRng, class IdxRng>
    auto operator()(DatRng&& dat_rng, IdxRng idx_rng) const {
        return AtView<DatRng, IdxRng>(
            FWD(dat_rng),
            std::move(idx_rng)
        );
    }
    template<class IdxRng>
    auto operator() ( IdxRng&& idx_rng ) const 
    {
        return ranges::make_pipeable([
            _param = nucleona::make_tuple(
                FWD(idx_rng), *this 
            )
        ]( auto&& dat_rng) mutable {
            auto& [_idx_rng, _func] = _param;
            return _func(
                FWD(dat_rng), 
                std::forward<IdxRng>(_idx_rng)
            );
        });
    }
    auto operator() ( const std::initializer_list<int>& idx_rng ) const {
        return operator()(std::vector<int>(idx_rng));
    }
};
RANGES_INLINE_VARIABLE(AtFn, at);
}
RANGES_SATISFY_BOOST_RANGE(::nucleona::range::AtView);