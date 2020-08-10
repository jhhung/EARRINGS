#pragma once
#include "common.h"
#include <Nucleona/tuple.hpp>
#include <Nucleona/algo/segment.hpp>
namespace nucleona::range {

template<class DatRng, class SegDesRng>
struct SegmentView 
: ranges::view_facade<
    SegmentView<DatRng, SegDesRng>,
    ranges::range_cardinality<SegDesRng>::value
>
{
    using This = SegmentView<DatRng, SegDesRng>;
    SegmentView() = default;
    SegmentView(DatRng&& d, SegDesRng&& i) 
    : dat_rng_ (std::forward<DatRng>(d))
    , seg_des_rng_ (std::forward<SegDesRng>(i))
    {}
private:
    friend struct ranges::range_access;
    using SegDesConstIter = typename std::decay_t<SegDesRng>::const_iterator;
    template<class SV>
    struct Cursor {
        // using Reference         = typename std::decay_t<DatRng>::reference;
        // using ConstReference    = typename std::decay_t<DatRng>::const_reference;
        Cursor() = default;
        Cursor( SV& rng, SegDesConstIter it )
        : rng_          (&rng)
        , seg_des_itr   (it)
        {}
        DEFAULT_COPY(Cursor)
        DEFAULT_MOVE(Cursor)

        decltype(auto) read() const{
            auto& dat_rng = rng_->dat_rng_.get();
            auto& seg_des = *seg_des_itr; 
            return ranges::view::slice(dat_rng, seg_des.pos, seg_des.pos + seg_des.size);
        }
        decltype(auto) read() {
            auto& dat_rng = rng_->dat_rng_.get();
            auto& seg_des = *seg_des_itr; 
            return ranges::view::slice(dat_rng, seg_des.pos, seg_des.pos + seg_des.size);
        }
        bool equal( const Cursor& c) const {
            return seg_des_itr == c.seg_des_itr;
        }

        void next() {
            ++seg_des_itr;
        }

        void prev() {
            --seg_des_itr;
        }

        std::ptrdiff_t distance_to( const Cursor& c ) const {
            return c.seg_des_itr - seg_des_itr;
        }

        void advance(std::ptrdiff_t n ) {
            seg_des_itr += n;
        }

        SV*                rng_            ;
        SegDesConstIter    seg_des_itr     ;
    };
    Cursor<This> begin_cursor() {
        return {*this, seg_des_rng_.begin()};
    }
    Cursor<This> end_cursor() {
        return {*this, seg_des_rng_.end()};
    }
    Cursor<const This> begin_cursor() const {
        return {*this, seg_des_rng_.begin()};
    }
    Cursor<const This> end_cursor() const {
        return {*this, seg_des_rng_.end()};
    }
    
    nucleona::util::AutoRefWrapper<DatRng> dat_rng_;
    SegDesRng seg_des_rng_;
};
struct SegmentFn
{
    template<
        class Rng, 
        class SegNum,
        CONCEPT_REQUIRES_(ranges::InputRange<Rng>())
    >
    auto operator()(
        Rng&& rng,
        SegNum&& seg_num
    ) const {
        auto size = ranges::distance(rng);
        auto segs = algo::segment(size, seg_num);
        return SegmentView<Rng, decltype(segs)>(
            FWD(rng), std::move(segs)
        );
        // return ranges::view::transform(std::move(segs), [
        //     _rng = std::forward<Rng>(rng)
        // ](auto&& seg){
        //     return ranges::view::slice(_rng, seg.pos, seg.pos + seg.size);
        // });
    }
    template<class SegNum>
    auto operator()(SegNum seg_num) const {
        return ranges::make_pipeable([
            ___this = this, seg_num
        ](auto&& rng){
            return (*___this)(FWD(rng), seg_num);
        });
    }
};
RANGES_INLINE_VARIABLE(SegmentFn, segment)

}