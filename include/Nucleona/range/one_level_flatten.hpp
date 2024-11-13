#pragma once
#include <Nucleona/range/core.hpp>
#include <Nucleona/util/auto_ref_wrapper.hpp>
#include <Nucleona/tuple.hpp>
#include <iostream>

namespace nucleona{ namespace range{
template< class CORE > struct OneLevelFlattenIterator;

namespace detail{

template<class T>
struct IteratorDegrade
{
    using Result = T;
};
template<>
struct IteratorDegrade< std::random_access_iterator_tag >
{
    using Result = std::bidirectional_iterator_tag;
};
template<>
struct IteratorDegrade< boost::random_access_traversal_tag >
{
    using Result = boost::bidirectional_traversal_tag;
};
template< class RNG >
struct OneLevelFlattenIteratorCore
{
    using BaseRange         = RNG;
    using BaseIterator      = typename std::remove_reference_t<RNG>::iterator;
    using InnerRangeRef     = typename std::iterator_traits<BaseIterator>::reference;
    using InnerRange        = std::remove_reference_t< InnerRangeRef >;
    using InnerIterator     = typename InnerRange::iterator;
    using InnerRangeAutoRef = nucleona::util::AutoRef< InnerRangeRef >;
    using AdaptIterator     = BaseIterator;
    using Reference         = typename std::iterator_traits< InnerIterator >::reference;
    using Value             = typename std::iterator_traits< InnerIterator >::value_type;
    using Difference        = typename std::iterator_traits< InnerIterator >::difference_type;
    using IteratorCategory  = boost::single_pass_traversal_tag;
    //     = typename IteratorDegrade< 
    //     typename std::iterator_traits< AdaptIterator >::iterator_category 
    // >::Result; 

    BaseRange base_range;
    // std::remove_reference_t<InnerRange> inner_range;

    // OneLevelFlattenIteratorCore( RNG&& rng )
    // : base_range ( FWD(rng) )
    // , inner_range( base_range.
};

template< class CORE >
using OneLevelFlattenIteratorProto = boost::iterator_adaptor<
      OneLevelFlattenIterator<CORE>
    , typename CORE::AdaptIterator
    , typename CORE::Value
    , typename CORE::IteratorCategory
    , typename CORE::Reference
    , typename CORE::Difference
>;

struct OneLevelFlattenParameter
{
    template< class RNG >
    auto requires_( RNG r ) NUCLEONA_EXPRS(
        r.begin()->begin()
    );
};


}
template< class CORE >
struct OneLevelFlattenIterator 
: public detail::OneLevelFlattenIteratorProto< CORE >
{
  private:
    friend class ::boost::iterator_core_access;
    using Base              = detail::OneLevelFlattenIteratorProto< CORE >      ;
    using InnerRangeAutoRef = typename CORE::InnerRangeAutoRef                  ;
    using InnerIterator     = typename CORE::InnerIterator                      ;
    using BaseRange         = std::decay_t<typename CORE::BaseRange>            ;

    InnerRangeAutoRef   inner_range_ref_;
    InnerIterator       inner_iter_;
    BaseRange*          p_base_range_;

    typename Base::reference dereference() const
    {
        return *inner_iter_;
    }
    void increment()
    {
        inner_iter_ ++;
        if ( inner_range_ref_.get().end() == inner_iter_ )
        {
            this->base_reference() ++;
            if ( this->base_reference() != p_base_range_->end() )
            {
                inner_range_ref_.set( *this->base_reference() );
                auto& inner_rng = inner_range_ref_.get();
                inner_iter_ = inner_rng.begin();
            }
            else
            {
                inner_range_ref_.set( nullptr );
                inner_iter_ = InnerIterator();
            }
        }
    }
    bool equal( const OneLevelFlattenIterator<CORE>& iter ) const
    {
        return 
               this->base_reference() == iter.base_reference() 
            && nucleona::util::equal( inner_range_ref_, iter.inner_range_ref_ )
        ;
    }

  public:
    OneLevelFlattenIterator( CORE& core )
    : Base              ( core.base_range.begin()           )
    , inner_range_ref_  ( FWD( *core.base_range.begin() )   )
    , inner_iter_       ( inner_range_ref_.get().begin()    )
    , p_base_range_     ( &(core.base_range)                )
    {}

    OneLevelFlattenIterator( nucleona::range::IteratorEnd, CORE& core )
    : Base              ( core.base_range.end()             )
    , inner_range_ref_  ( nullptr                           )
    , p_base_range_     ( &core.base_range                  )
    {}
    template< class RNG >
    static auto make_range_core( RNG&& rng )
    {
        return detail::OneLevelFlattenIteratorCore< RNG >{ FWD(rng) };
    }
};
constexpr struct OneLevelFlatten
{
    template< class RNG >
    decltype(auto) make_one_level_flatten_range( RNG&& rng ) const 
    {
        return nucleona::range::make< 
            OneLevelFlattenIterator<
                detail::OneLevelFlattenIteratorCore< RNG >
            >
        >( FWD( rng ) ); 
    }

    template<class RNG, NUCLEONA_CONCEPT_REQUIRE_(detail::OneLevelFlattenParameter, RNG )>
    decltype(auto) operator()(RNG&& rng ) const
    {
        // return int();
        return make_one_level_flatten_range( FWD(rng) );
    } 
} one_level_flatten;

constexpr struct OneLevelFlattenPipeOp
{
    decltype(auto) operator()( ) const
    {
        return nucleona::make_tuple( one_level_flatten );
    }
} one_level_flattened;

}}
