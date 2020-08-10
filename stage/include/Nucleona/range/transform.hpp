#pragma once
#include <tuple>
#include <Nucleona/tuple.hpp>
#include <boost/iterator_adaptors.hpp>
#include <Nucleona/range/core.hpp>
#include <Nucleona/mpl/integer_sequence.hpp>
#include <Nucleona/range/iterator_tag_order.hpp>
namespace nucleona{ namespace range{

template< class CORE > struct TransformIterator; // forward declartion

template< class CORE >
using TransformIteratorProto = boost::iterator_adaptor<
      TransformIterator<CORE>
    , typename CORE::BaseIterator
    , typename CORE::Value
    , typename CORE::IteratorCategory
    , typename CORE::Reference
    , typename CORE::Difference
>;

template<class RNG, class FUNC, class TAG>
struct TransformIteratorCore
{
    RNG base_range;
    FUNC function;
    using BaseRange         = RNG;
    using Function          = FUNC;
    using BaseIterator      = typename std::remove_reference_t<RNG>::iterator;
    using Reference         = decltype(function(
        std::declval<
            typename std::iterator_traits< 
                BaseIterator 
            >::reference
        >()
    ));
    using IteratorCategory  = LowerOrder< 
          typename std::iterator_traits< BaseIterator >::iterator_category
        , TAG
    >;
    using Value             = typename std::remove_reference<Reference>::type; // TODO care this type
    using Difference        = typename std::iterator_traits<BaseIterator>::difference_type;
};

template< class CORE >
struct TransformIterator : public TransformIteratorProto< CORE >
{
private:
    friend class boost::iterator_core_access;
    using Base          = TransformIteratorProto< CORE >;
    using Function      = typename CORE::Function;
    using BaseRange     = typename CORE::BaseRange;
    typename Base::reference dereference() const
    {
        return (*func_) ( *this->base_reference() );
    }
    Function* func_ { nullptr };
public:
    TransformIterator( CORE& core )
    : Base  ( core.base_range.begin()   )
    , func_ ( &core.function            )
    {}
    TransformIterator( nucleona::range::IteratorEnd, CORE& core )
    : Base  ( core.base_range.end() )
    , func_ ( nullptr )
    {}

    TransformIterator(){}

    template<class RNG, class FUNC, class TAG>
    static auto make_range_core( RNG&& rng, FUNC&& func, TAG tag )
    {
        return TransformIteratorCore< RNG, FUNC, TAG >{ FWD(rng), FWD(func) };
    }
};

constexpr struct Transform
{
    template<class RNG, class FUNC, class TAG = std::random_access_iterator_tag>
    decltype(auto) make_transform_range( RNG&& rng, FUNC&& func, TAG tag = std::random_access_iterator_tag() ) const 
    {
        return nucleona::range::make< 
            TransformIterator<
                TransformIteratorCore< RNG, FUNC, TAG >
            >
        >( FWD( rng ), FWD( func ), tag ); 
    }
    template<class RNG, class FUNC, class TAG = std::random_access_iterator_tag>
    decltype(auto) operator()(RNG&& rng, FUNC&& func, TAG tag = std::random_access_iterator_tag{} ) const
    {
        return make_transform_range( FWD(rng), FWD(func), tag );
    } 
} transform;

constexpr struct TransformPipeOp
{
    template<class FUNC, class TAG = std::random_access_iterator_tag>
    decltype(auto) operator()( FUNC&& func, TAG tag = std::random_access_iterator_tag() ) const
    {
        return nucleona::make_tuple( transform, FWD(func), tag );
    }
} transformed;

}}
