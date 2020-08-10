#pragma once
#include <tuple>
#include <Nucleona/tuple.hpp>
#include <boost/iterator_adaptors.hpp>
#include <Nucleona/range/core.hpp>
#include <Nucleona/mpl/integer_sequence.hpp>
#include <Nucleona/util/auto_ref_wrapper.hpp>
namespace nucleona{ namespace range{
namespace detail{
template<class T>
using FilterStorage = nucleona::util::AutoRefWrapper<T>;

struct FilterParameter
{
    template< class RNG, class FUNC >
    auto requires( RNG r, FUNC func ) NUCLEONA_EXPRS(
        (bool)func(*r.begin())
    );
};

}
template< class CORE > struct FilterIterator; // forward declartion

template< class CORE >
using FilterIteratorProto = boost::iterator_adaptor<
      FilterIterator<CORE>
    , typename CORE::BaseIterator
    , typename CORE::Value
    , typename CORE::IteratorCategory
    , typename CORE::Reference
    , typename CORE::Difference
>;
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
template<class T> 
struct FilterReference
{
    // if a reference come, then return reference, if a value come, then return const reference
    using Result = const T&;
};
template<class T>
struct FilterReference<T&>
{
    using Result = T&;
};
template<class RNG, class FUNC>
struct FilterIteratorCore
{
    RNG     base_range;
    FUNC    function;
    using BaseRange         = RNG;
    using Function          = FUNC;
    using BaseIterator      = typename std::remove_reference_t<RNG>::iterator;
    using BaseReference     = typename std::iterator_traits< BaseIterator >::reference;
    using Reference         = typename FilterReference< BaseReference >::Result;
    using IteratorCategory  = typename IteratorDegrade< 
        typename std::iterator_traits< BaseIterator >::iterator_category 
    >::Result; 
    using Value             = typename std::iterator_traits< BaseIterator >::value_type;
    using Difference        = typename std::iterator_traits< BaseIterator >::difference_type;
};
template< class CORE >
struct FilterIterator 
: public FilterIteratorProto< CORE >
, public detail::FilterStorage< typename CORE::BaseReference >
{
  private:
    friend class boost::iterator_core_access;
    using Base          = FilterIteratorProto< CORE >;
    using Function      = typename CORE::Function;
    using BaseRange     = typename CORE::BaseRange;
    using Storage       = detail::FilterStorage< typename CORE::BaseReference >;
    typename Base::reference dereference() const
    {
        return const_cast<typename Base::reference>(Storage::get());
    }

    void increment()
    {
        do
        {
            this->base_reference() ++;
            Storage::set( *this->base_reference() );
        }
        while( 
               !(*func_)( Storage::get() ) 
            && this->base_reference() != end_itr
        );
    }

    void decrement()
    {
        do
        {
            this->base_reference() --;
            Storage::set( *this->base_reference() );
        }
        while( 
               !(*func_)( Storage::get() ) 
            && this->base_reference() >= begin_itr
        );
    }

    Function* func_ { nullptr };
    typename CORE::BaseIterator end_itr;
    typename CORE::BaseIterator begin_itr;
  public:
    FilterIterator( CORE& core )
    : Base      ( core.base_range.begin()   )
    , Storage   ( *this->base_reference()   )
    , func_     ( &core.function            )
    , end_itr   ( core.base_range.end()     )
    , begin_itr ( core.base_range.begin()   ) 
    {
        while(true)
        {
            if ( (*func_)( Storage::get() ) )
                break;
            if ( this->base_reference() == end_itr )
                break;
            this->base_reference()++;
            Storage::set( *this->base_reference() );
        }
    }
    FilterIterator( nucleona::range::IteratorEnd, CORE& core )
    : Base      ( core.base_range.end()     )
    , func_     ( nullptr                   )
    , end_itr   ( core.base_range.end()     )
    , begin_itr ( core.base_range.begin()   ) 
    {}

    FilterIterator(){}

    template<class RNG, class FUNC, NUCLEONA_CONCEPT_REQUIRE_(detail::FilterParameter, RNG, FUNC)>
    static auto make_range_core( RNG&& rng, FUNC&& func )
    {
        return FilterIteratorCore< RNG, FUNC >{ FWD(rng), FWD(func) };
    }
};


constexpr struct Filter
{
    template<class RNG, class FUNC>
    decltype(auto) make_filter_range( RNG&& rng, FUNC&& func ) const 
    {
        return nucleona::range::make< 
            FilterIterator<
                FilterIteratorCore< RNG, FUNC>
            >
        >( FWD( rng ), FWD( func ) ); 
    }

    template<class RNG, class FUNC, NUCLEONA_CONCEPT_REQUIRE_(detail::FilterParameter, RNG, FUNC)>
    decltype(auto) operator()(RNG&& rng, FUNC&& func) const
    {
        return make_filter_range( FWD(rng), FWD(func) );
    } 
} filter;

constexpr struct FilterPipeOp
{
    template<class FUNC>
    decltype(auto) operator()( FUNC&& func ) const
    {
        return nucleona::make_tuple( filter, FWD(func) );
    }
} filtered;

}}
