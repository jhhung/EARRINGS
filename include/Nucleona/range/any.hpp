#pragma once
#include <memory>
#include <boost/range/any_range.hpp>
#include <Nucleona/language.hpp>
#include <Nucleona/tuple.hpp>

namespace nucleona{ namespace range{
namespace detail{

template<class RNG>
struct AnyRangeStorage
{
    AnyRangeStorage( RNG&& rng )
    : rng_ ( FWD( rng ) )
    {}
  protected:
    RNG rng_;
};
template<
    class Value
  , class Traversal
  , class Reference
  , class Difference
  , class Buffer = boost::any_iterator_default_buffer
>
struct AnyInterface
{
    using iterator = boost::range_detail::any_iterator<
        Value
      , Traversal
      , Reference
      , Difference
      , typename boost::range_detail::any_range_default_help<
            Buffer
          , boost::mpl::identity<boost::any_iterator_default_buffer>
        >::type
    >;
    using const_iterator = iterator;
    using reference = Reference;
    using difference_type = Difference;

    virtual iterator begin()        = 0;
    virtual iterator end()          = 0;
    // virtual const_iterator begin() const  = 0;
    // virtual const_iterator end()   const  = 0;
};
template<class RNG, class INTERFACE>
struct AnyImpl 
: public AnyRangeStorage< RNG >
, public INTERFACE
{
    using Interface = INTERFACE;
    using iterator          = typename Interface::iterator;
    using const_iterator    = typename Interface::const_iterator;
    AnyImpl( RNG&& rng )
    : AnyRangeStorage< RNG > ( FWD( rng ) )
    {}
    virtual iterator begin()    override { return iterator( this->rng_.begin() ); }
    virtual iterator end()      override { return iterator( this->rng_.begin() ); }

};

}
template<
    class Value
  , class Traversal
  , class Reference
  , class Difference
  , class Buffer = boost::any_iterator_default_buffer
>
struct AnyWrapper
{
    using Interface = detail::AnyInterface< 
          Value
        , Traversal
        , Reference
        , Difference
        , Buffer
    >;
    using iterator          = typename Interface::iterator; 
    using const_iterator    = typename Interface::const_iterator; 
    using reference         = typename Interface::reference;

    template<class RNG >
    AnyWrapper( RNG&& rng )
    :  impl_( new detail::AnyImpl<RNG, Interface>( FWD(rng) ) )
    {}
    // auto begin()    const { return impl_->begin();  }
    // auto end()      const { return impl_->end();    }
    auto begin()    { return impl_->begin();  }
    auto end()      { return impl_->end();    }
    std::unique_ptr<detail::AnyInterface<
          Value
        , Traversal
        , Reference
        , Difference
        , Buffer
    >> impl_;
};

template<class RNG>
auto make_any_range( RNG&& rng )
{
    return AnyWrapper<
          typename RNG::value_type
        , typename boost::iterator_traversal<
            typename boost::range_iterator<RNG>::type 
        >::type
        , typename RNG::reference
        , typename RNG::difference_type
    >( FWD(rng) );
}
constexpr struct DefaultEraseType
{
    template<class RNG>
    decltype(auto) operator()( RNG&& rng ) const 
    {
        return AnyWrapper<
              typename RNG::value_type
            , typename boost::iterator_traversal<
                typename boost::range_iterator<RNG>::type 
            >::type
            , typename RNG::reference
            , typename RNG::difference_type
        >( FWD(rng) );
        // return make_any_range( FWD( rng ) );
    }
} default_erase_type;
template<
      class Value       
    , class Traversal   
    , class Reference   
    , class Difference  
    , class Buffer      
>
struct EraseType
{
    template<class RNG>
    decltype(auto) operator()( RNG&& rng ) const 
    {
        return AnyWrapper<
              typename RNG::value_type
            , typename boost::iterator_traversal<
                typename boost::range_iterator<RNG>::type 
            >::type
            , typename RNG::reference
            , typename RNG::difference_type
        >( FWD(rng) );
    }
}; 

template<
      class Value       
    , class Traversal   
    , class Reference   
    , class Difference  
    , class Buffer      
>
constexpr EraseType<
      Value       
    , Traversal   
    , Reference   
    , Difference  
    , Buffer      
> erase_type;

template<
      class Value       
    , class Traversal   
    , class Reference   
    , class Difference  
    , class Buffer      
>
struct EraseTypePipeOp
{
    decltype(auto) operator()( ) const 
    {
        return nucleona::make_tuple( erase_type<
              Value       
            , Traversal   
            , Reference   
            , Difference  
            , Buffer      
        > );    
    }
};

template<
      class Value       
    , class Traversal   
    , class Reference   
    , class Difference  
    , class Buffer      
>
constexpr EraseTypePipeOp<
      Value       
    , Traversal   
    , Reference   
    , Difference  
    , Buffer      
> type_erased;


}

using range::make_any_range;

}
