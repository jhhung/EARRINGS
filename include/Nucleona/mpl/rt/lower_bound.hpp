#pragma once
#include <boost/mpl/int.hpp>
#include <type_traits>
#include <boost/mpl/divides.hpp>
#include <boost/mpl/minus.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/arithmetic.hpp>
#include <boost/mpl/comparison.hpp>
#include <utility>
#include <Nucleona/language.hpp>

namespace nucleona { namespace mpl { 
template<class SEQ, int idx >
struct ValueIterator
{
    constexpr static int index = idx;
};
template<class SEQ>
struct LengthImpl
{};
template<class T, T... n>
struct LengthImpl< std::integer_sequence< T, n... > >
{
    constexpr static std::size_t value = sizeof...(n) ;
};

template<class SEQ>
struct ValueBeginImpl
{
    using Result = ValueIterator< SEQ, 0 >;
};
template<class SEQ>
constexpr static std::size_t length = LengthImpl< SEQ >::value;

template<class SEQ>
struct ValueEndImpl
{
    using Result = ValueIterator< SEQ, length< SEQ > - 1>;
};
template<class SEQ>
using ValueBegin = typename ValueBeginImpl<SEQ>::Result;

template<class SEQ>
using ValueEnd = typename ValueEndImpl<SEQ>::Result;


template<class BEGIN, class END>
struct DistanceImpl
{
};

template<class S1, int idx1, class S2, int idx2 >
struct DistanceImpl< ValueIterator< S1, idx1>, ValueIterator<S2, idx2> > 
{
    using Result = boost::mpl::int_< idx2 - idx1 >;
};


template<class BEGIN, class END>
using Distance = typename DistanceImpl<BEGIN, END>::Result;

template<class BEGIN, class END>
constexpr auto distance = DistanceImpl< BEGIN, END >::Result::value;

template<class T1, class T2>
struct AddImpl
{
    using Result = typename boost::mpl::plus< T1, T2 >::type;
};

template<class S, int EN, class N>
struct AddImpl<  ValueIterator< S, EN >, N >
{
    using Result = ValueIterator< S, EN + N::value >;
};
template<class T1, class T2>
using Add = typename AddImpl<T1, T2>::Result;


template<class T1, class T2>
struct DivImpl
{
    using Result = typename boost::mpl::divides< T1, T2 >::type;
};
template<class T1, class T2 >
using Div = typename DivImpl< T1, T2 >::Result;

template<class T1, class T2>
struct SubImpl
{
    using Result = typename boost::mpl::minus< T1, T2 >::type;
};
template<class T1, class T2 >
using Sub = typename SubImpl< T1, T2 >::Result;

template<class S, int i>
struct AtImpl
{
};
template<class T, T n1, T... n, int i >
struct AtImpl< std::integer_sequence< T, n1, n... >, i >
{
    constexpr static T value 
        = AtImpl< std::integer_sequence< T, n... >, i - 1 >::value;
};
template<class T, T n1, T... n >
struct AtImpl< std::integer_sequence< T, n1, n... >, 0>
{
    constexpr static T value = n1;
};
template<class S, int i>
constexpr auto at = AtImpl< S, i >::value ;


template<class MPL_ITER>
struct DerefImpl
{
};
template<class S, int i>
struct DerefImpl< ValueIterator< S, i > >
{
    constexpr static auto value = at<S, i>;
};

template<class ITER >
constexpr auto deref = DerefImpl<ITER>::value;

template<class I1, class I2>
struct GreaterThanImpl
{
    using Result = typename boost::mpl::greater< I1, I2 >::type;
};
template<class I1, class I2>
using GreaterThan = typename GreaterThanImpl< I1, I2 >::Result;


template<class ITER, class N>
using Advance = Add< ITER, N >;

    
namespace rt { 


template<
      class BEGIN
    , class END
    , class COUNT   = Distance< BEGIN, END >
    , class COND    = GreaterThan< 
          COUNT
        , boost::mpl::int_<0> 
    > 
    , class STEP    = Div< COUNT, boost::mpl::int_<2> >
    , class IT      = Advance< BEGIN, STEP >
>
struct LowerBound // first great equal than
{
    template<class VAL, class FUNC>
    inline static void run( const VAL& value, const FUNC& handle_return )
    {
        using ITAddOne = Add< IT, boost::mpl::int_<1> >;
        if ( deref<IT> < value )
        {
            LowerBound< 
                  ITAddOne
                , END 
                , Sub<COUNT, Add< STEP, boost::mpl::int_<1> > >
            >::run( value, handle_return );
        }
        else
        {
            LowerBound<
                  BEGIN
                , END
                , STEP
            >::run( value, handle_return );
        }
    };
};
template<class BEGIN, class END, class COUNT>
struct LowerBound< BEGIN, END, COUNT, boost::mpl::bool_<false>>
{
    template<class VAL, class FUNC>
    inline static void run( const VAL& value, const FUNC& handle_return )
    {
        handle_return.template operator()(BEGIN());
    }
};

template<class BEGIN, class END>
constexpr LowerBound<BEGIN, END> lower_bound;


} /* rt */ 
} /* mpl */ 
} /* nucleona */ 





