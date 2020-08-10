#pragma once
#include <iterator>
#include <boost/iterator.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/map.hpp>
#include <boost/range/iterator.hpp>
#include <boost/iterator/iterator_categories.hpp>
#include <boost/mpl/comparison.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/at_fwd.hpp>
#include <boost/mpl/less_equal.hpp>

namespace nucleona{
namespace range{
using IteratorTagOrder = boost::mpl::map<
      boost::mpl::pair< std::input_iterator_tag             , boost::mpl::int_<0> >
    , boost::mpl::pair< std::output_iterator_tag            , boost::mpl::int_<0> >
    , boost::mpl::pair< boost::single_pass_traversal_tag    , boost::mpl::int_<0> >
    , boost::mpl::pair< std::forward_iterator_tag           , boost::mpl::int_<1> >
    , boost::mpl::pair< boost::forward_traversal_tag        , boost::mpl::int_<1> >
    , boost::mpl::pair< std::bidirectional_iterator_tag     , boost::mpl::int_<2> >
    , boost::mpl::pair< boost::bidirectional_traversal_tag  , boost::mpl::int_<2> >
    , boost::mpl::pair< std::random_access_iterator_tag     , boost::mpl::int_<3> >
    , boost::mpl::pair< boost::random_access_traversal_tag  , boost::mpl::int_<3> >
>;
namespace detail{
template<class TAG1, class TAG2>
struct LowerOrderProto
{
    using T1N = typename boost::mpl::at< IteratorTagOrder, TAG1 >::type;
    using T2N = typename boost::mpl::at< IteratorTagOrder, TAG2 >::type;
    using LE12 = typename boost::mpl::less_equal< T1N, T2N >::type;
    constexpr static bool le = LE12::value;
    using Result
        = std::conditional_t<
              le
            , TAG1
            , TAG2
        >
    ;
};
}

template<class TAG1, class TAG2>
using LowerOrder = typename detail::LowerOrderProto<TAG1, TAG2>::Result;

}
}



