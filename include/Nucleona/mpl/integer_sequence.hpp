#pragma once
#include <type_traits>
#include <utility>
namespace nucleona{ namespace mpl{

namespace detail {

template<class END, END end>
struct BoundEnd
{
    template<class INT, INT ini >
    constexpr static bool value = ini < end;
};
template<class STEP, STEP step>
struct StepNext
{
    template<class INT, INT ini >
    constexpr static INT value = ini + step;
};

}

template<class INT, INT ini, class COND, class NEXT, bool b, INT... n>
struct MakeIndexForProto
{
};
template<class INT, INT ini, class COND, class NEXT, INT... n>
struct MakeIndexForProto< INT, ini, COND, NEXT, true, n... >
{
    constexpr static INT next = NEXT:: template value< INT, ini>;
    using Result 
        = typename MakeIndexForProto<
              INT
            , next
            , COND
            , NEXT
            , COND:: template value< INT, next>
            , n...
            , ini
        >::Result
    ;
};
template<class INT, INT ini, class COND, class NEXT, INT... n>
struct MakeIndexForProto< INT, ini, COND, NEXT, false, n... >
{
    using Result 
        = std::integer_sequence<
              INT
            , n...
        >
    ;

};
template<class INT, INT ini, class COND, class NEXT>
using MakeIndexFor = typename MakeIndexForProto<
      INT
    , ini
    , COND
    , NEXT
    , COND:: template value<INT, ini>
>::Result;

template<std::size_t start, std::size_t end, std::size_t step = 1>
constexpr auto make_index_for() 
{
    return MakeIndexFor< 
          std::size_t
        , start
        , detail::BoundEnd<std::size_t, end>
        , detail::StepNext<std::size_t, step>
    >{};
}

}}
