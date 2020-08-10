#pragma once
#if defined(NUCLEONA_RANGE_USE_V3)
#   include <Nucleona/range/v3_impl.hpp>
#else
#   include <Nucleona/range/core.hpp>
#   include <Nucleona/range/concept.hpp>
#   include <Nucleona/range/adapter.hpp>
#   include <Nucleona/range/pipe_operator.hpp>
#   include <Nucleona/range/transform.hpp>
#   include <Nucleona/range/print.hpp>
#   include <Nucleona/range/endp.hpp>
#   include <Nucleona/range/filter.hpp>
#   include <Nucleona/range/one_level_flatten.hpp>
#   include <Nucleona/range/indexed.hpp>
#   include <Nucleona/range/at.hpp>
#   include <boost/range/adaptor/type_erased.hpp>
#endif

#include <Nucleona/range/utils.hpp>