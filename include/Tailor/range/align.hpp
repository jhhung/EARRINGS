#pragma once

#include <range/v3/all.hpp>
#include <type_traits>

namespace pipeline::range {

template <class Table>
auto align(Table&& table)
{
  return ranges::view::transform (
    [&table](auto&& query)
    {
      if constexpr (
        std::is_same_v<
          std::remove_const_t<std::remove_reference_t<decltype(query)>>
        , typename std::remove_reference_t<decltype(table)>::FASTQ
        >
      )
      {
        return table.align(query);
      }
      else
      {
        return table.align(query.first, query.second);
      }
    }
  );
}

}

