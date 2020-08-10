
#pragma once
#include <range/v3/all.hpp>

namespace pipeline::range {

auto comp = [](auto&& compared, auto&& smallest) {
  if (compared.tail_pos_ > smallest.tail_pos_)
    return false;
  else if (compared.tail_pos_ < smallest.tail_pos_)
    return true;
  else 
    return compared.flag_ == 0; 
};

template <class Fun>
auto merge_results(Fun comp)
{
  /*
  return ranges::view::filter ([](auto&& AlignedReads_v){ return !AlignedReads_v.empty(); }) 
       | ranges::view::transform (
           [comp](auto&& AlignedReads_v) {
             auto iter = std::min_element (
               AlignedReads_v.cbegin(),
               AlignedReads_v.cend(),
               comp
             );
             return *iter;  
         });
  */
  return ranges::view::transform (
    [comp](auto&& AlignedReads_v) { 
      using AlignedReads = typename std::remove_reference_t<decltype(AlignedReads_v)>::value_type;
      using Fastq = typename AlignedReads::fastq;
      
      if (AlignedReads_v.empty())
        return AlignedReads(Fastq{});
      auto iter = std::min_element (
        AlignedReads_v.cbegin(),
        AlignedReads_v.cend(),
        comp
      );
      return *iter;
  });
}
}
