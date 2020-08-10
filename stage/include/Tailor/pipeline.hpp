#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <mutex>
#include <map>
#include <atomic>
#include <Tailor/tailor.hpp>
#include <range.hpp>
#include <Nucleona/range/v3_impl.hpp>
#include <Nucleona/parallel/thread_pool.hpp>
#include <unistd.h>


namespace pipeline {


// template <class Str, class TailorMain>
// void pipeline(Str&& input_name, TailorMain&& tailor_mapping, std::ostream& out)
// {
//   using trueType = std::bool_constant<true>;
//   using falseType = std::bool_constant<false>;
//   std::ifstream input(input_name);  
//   if (isFasta)
//   {
//     pipeline::range::Reader<falseType> reader(input);
//     input.close();
//     auto& buf = reader.get_buf();
//     pipeline_exec(tailor_mapping, buf, out);
//   }
//   else
//   {
//     pipeline::range::Reader<trueType> reader(input);
//     input.close();
//     auto& buf = reader.get_buf();
//     pipeline_exec(tailor_mapping, buf, out);
//   }
// }

template <class Str, class TailorMain>
void pipeline(Str&& input_name, TailorMain&& tailor_mapping, std::ostream& out)
{
  const auto& aligner = tailor_mapping.get_table();
  const auto& paras = tailor_mapping.get_paras();
  auto tp = nucleona::parallel::make_thread_pool(paras.n_thread);
  using IndexType = typename std::remove_const_t<std::remove_reference_t<decltype(aligner)>>::IndexType;
  using DataType = typename std::remove_const_t<std::remove_reference_t<decltype(aligner)>>::FASTQ;
  constexpr pipeline::range::format_reader_fn<DataType> format_reader{};  
  std::ifstream input(input_name);  
  // pipeline::range::Reader<DataType> reader(input);
  // input.close();
  // auto& buf = reader.get_buf();

  // buf 
  // | pipeline::range::align(aligner)
  // | ranges::view::transform(
  //     [](auto&& AlignedReads_v) {
  //       for (auto&& i : AlignedReads_v)
  //       {
  //         i.print();
  //       }
  //       return AlignedReads_v;
  //     }
  //   ) 
  // | nucleona::range::p_endp(tp);
  
  input
  | format_reader()
  | pipeline::range::align(aligner)
  | ranges::view::transform(
      [](auto&& AlignedReads_v) {
        for (auto&& i : AlignedReads_v)
        {
          i.print();
        }
        return AlignedReads_v;
      }
    ) 
  | nucleona::range::endp;
}

}
