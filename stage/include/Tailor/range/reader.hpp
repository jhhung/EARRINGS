#pragma once

#include <Tailor/range/get_fastq_view.hpp>
#include <Tailor/range/get_fasta_view.hpp>
#include <Tailor/range/fastq_reader.hpp>
#include <Tailor/tailor/tailor_main.hpp>
#include <istream>
#include <unordered_map>
#include <functional>

namespace pipeline::range {

template <class DataType>
class Reader {};

template <>
class Reader <tailor::Fastq>
{
  public:
  std::unordered_map<
    tailor::Fastq
  , std::uint32_t
  , pipeline::range::Hash<tailor::Fastq>
  , pipeline::range::Equal<tailor::Fastq>
  > buf;
 
  Reader(){}; 

  explicit Reader(std::istream& input)
  {
    for (auto&& i : input | pipeline::range::get_fastq(buf)) {} 
  }  

  auto& get_buf() noexcept
  {
    return buf;
  }
};

template <>
class Reader <tailor::Fasta>
{
  public:
  std::unordered_map<
    tailor::Fasta
  , std::uint32_t
  , pipeline::range::Hash<tailor::Fasta>
  , pipeline::range::Equal<tailor::Fasta>
  > buf;
 
  Reader(){}; 

  explicit Reader(std::istream& input)
  {
    for (auto&& i : input | pipeline::range::get_fasta(buf)) {} 
  }  

  auto& get_buf() noexcept
  {
    return buf;
  }
};

}
