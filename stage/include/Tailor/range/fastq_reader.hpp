#pragma once

#include <Tailor/range/get_fastq_view.hpp>
#include <istream>
#include <unordered_map>
#include <functional>

namespace pipeline::range {

enum ContainerType
{
  MAP,
  VEC,
  UNORDERED_MAP
};

template <class Fastq>
struct Hash
{
  auto operator()(const Fastq& fq) const noexcept
  {
    return std::hash<std::string>{}(fq.seq);
  }
};

template <class Fastq>
struct Equal
{
  auto operator()(const Fastq& fq1, const Fastq& fq2) const noexcept
  {
    return fq1.seq == fq2.seq;
  }
};

template <class Fastq>
struct Less
{
  auto operator()(const Fastq& fq1, const Fastq& fq2) const noexcept
  {
    return fq1.seq < fq2.seq;
  }
};

template <class Fastq, int ContainerType>
class FastqReader {};

template <class Fastq>
class FastqReader <Fastq, ContainerType::MAP>
{
  public:
  std::map<Fastq, std::uint32_t, Less<Fastq>> buf;
  FastqReader(){};
  explicit FastqReader(std::istream& input)
  {
    for (auto&& i : input | pipeline::range::get_fastq(buf)) {} 
  }  
  /*
  explicit FastqReader (
      std::istream& input
    , std::map<Fastq, std::uint32_t, Less<Fastq>>& buf
  )
  {
    for (auto&& i : input | pipeline::range::get_fastq(buf)) {}
  }
  */
  auto& get_buf() noexcept
  {
    return buf;
  }
};


template <class Fastq>
class FastqReader <Fastq, ContainerType::UNORDERED_MAP>
{
  public:
  std::unordered_map<Fastq, std::uint32_t, Hash<Fastq>, Equal<Fastq>> buf;
 
  FastqReader(){}; 
  /*explicit FastqReader(const std::string& filename)
  {
    std::ifstream input(filename);
    for (auto&& i : input | pipeline::range::get_fastq(buf)) {}
  }*/
  explicit FastqReader(std::istream& input)
  {
    for (auto&& i : input | pipeline::range::get_fastq(buf)) {} 
  }  
  /*
  explicit FastqReader (
      std::istream& input
    , std::unordered_map<Fastq, std::uint32_t, Hash<Fastq>, Equal<Fastq>>& buf
  )
  {
    for (auto&& i : input | pipeline::range::get_fastq(buf)) {}
  }
  */
  auto& get_buf() noexcept
  {
    return buf;
  }
};

template <class Fastq>
class FastqReader <Fastq, ContainerType::VEC>
{
  public:
  std::vector<Fastq> buf;
 
  FastqReader(){}; 
  /*explicit FastqReader(const std::string& filename)
  {
    std::ifstream input(filename);
    for (auto&& i : input | pipeline::range::get_fastq(buf)) {}
  }*/
  explicit FastqReader(std::istream& input)
  {
    buf.reserve(500000);
    for (auto&& i : input | pipeline::range::get_fastq(buf)) {}
  }

  auto& get_buf() noexcept
  {
    return buf;
  }
};
}
