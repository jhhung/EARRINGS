#pragma once

#include <string>
#include <Tailor/tailor/paras.hpp>
#include <Tailor/tailor/aligned_reads.hpp>
#include <boost/filesystem.hpp>
#include <Biovoltron/format/fastq.hpp>
#include <Biovoltron/format/newlineSeperatedFasta.hpp>
#include <Biovoltron/indexer/FMIndex.hpp>
#include <Biovoltron/string_sorter/RadixSort.hpp>
#include <Biovoltron/string_sorter/SBWT.hpp>


namespace tailor {

constexpr auto Interval = 9;
constexpr auto CharTypeNum = 4;
constexpr auto LookupStrLen = 12;
constexpr auto PrefixLen = 256;
constexpr auto IsSBWT = true;
constexpr auto ASCIISize = 256;

using RadixSort = biovoltron::string_sorter::RadixSort<>;
using IntType = std::uint32_t;
using SuffixArrayType = std::vector<IntType>;
using SBWT = biovoltron::string_sorter::SBWT<SuffixArrayType, RadixSort>;
using Seq2bits = biovoltron::vector<biovoltron::char_type>; 
using Seq8bits = std::string; 

using Seq = Seq8bits;
using Fastq = biovoltron::format::FASTQ<Seq>;
using Fasta = biovoltron::format::FASTA<Seq>;

using Indexer = biovoltron::indexer::FMIndex<
  Seq, SuffixArrayType, SBWT, Interval, 
  CharTypeNum, LookupStrLen, PrefixLen, IsSBWT>; 

using FastqReads = AlignedReads<Fastq, typename Indexer::IndexType>;
using FastaReads = AlignedReads<Fasta, typename Indexer::IndexType>;
  
unsigned long fast_rand()
{
  unsigned long t;
  fast_rand_seed[0] ^= fast_rand_seed[0] << 16;
  fast_rand_seed[0] ^= fast_rand_seed[0] >> 5;
  fast_rand_seed[0] ^= fast_rand_seed[0] << 1;
  
  t = fast_rand_seed[0];
  fast_rand_seed[0] = fast_rand_seed[1];
  fast_rand_seed[1] = fast_rand_seed[2];
  fast_rand_seed[2] = t ^ fast_rand_seed[0] ^ fast_rand_seed[1];
  
  return fast_rand_seed[2];
}

constexpr std::array<std::uint32_t, ASCIISize> fill_array()
{
	std::array<std::uint32_t, ASCIISize> res{0};
  for (int i = 0; i < ASCIISize; ++i)
  {
    res[i] = 0;
  }
  return res;
}

template <typename OArray>
constexpr auto char_to_order_init(OArray&& o2c)
{
  std::array<std::uint32_t, ASCIISize> c2o = fill_array();
  auto count = 0;
  for (std::size_t i = 0; i < o2c.size(); ++i)
  {
    c2o[ o2c.at(i) ] = count;
    ++count;
  }
  return c2o;
}

constexpr std::array<std::uint32_t, CharTypeNum> order_to_char{'A', 'C', 'G', 'T'};
constexpr std::array<std::uint32_t, ASCIISize> char_to_order = char_to_order_init(order_to_char);
  
template <typename SEQ>
void reverse_c(SEQ& seq) {}

// the base in the middle wouldn't be filpped when seq size is odd.
template <>
void reverse_c<Seq8bits>(Seq8bits& seq)
{
  auto first = seq.begin();
  auto last = seq.end(); 
  while ((first != last) && (first != --last)) 
  {
        *first = order_to_char[~(char_to_order[*first]) & 3]; 
        *last = order_to_char[~(char_to_order[*last]) & 3];
        std::iter_swap(first++, last);
  }
  if (seq.size() & 1)
  {
    auto mid_idx = seq.size() >> 1;
    seq[mid_idx] = order_to_char[~(char_to_order[seq[mid_idx]]) & 3];
  }
}

template <>
void reverse_c<Seq2bits>(Seq2bits& seq)
{    
  seq = Seq2bits{seq.rbegin(), seq.rend()};
  seq.flip();
}

template <bool boolType>
class TailorMain {};

template <>
class TailorMain <true>
{
  TailParas pt;
  TailorSearcher<Fasta, Indexer, TailParas, FastaReads> table;
  
  public:
  
  explicit TailorMain(
      const std::size_t n_thread
    , const uint32_t min_prefix_len
    , const uint32_t min_multi
    , const std::string& prefix_name
    , bool allow_mm = true 
    )
    : pt(n_thread, min_prefix_len, min_multi, allow_mm)
    , table(char_to_order, order_to_char, prefix_name, pt)
  {}
  auto& get_table() const noexcept
  {
    return table;
  }
  auto& get_paras() const noexcept
  {
    return pt;
  }
};

template <>
class TailorMain <false>
{
  TailParas pt;
  TailorSearcher<Fastq, Indexer, TailParas, FastqReads> table;
  
  public:
  
  explicit TailorMain(
      const std::size_t n_thread
    , const uint32_t min_prefix_len
    , const uint32_t min_multi
    , const std::string& prefix_name
    , bool allow_mm = true 
    )
    : pt(n_thread, min_prefix_len, min_multi, allow_mm)
    , table(char_to_order, order_to_char, prefix_name, pt)
  {}
  auto& get_table() const noexcept
  {
    return table;
  }
  auto& get_paras() const noexcept
  {
    return pt;
  }
};

class build_tables 
{
  private:
  using IndexType = typename Indexer::IndexType;

  template <
    typename IStream, 
    typename SEQ,
    typename Vec, 
    typename VecP, 
    typename Arr
  >
  void read_file (
    IStream&& is, 
    SEQ&& seq,
    Vec&& chr_v, 
    VecP&& seg_info, 
    VecP&& n_table, 
    Arr&& order_to_char_
  ) const
  {
      std::string buf;
      std::string seqs;
      std::string chr_name;
      auto seg_pos(0);
  
      while (std::getline(is, chr_name))
      //for (int i = 0; std::getline(is, tmp2); ++i)
      {
        if (chr_name.front() != '>')
          throw std::runtime_error(
            "ERROR: Input file format error\n"
          );
        chr_v.emplace_back(chr_name.substr(1));
        while (is.peek() != '>' && std::getline(is, seqs))
        {
          buf += seqs;
        }
        //if (i == 0) continue;
        seg_info.emplace_back(seg_pos, buf.size());
        seg_pos += buf.size();
        //std::cout << "before: " << seq.capacity() << " " << buf.size() << std::endl;
        seq.reserve(seq.size() + buf.size());
        //std::cout << "after: " << seq.capacity() << " " << buf.size() << std::endl;
        for (auto it(buf.cbegin()); it != buf.cend(); it++)
        {
          if (*it == 'N' || *it == 'n')
          {
            n_table.emplace_back(
              std::distance(buf.cbegin(), it) + seg_info.back().first, 1);
            seq.push_back(order_to_char_[fast_rand() % CharTypeNum]);
  
            auto& n_count(n_table.back().second);
            for (it++; it != buf.cend(); it++)
            {
  						if (*it == 'N' || *it == 'n')
  						{
  							n_count++;
  							seq.push_back(
  								order_to_char_[fast_rand() % CharTypeNum]
  							);
  						}
  						else
  							break;
            }
          }
          if (it == buf.cend()) 
              break;
          if (*it > 'T')
            seq.push_back(*it - ' ');
          else
            seq.push_back(*it);
        }
        buf.clear();
      }
  }
  
  template <typename Vec>
  void reverse_seg_n_table (
    Vec&& seg_info,
    Vec&& n_table,
    Vec&& rc_seg_info,
    Vec&& rc_n_table
  ) const
  {
    std::uint32_t total = 0;
    auto r_seg_iter = seg_info.crbegin();
    auto r_n_iter = n_table.crbegin();
  
    while (r_seg_iter != seg_info.crend())
    {
      rc_seg_info.emplace_back(total, (*r_seg_iter).second);
      total += (*r_seg_iter).second;
      ++r_seg_iter;
    }
    
    while (r_n_iter != n_table.crend())
    {
      rc_n_table.emplace_back(total - (*r_n_iter).first - (*r_n_iter).second, (*r_n_iter).second);
      ++r_n_iter;
    }
  }
  
  public:
  void operator()(const std::string& filename, const std::string& prefix_name) const
  {
    std::ifstream input(filename);
    if (!input.is_open())
      throw std::runtime_error("Can't open input reference file\n");
    Seq seq;
  
    Indexer index_table(char_to_order, order_to_char);
    Indexer rc_index_table(char_to_order, order_to_char);

    read_file(input, seq, index_table.chr_names, index_table.seg_info, index_table.n_table, order_to_char);
    reverse_seg_n_table(index_table.seg_info, index_table.n_table, rc_index_table.seg_info, rc_index_table.n_table);

    SBWT sbwt(90*1024*1024, 20*1024*1024);
    SBWT rc_sbwt(90*1024*1024, 20*1024*1024);
    
    index_table.build(seq, sbwt);
    
    reverse_c(seq);
    
    rc_index_table.build(seq, rc_sbwt);
    
    index_table.save(prefix_name + ".table");
    rc_index_table.save(prefix_name + ".rc_table");
  }
};

class build_tables_long 
{
  private:
  using IndexType = typename Indexer::IndexType;

  template <
    typename IStream, 
    typename SEQ,
    typename Vec, 
    typename VecP, 
    typename Arr
  >
  void read_file (
    IStream&& is, 
    SEQ&& seq,
    Vec&& chr_v, 
    VecP&& seg_info, 
    VecP&& n_table, 
    Arr&& order_to_char_
  ) const
  {
      std::string buf;
      std::string seqs;
      std::string chr_name;
      auto seg_pos(0);
  
      while (std::getline(is, chr_name))
      //for (int i = 0; std::getline(is, tmp2); ++i)
      {
        if (chr_name.front() != '>')
          throw std::runtime_error(
            "ERROR: Input file format error\n"
          );
        chr_v.emplace_back(chr_name.substr(1));
        while (is.peek() != '>' && std::getline(is, seqs))
        {
          buf += seqs;
        }
        //if (i == 0) continue;
        seg_info.emplace_back(seg_pos, buf.size());
        seg_pos += buf.size();
        //std::cout << "before: " << seq.capacity() << " " << buf.size() << std::endl;
        seq.reserve(seq.size() + buf.size());
        //std::cout << "after: " << seq.capacity() << " " << buf.size() << std::endl;
        for (auto it(buf.cbegin()); it != buf.cend(); it++)
        {
          if (*it == 'N' || *it == 'n')
          {
            n_table.emplace_back(
              std::distance(buf.cbegin(), it) + seg_info.back().first, 1);
            seq.push_back(order_to_char_[fast_rand() % CharTypeNum]);
  
            auto& n_count(n_table.back().second);
            for (it++; it != buf.cend(); it++)
            {
  						if (*it == 'N' || *it == 'n')
  						{
  							n_count++;
  							seq.push_back(
  								order_to_char_[fast_rand() % CharTypeNum]
  							);
  						}
  						else
  							break;
            }
          }
          if (it == buf.cend()) 
              break;
          if (*it > 'T')
            seq.push_back(*it - ' ');
          else
            seq.push_back(*it);
        }
        buf.clear();
      }
  }
  
  template <typename Vec>
  void reverse_seg_n_table (
    Vec&& seg_info,
    Vec&& n_table,
    Vec&& rc_seg_info,
    Vec&& rc_n_table
  ) const
  {
    std::uint32_t total = 0;
    auto r_seg_iter = seg_info.crbegin();
    auto r_n_iter = n_table.crbegin();
  
    while (r_seg_iter != seg_info.crend())
    {
      rc_seg_info.emplace_back(total, (*r_seg_iter).second);
      total += (*r_seg_iter).second;
      ++r_seg_iter;
    }
    
    while (r_n_iter != n_table.crend())
    {
      rc_n_table.emplace_back(total - (*r_n_iter).first - (*r_n_iter).second, (*r_n_iter).second);
      ++r_n_iter;
    }
  }
  
  public:
  void operator()(const std::string& filename, const std::string& prefix_name) const
  {
    std::ifstream input(filename);
    Seq seq;
  
    Indexer index_table(char_to_order, order_to_char);
    Indexer rc_index_table(char_to_order, order_to_char);

    read_file(input, seq, index_table.chr_names, index_table.seg_info, index_table.n_table, order_to_char);
    reverse_seg_n_table(index_table.seg_info, index_table.n_table, rc_index_table.seg_info, rc_index_table.n_table);

    SBWT sbwt(90*1024*1024, 20*1024*1024);
    SBWT rc_sbwt(90*1024*1024, 20*1024*1024);
    auto tmp = seq; 
    
    index_table.build(seq, sbwt);
    index_table.save(prefix_name + "table");
    
    std::cout << "starting reverse..." << std::endl;
    reverse_c(seq);
    std::cout << "reverse done!" << std::endl;
    std::cout << "starting cancat..." << std::endl;
    tmp += seq;
    std::cout << "cancat done" << std::endl;
    
    std::cout << "starting build reverse table..." << std::endl;
    rc_index_table.build(tmp, rc_sbwt);
    std::cout << "build reverse table done" << std::endl;
    
    rc_index_table.save(prefix_name + "rc_table");
  }
};
}

bool checkIndexIntact (const std::string& prefixname) 
{
  if (boost::filesystem::exists (prefixname + "table") &&
      boost::filesystem::exists (prefixname + "rc_table"))
  {
    return true;
  } 
  else 
  {
    return false;
  }
}

