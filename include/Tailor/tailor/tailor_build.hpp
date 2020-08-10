/*
# Tailor, a BWT-based aligner for non-templated RNA tailing
# Copyright (C) 2014 Min-Te Chou, Bo W Han, Jui-Hung Hung
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
#pragma once
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <Biovoltron/base_vector.hpp>
#include <Biovoltron/indexer/FMIndex.hpp>
//#include <Biovoltron/string_sorter/BucketSort.hpp>
#include <Biovoltron/string_sorter/RadixSort.hpp>
#include <Biovoltron/string_sorter/SBWT.hpp>
//#include <Tailor/tailor/tailor_index.hpp>
#include <Tailor/tailor/tailor_search.hpp>
#include <chrono>

static unsigned long fast_rand()
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

template <typename CArray, typename OArray>
void char_to_order_init(CArray&& c2o, OArray&& o2c)
{
  c2o.fill(0);
  auto count = 0;
  for (std::size_t i = 0; i < o2c.size(); ++i)
  {
    c2o[ o2c.at(i) ] = count;
    ++count;
  }
}

namespace tailor{

// add inline to avoid name collision
// use "uint32_t" to declare is safer than auto, 
// as auto is set to "int" by compiler and pass "uint32_t" as template para 
// so, the process involve in type conversion.
constexpr auto Interval = 9;
constexpr auto CharTypeNum = 4;
constexpr auto LookupStrLen = 12;
constexpr auto PrefixLen = 256;
constexpr auto IsSBWT = true;
constexpr auto ASCIISize = 256;

using RadixSort = biovoltron::string_sorter::RadixSort<>;
//using BucketSort = biovoltron::string_sorter::BucketSort<>;
using IntType = std::uint32_t;
using SuffixArrayType = std::vector<IntType>;
using SBWT = biovoltron::string_sorter::SBWT<SuffixArrayType, RadixSort>;

using Seq2bits = biovoltron::vector<biovoltron::char_type>; 
using Seq8bits = std::string; 

using Seq = Seq8bits;

/*
  using Indexer = tailor::TailorIndexer<
  Seq, SuffixArrayType, SBWT, Interval, 
  CharTypeNum, LookupStrLen, ASCIISize>;
*/

using Indexer = biovoltron::indexer::FMIndex<
  Seq, SuffixArrayType, SBWT, Interval, 
  CharTypeNum, LookupStrLen, PrefixLen, IsSBWT>;

std::array<std::uint32_t, ASCIISize> char_to_order;
std::array<std::uint32_t, CharTypeNum> order_to_char{'A', 'C', 'G', 'T'};

//template <typename CArray, typename OArray>
//char_to_order_init(CArray &&, OArray &&) -> char_to_order_init<CArray, OArray>;

//char_to_order_init(char_to_order, order_to_char);

//using IndexType = typename SuffixArrayType::value_type;

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
}

template <>
void reverse_c<Seq2bits>(Seq2bits& seq)
{
  seq = Seq2bits{seq.rbegin(), seq.rend()};
  seq.flip();
}

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
)
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
          seq.push_back(order_to_char_[::fast_rand() % CharTypeNum]);

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
        seq.push_back(*it);
      }
      buf.clear();
    }
    //buf.clear();
    //buf.shrink_to_fit();
}

template <typename Vec>
void reverse_seg_n_table (
  Vec&& seg_info,
  Vec&& n_table,
  Vec&& rc_seg_info,
  Vec&& rc_n_table
)
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


void build_tables (const std::string& filename, const std::string& prefix_name)
{

  using IndexType = typename Indexer::IndexType;

  ::char_to_order_init(char_to_order, order_to_char);
  //std::vector<std::pair<IndexType, IndexType>> seg_info, n_table;
  
  std::ifstream input(filename);
  Seq seq;

  Indexer index_table(char_to_order, order_to_char);
  Indexer rc_index_table(char_to_order, order_to_char);

  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  
  read_file(input, seq, index_table.chr_names, index_table.seg_info, index_table.n_table, order_to_char);
  std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  std::cout << "reading time: " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " s" << std::endl;
  
  std::chrono::steady_clock::time_point n_table1 = std::chrono::steady_clock::now();
  reverse_seg_n_table(index_table.seg_info, index_table.n_table, rc_index_table.seg_info, rc_index_table.n_table);
  std::chrono::steady_clock::time_point n_table2 = std::chrono::steady_clock::now();
  std::cout << "reverse table time: " << std::chrono::duration_cast<std::chrono::seconds>(n_table2 - n_table1).count() << " s" << std::endl;


  // n_table testing //

  /*std::cout << "seg_info: " << "\n";
  for (auto& i : index_table.seg_info)
  {
    std::cout << "pos: " << i.first << "\n";
    std::cout << "size: " << i.second << "\n";
  }
  std::cout << "n_table: " << "\n";
  for (auto& i : index_table.n_table)
  {
    std::cout << "pos: " << i.first << "\n";
    std::cout << "size: " << i.second << "\n"; 
  }
  std::cout << "=========== reverse ==========\n";
  std::cout << "seg_info: " << "\n";
  for (auto& i : rc_index_table.seg_info)
  {
    std::cout << "pos: " << i.first << "\n";
    std::cout << "size: " << i.second << "\n";
  }
  std::cout << "n_table: " << "\n";
  for (auto& i : rc_index_table.n_table)
  {
    std::cout << "pos: " << i.first << "\n";
    std::cout << "size: " << i.second << "\n"; 
  }*/
  
  //rc_index_table.seg_info = index_table.seg_info;
  //rc_index_table.n_table = index_table.n_table;
  

  SBWT sbwt(90*1024*1024, 20*1024*1024);
  SBWT rc_sbwt(90*1024*1024, 20*1024*1024);
  
  std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
  //Indexer index_table(prefixname, seq, char_to_order, order_to_char, seg_info, n_table, sbwt, 0);
  index_table.build(seq, sbwt);
  std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
  std::cout << "building fm table time: " << std::chrono::duration_cast<std::chrono::seconds>(t4 - t3).count() << " s" << std::endl;
  
  //auto antiseq = seq;
  std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();
  //reverse_c(antiseq);
  reverse_c(seq);
  std::chrono::steady_clock::time_point t6 = std::chrono::steady_clock::now();
  std::cout << "reverse time: " << std::chrono::duration_cast<std::chrono::seconds>(t6 - t5).count() << " s" << std::endl;
 
  // TO TEST CORRECTNESS WITH TAILOR // 
  /*auto antiseq = seq;
  reverse_c(antiseq);
  seq += antiseq;*/

  std::chrono::steady_clock::time_point t7 = std::chrono::steady_clock::now();
  //Indexer rc_index_table(prefixname+"rc_", seq, char_to_order, order_to_char, seg_info, n_table, rc_sbwt, 0);
  rc_index_table.build(seq, rc_sbwt);
  std::chrono::steady_clock::time_point t8 = std::chrono::steady_clock::now();
  std::cout << "building rc fm table time: " << std::chrono::duration_cast<std::chrono::seconds>(t8 - t7).count() << " s" << std::endl;
  

  index_table.save(prefix_name + "table");
  rc_index_table.save(prefix_name + "rc_table");

  //seq.clear();
  //seq.shrink_to_fit(); 
  
  //std::cout << "fm_index: " << index_table.get_bwt_seq() << std::endl;
  //std::cout << "rc_fm_index: " << rc_index_table.get_bwt_seq() << std::endl;
}

template <typename OStream, typename TailParas>
void tailor_mapping (const std::string& prefix_name, const std::string& input_name, OStream&& output, TailParas&& para_pack)
{
  ::char_to_order_init(char_to_order, order_to_char);
  output << "@HD" << '\t' << "VN:1.0" << '\t' << "SO:unsorted\n";
  Indexer fm_index(char_to_order, order_to_char);
  Indexer rc_fm_index(char_to_order, order_to_char);
  fm_index.load(prefix_name + "table");
  rc_fm_index.load(prefix_name + "rc_table");
 
  auto chr = fm_index.chr_names.begin();
  for (const auto& seg : fm_index.seg_info)
  {
    output << "@SQ\tSN:" << *chr << "\tLN:" << seg.second << '\n';
    ++chr;
  }
  std::ifstream input(input_name);
  // TODO: n threads for searcher
  //TailorSearcher searcher(std::move(fm_index), std::move(rc_fm_index), std::move(input), std::move(output), para_pack);
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  searcher::TailorSearcher searcher(fm_index, rc_fm_index, std::move(input), std::move(output), para_pack);
  std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  std::cout << "search time: " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " s" << std::endl;

  //output << "fm_index: " << fm_index.get_bwt_seq() << std::endl;
  //output << "rc_fm_index: " << rc_fm_index.get_bwt_seq() << std::endl;
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

}

