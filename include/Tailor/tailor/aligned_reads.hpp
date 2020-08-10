#pragma once
#include <iostream>
#include <vector>

namespace tailor {

template <class Fastq, class IndexType>
struct AlignedReads
{
  using fastq = Fastq;
  using seq_type = typename fastq::SEQ;
  using index_type = IndexType;
  
  bool has_data_;
  std::uint32_t reads_count_;
  const Fastq& fq_;
  seq_type tail_seq_;
  std::uint16_t flag_;
  int tail_pos_;
  std::pair<IndexType, IndexType> range_;
  
  // 1-to-1 mapping //
  std::vector<IndexType> real_pos_;
  std::vector<std::string> chrs_;

  IndexType NH_num_;
  IndexType prefix_match_len_;

  bool has_seed_mismatch_;
  bool has_tail_mismatch_;
  
  std::pair<int, char> seed_mismatch_;
  std::pair<int, char> tail_mismatch_;
  std::string MM_tag_;
  
  AlignedReads (
    bool has_data,
    const Fastq& fq,
    std::uint16_t flag,
    int tail_pos,
    std::pair<IndexType, IndexType> range,
    const std::vector<IndexType>& real_pos,
    const std::vector<std::string>& chrs,
    IndexType NH_num,
    IndexType prefix_match_len,
    bool has_seed_mismatch,
    bool has_tail_mismatch,
    std::pair<int, char> seed_mismatch,
    std::pair<int, char> tail_mismatch
  )
    : has_data_(has_data)
    , reads_count_(0)
    , fq_(fq)
    , flag_(flag)
    , tail_pos_(tail_pos)
    , range_(range)
    , real_pos_(real_pos)
    , chrs_(chrs)
    , NH_num_(NH_num)
    , prefix_match_len_(prefix_match_len)
    , has_seed_mismatch_(has_seed_mismatch)
    , has_tail_mismatch_(has_tail_mismatch)
    , seed_mismatch_(seed_mismatch)
    , tail_mismatch_(tail_mismatch)
    , MM_tag_("")
  {}

  AlignedReads (Fastq&& fq): fq_(fq), flag_(-1), has_data_(false) {}

  void set_reads_count(std::uint32_t count)
  {
    reads_count_ = count;
  } 

  void print(std::ostream& out = std::cout) const 
  {
    if (flag_ != 0 && flag_ != 16) return;
    out << std::endl;
    out << "******************\n";
    out << fq_ << std::endl;
    
    out << "flag: " << flag_ << "\n";
    out << "real_pos: [ ";
    for (auto& i : real_pos_)
      out << i << " ";
    out << "]\n";
    out << "chrs: [ ";
    for (auto& i : chrs_)
      out << i << " ";
    out << "]\n";
    
    out << "prefix_len: " << prefix_match_len_ << std::endl;
    out << "tail_pos: " << tail_pos_ << std::endl;
    
    if (has_seed_mismatch_)
      out << "seed_mismatch: " << seed_mismatch_.first << ", " << seed_mismatch_.second << std::endl;
    if (has_tail_mismatch_)
      out << "tail_mismatch: " << tail_mismatch_.first << ", " << tail_mismatch_.second << std::endl;
    out << "******************\n";   
  }
};
}
