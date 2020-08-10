
#pragma once
#include <range/v3/all.hpp>
#include <algorithm>
#include <atomic>
#include <mutex>

namespace pipeline::range {

std::mutex output_mutex;
//std::mutex adding_mutex;

enum SAM_INDEX
{
    QNAME
  , FLAG
  , RNAME
  , POS
  , MAPQ
  , CIGAR
  , RNEXT
  , PNEXT
  , TLEN
  , SEQ
  , QUAL
  , TAGS
};

template <class IndexType, class SEQ>
using SamType = std::tuple<
    std::string // QNAME
  , std::uint16_t // FLAG
  , std::string // RNAME
  , IndexType // POS
  , std::uint32_t // MAPQ
  , std::string // CIGAR
  , std::string // RNEXT
  , std::uint32_t // PNEXT
  , std::int32_t // TLEN
  , SEQ // SEQ
  , std::string // QUAL
  , std::string // TAGS
>;

template <class SamType, std::size_t N>
void dump (SamType&& sam, std::ostream& out)
{
  if constexpr (N == TAGS)
  {
    out << std::get<N>(sam) << "\n";
  }
  else
  {
    out << std::get<N>(sam) << "\t";
    dump<SamType, N+1>(std::move(sam), out);
  }
}

auto ToSAM = [](auto&& AlignedReads, std::ostream& out) {
   
    using SEQ = typename std::remove_reference_t<decltype(AlignedReads)>::seq_type;
    using IndexType = typename std::remove_reference_t<decltype(AlignedReads)>::index_type;
    
    std::uint16_t flag = AlignedReads.flag_;
    auto tail_pos = AlignedReads.tail_pos_;
    
    auto reads_count = AlignedReads.reads_count_;
    const auto& fq = AlignedReads.fq_;
    auto& real_pos = AlignedReads.real_pos_;
    auto& chrs = AlignedReads.chrs_;
    auto& MM_tag = AlignedReads.MM_tag_;
    auto multi_size = real_pos.size();
    auto NH_num = std::to_string(AlignedReads.NH_num_);
    auto prefix_match_len = AlignedReads.prefix_match_len_; // same as origin_tail_pos
    auto seq_size = static_cast<int>(prefix_match_len) + tail_pos;

    auto has_seed_mismatch = AlignedReads.has_seed_mismatch_;
    auto has_tail_mismatch = AlignedReads.has_tail_mismatch_;
    
    auto seed_mm_pos = AlignedReads.seed_mismatch_.first;
    auto seed_mm_char = AlignedReads.seed_mismatch_.second;
    
    auto tail_mm_pos = AlignedReads.tail_mismatch_.first;
    auto tail_mm_char = AlignedReads.tail_mismatch_.second;

    auto tag = "NH:i:" + NH_num;
    auto MAPQ = 255 - tail_pos - 1;

    if (flag == 0)
		{
      auto CIGAR = 
          tail_pos == -1
        ? std::to_string(prefix_match_len) + 'M'
        : std::to_string(prefix_match_len) + 'M' 
        + std::to_string(tail_pos + 1) + 'S';
      
      if (has_seed_mismatch || has_tail_mismatch)
      {
        tag += "\tMD:Z:";
        if (has_seed_mismatch && !has_tail_mismatch)
        {
          switch (seed_mm_char)
          {
            case 'A': seed_mm_char = 'T'; break;
            case 'T': seed_mm_char = 'A'; break;
            case 'C': seed_mm_char = 'G'; break;
            case 'G': seed_mm_char = 'C'; break;
          }
          seed_mm_pos = seq_size - seed_mm_pos;
          auto back_len = seq_size - tail_pos - seed_mm_pos - 1;
          tag += (seed_mm_pos ? std::to_string(seed_mm_pos) : "") + seed_mm_char + std::to_string(back_len);
          MM_tag = (seed_mm_pos ? std::to_string(seed_mm_pos) : "") + seed_mm_char + std::to_string(back_len);
        }
        else if (!has_seed_mismatch && has_tail_mismatch)
        {
          switch (tail_mm_char)
          {
            case 'A': tail_mm_char = 'T'; break;
            case 'T': tail_mm_char = 'A'; break;
            case 'C': tail_mm_char = 'G'; break;
            case 'G': tail_mm_char = 'C'; break;
          }
          tail_mm_pos = seq_size - tail_mm_pos;
          auto back_len = seq_size - tail_pos - tail_mm_pos - 1;
          tag += std::to_string(tail_mm_pos) + tail_mm_char + std::to_string(back_len);
          MM_tag = std::to_string(tail_mm_pos) + tail_mm_char + std::to_string(back_len);
        }
        else if (has_seed_mismatch && has_tail_mismatch)
        {
          switch (seed_mm_char)
          {
            case 'A': seed_mm_char = 'T'; break;
            case 'T': seed_mm_char = 'A'; break;
            case 'C': seed_mm_char = 'G'; break;
            case 'G': seed_mm_char = 'C'; break;
          } 
          switch (tail_mm_char)
          {
            case 'A': tail_mm_char = 'T'; break;
            case 'T': tail_mm_char = 'A'; break;
            case 'C': tail_mm_char = 'G'; break;
            case 'G': tail_mm_char = 'C'; break;
          }
          seed_mm_pos = seq_size - seed_mm_pos;
          tail_mm_pos = seq_size - tail_mm_pos;
          auto mid_len = tail_mm_pos - seed_mm_pos - 1;
          auto back_len = seq_size - tail_pos - tail_mm_pos - 1;
          tag += (seed_mm_pos ? std::to_string(seed_mm_pos) : "") + seed_mm_char + std::to_string(mid_len) + tail_mm_char + std::to_string(back_len);
          MM_tag = (seed_mm_pos ? std::to_string(seed_mm_pos) : "") + seed_mm_char + std::to_string(mid_len) + tail_mm_char + std::to_string(back_len);
        }
      }

      if (tail_pos != -1)
      {
        SEQ tail_seq;
        for (auto i = prefix_match_len; i <= seq_size; ++i)
        {
          tail_seq.push_back(fq.seq[i]);
        }
        AlignedReads.tail_seq_ = tail_seq; 
        tag += "\tTL:Z:" + tail_seq;
      }
      
      std::scoped_lock<std::mutex> sl(output_mutex);

      for (std::uint32_t n = 0; n < reads_count; ++n)
      {
        for (int i = 0; i < multi_size; ++i)
        {
          dump<SamType<IndexType, SEQ>, 0>
          ( SamType<IndexType, SEQ>{
              fq.name.substr(0, fq.name.find(" "))
            , flag
            , chrs[i].substr(0, chrs[i].find(" "))
            , real_pos[i]
            , MAPQ
            , CIGAR 
            , "*"
            , 0
            , 0
            , fq.seq
            , fq.seq_qual
            , tag
          }, out);
        }
      }
    }
    else if (flag == 16)
    {
      auto CIGAR = 
          tail_pos == -1
        ? std::to_string(prefix_match_len) + 'M'
        : std::to_string(tail_pos + 1) + 'S' 
        + std::to_string(prefix_match_len) + 'M';
      
      if (has_seed_mismatch || has_tail_mismatch)
      {
        tag += "\tMD:Z:";
        if (has_seed_mismatch && !has_tail_mismatch)
        {
          seed_mm_pos = seq_size - seed_mm_pos;
          auto back_len = seq_size - tail_pos - seed_mm_pos - 1;
          tag += std::to_string(back_len) + seed_mm_char + (seed_mm_pos ? std::to_string(seed_mm_pos) : "");
          MM_tag = std::to_string(back_len) + seed_mm_char + (seed_mm_pos ? std::to_string(seed_mm_pos) : "");
        }
        else if (!has_seed_mismatch && has_tail_mismatch)
        {
          tail_mm_pos = seq_size - tail_mm_pos;
          auto back_len = seq_size - tail_pos - tail_mm_pos - 1;
          tag += std::to_string(back_len) + tail_mm_char + std::to_string(tail_mm_pos);
          MM_tag = std::to_string(back_len) + tail_mm_char + std::to_string(tail_mm_pos);
        }
        else if (has_seed_mismatch && has_tail_mismatch)
        {
          seed_mm_pos = seq_size - seed_mm_pos;
          tail_mm_pos = seq_size - tail_mm_pos;
          auto mid_len = tail_mm_pos - seed_mm_pos - 1;
          auto back_len = seq_size - tail_pos - tail_mm_pos - 1;
          tag += std::to_string(back_len) + tail_mm_char + std::to_string(mid_len) + seed_mm_char + (seed_mm_pos ? std::to_string(seed_mm_pos) : "");
          MM_tag = std::to_string(back_len) + tail_mm_char + std::to_string(mid_len) + seed_mm_char + (seed_mm_pos ? std::to_string(seed_mm_pos) : "");
        }
      }
      
      if (tail_pos != -1)
      {
        SEQ tail_seq;
        for (auto i = prefix_match_len; i <= seq_size; ++i)
        {
          tail_seq.push_back(fq.seq[i]);
        }
        AlignedReads.tail_seq_ = tail_seq;
        tag += "\tTL:Z:" + tail_seq;
      }
     
      auto seq_qual = fq.seq_qual; 
      std::reverse(seq_qual.begin(), seq_qual.end());
      std::scoped_lock<std::mutex> sl(output_mutex);

      for (std::uint32_t n = 0; n < reads_count; ++n)
      {
        for (int i = 0; i < multi_size; ++i)
        {  
          dump<SamType<IndexType, SEQ>, 0>
          ( SamType<IndexType, SEQ>{
              fq.name.substr(0, fq.name.find(" "))
            , flag
            , chrs[i].substr(0, chrs[i].find(" "))
            , real_pos[i]
            , MAPQ
            , CIGAR 
            , "*"
            , 0
            , 0
            , fq.get_antisense()
            , seq_qual
            , tag
          }, out);
        }
      }
    }
};

auto output_sam(std::atomic_uint32_t& mappable_count, std::ostream& out = std::cout)
{
  return ranges::view::transform (
    [&mappable_count, &out](auto&& AlignedReads) {
      ToSAM(std::move(AlignedReads), out);
      
      //std::scoped_lock<std::mutex> sl(adding_mutex);  
      if (AlignedReads.has_data_)
        mappable_count.fetch_add(AlignedReads.reads_count_, std::memory_order_relaxed);
        //mappable_count += AlignedReads.reads_count_;
      
      // have to test if AlignedReads still have data after move to ToSAM function
      return AlignedReads;
  });
}

auto output_sam_set(std::atomic_uint32_t& mappable_count, std::ostream& out = std::cout)
{
  return ranges::view::transform (
    [&mappable_count, &out](auto&& AlignedReads_v) {
      
      if (AlignedReads_v.empty())
        return AlignedReads_v;
      else {
        mappable_count.fetch_add(AlignedReads_v.front().reads_count_, std::memory_order_relaxed);
        
        //if (AlignedReads_v.front().tail_pos_ != AlignedReads_v.back().tail_pos_) {
          auto min_tail_pos = std::min_element (
              AlignedReads_v.cbegin()
            , AlignedReads_v.cend()
            , [](auto&& a, auto&& b) { return a.tail_pos_ < b.tail_pos_; } 
          );
          //auto index_of_min = std::distance(AlignedReads_v.cbegin(), min_tail_pos);
          std::remove_reference_t<decltype(AlignedReads_v)> results;
          
          for (auto&& reads : AlignedReads_v)
          {
            if (reads.tail_pos_ <= (*min_tail_pos).tail_pos_)
              results.emplace_back(reads);
          }
          for (auto&& i : results)
          {
            ToSAM(std::move(i), out);
          }
          return results;
        //}
        /*else {
          for (auto&& i : AlignedReads_v)
          {
            ToSAM(std::move(i), out);
          }
          return AlignedReads_v;
        }*/
      }
  });
}



}
