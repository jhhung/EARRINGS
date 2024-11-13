#pragma once

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <array>
#include <tuple>

namespace tailor {
template <
  class Fastq,
  class Indexer,
  class TailParas,
  class AlignedReads
>
class TailorSearcher
{

public:
  using IndexType = typename Indexer::IndexType;
  using SEQ = typename Indexer::SeqType;
  using FMIdxRange = std::pair<IndexType, IndexType>;
  using Mismatchtype = std::pair<int, char>;
  using SeedMMType = std::tuple<IndexType, int, char>;
  using FASTQ = Fastq;

  static constexpr std::array<std::uint8_t, 4> ch_set {'A', 'C', 'G', 'T'};
  static constexpr bool IsRC = true;

private:
  Indexer fm_index;
  Indexer rc_fm_index;
  const TailParas& para_pack;
  
public:
  explicit TailorSearcher(
    const std::array<IndexType, 256>& c_to_o,
    const std::array<IndexType, 4>& o_to_c,
    const std::string& prefix_name,
    const TailParas& pack
  )
    : fm_index(c_to_o, o_to_c)
    , rc_fm_index(c_to_o, o_to_c)
    , para_pack(pack)
  {
    fm_index.load(prefix_name + ".table");
    rc_fm_index.load(prefix_name + ".rc_table");
    // Test //
    std::cout << "load table success" << std::endl;
    pack.print();   
  }

  std::vector<AlignedReads> align(const Fastq& fq, const std::uint32_t& reads_count = 1) const
  {
    std::vector<AlignedReads> res;
    if (fq.seq.size() < para_pack.min_prefix_len || !fq.n_base_info_table.empty())
      return res;

    auto rc_query = fq.get_antisense();

    std::vector<FMIdxRange> sense_v;
    std::vector<FMIdxRange> antisense_v;
    
    auto fm_mm_idx = fm_index.sbwt_exact_match_by_base(fq.seq, 0, sense_v);
    auto rc_fm_mm_idx = rc_fm_index.sbwt_exact_match_by_base(fq.seq, 0, antisense_v);

    auto fm_multi_align = fm_mm_idx.second.second - fm_mm_idx.second.first;
    auto rc_fm_multi_align = rc_fm_mm_idx.second.second - rc_fm_mm_idx.second.first;
    
    if ( fm_multi_align + rc_fm_multi_align > para_pack.min_multi )
    {
      return res;
    }
    
    search(fq, sense_v, rc_query, fm_mm_idx, res);
    search<IsRC>(fq, antisense_v, rc_query, rc_fm_mm_idx, res);
    for (auto&& r : res)
    {
      r.set_reads_count(reads_count);
    }
    if (res.size() > 1)
    {
      if (res.front().real_pos_.size() + res.back().real_pos_.size() > para_pack.min_multi)
      {
        res.clear();
      }
    }
    return res;
  }

private:
  template <bool IsRC = false>
  void search (
    const Fastq& fq,
    const std::vector<FMIdxRange>& match_pos_record,
    const SEQ& query,
    const std::pair<IndexType, FMIdxRange>& mm_idx,
    std::vector<AlignedReads>& res
  ) const
  {
    IndexType prefix_match_len = mm_idx.first;
    if (prefix_match_len == -2)
        return;
    
    if (prefix_match_len == -1)
    {
      record_res<IsRC>(fq, match_pos_record.back(), -1, res);
    }
    else if (prefix_match_len < para_pack.min_prefix_len && para_pack.allow_mm)
    {
      std::vector<SeedMMType> seed_mismatch_record;
      std::vector<FMIdxRange> last_match_pos_record;
      IndexType mismatch_idx = query.size() - prefix_match_len - 1;
      std::vector<int> tail_pos_record;
      
      check_seed_mismatch<IsRC>(match_pos_record, query, seed_mismatch_record, last_match_pos_record, mismatch_idx, tail_pos_record);
      
      if (!seed_mismatch_record.empty())
      {
        // take the mismatch index that has the shortest tail // 
        auto min_tail_pos = std::min_element(tail_pos_record.cbegin(), tail_pos_record.cend(), [](auto&& l, auto&& r){ return l <= r; });
        auto index_of_min = std::distance(tail_pos_record.cbegin(), min_tail_pos);
        //for (std::size_t i = 0; i < seed_mismatch_record.size(); ++i)
        //{
        tail_pos_record = { tail_pos_record[index_of_min] };
        seed_mismatch_record = { seed_mismatch_record[index_of_min] };
        last_match_pos_record = { last_match_pos_record[index_of_min] };
        //*******************//

        std::vector<Mismatchtype> tail_mismatch_record;
        std::vector<FMIdxRange> last_search_range;
        auto tail_pos_tmp = tail_pos_record.front();
        check_tail_mismatch<IsRC>(query, last_match_pos_record.front(), tail_pos_tmp, tail_mismatch_record, last_search_range, tail_pos_record);
        
        if (tail_mismatch_record.empty())
        {
          record_res_seed_mm<IsRC>(
            fq, 
            last_match_pos_record.front(), 
            tail_pos_tmp, 
            res, 
            seed_mismatch_record.front()
          );
        }
        
        for (std::size_t i = 0; i < tail_mismatch_record.size(); ++i)
        {
          record_res_seed_mm<IsRC>(
            fq,
            last_search_range[i],
            tail_pos_record[i],
            res,
            seed_mismatch_record.front(),
            tail_mismatch_record[i]
          );
        }
      }
    }
    else
    {
      int tail_pos = query.size() - prefix_match_len - 1;
      std::vector<int> tail_pos_record;
      std::vector<Mismatchtype> tail_mismatch_record;
      std::vector<FMIdxRange> last_search_range;
      
      if (para_pack.allow_mm)
        check_tail_mismatch<IsRC>(query, match_pos_record.back(), tail_pos, tail_mismatch_record, last_search_range, tail_pos_record);
      
      if (tail_mismatch_record.empty())
      {
        record_res<IsRC>(
          fq,
          match_pos_record.back(),
          tail_pos,
          res
        ); 
      }

      for (std::size_t i = 0; i < tail_mismatch_record.size(); ++i)
      {
        record_res<IsRC>(
          fq,
          last_search_range[i],
          tail_pos_record[i],
          res,
          tail_mismatch_record[i]
        );
      }
    }
    return;
  }

  template <bool IsRC>
  void check_seed_mismatch (
    const std::vector<FMIdxRange>& match_record,
    const SEQ& query,
    std::vector<SeedMMType>& seed_mismatch_record,
    std::vector<FMIdxRange>& last_match_pos_record,
    IndexType mismatch_idx,
    std::vector<int>& tail_pos_record
  ) const
  {
    auto match_rng = match_record.rbegin();
    FMIdxRange rng;
    for (auto i = mismatch_idx; i < query.size(); ++i)
    {
      if (i == 0)
          continue;
      char c;
      for (auto j = 0; j < ch_set.size(); ++j)
      {
        c = ch_set[j];
        
        if (query[i] == c) 
            continue;
        // third para pass rng and return bool //
        if constexpr (IsRC)
        {
          rng.first = rc_fm_index.lf_mapping ((*match_rng).first, c);
          rng.second = rc_fm_index.lf_mapping ((*match_rng).second, c);   
        }
        else
        {
          rng.first = fm_index.lf_mapping ((*match_rng).first, c);
          rng.second = fm_index.lf_mapping ((*match_rng).second, c);    
        }
        if (rng.first >= rng.second) 
            continue;
        
        base_search<IsRC>(query, i, rng, c, seed_mismatch_record, last_match_pos_record, tail_pos_record);
      }
      ++match_rng;
    }
  }

  template <bool IsRC>
  void base_search(
    const SEQ& query,
    std::uint32_t start_pos,
    FMIdxRange& rng,
    char mismatch_char,
    std::vector<SeedMMType>& seed_mismatch_record,
    std::vector<FMIdxRange>& last_match_pos_record,
    std::vector<int>& tail_pos_record
  ) const
  {
    auto idx = start_pos - 1;    
    FMIdxRange previous_rng;
    std::vector<IndexType> real_pos;
    while (1)
    {
      previous_rng = rng;
            
      if constexpr (IsRC)
      {
        rng.first = rc_fm_index.lf_mapping (rng.first, query[idx]);
        rng.second = rc_fm_index.lf_mapping (rng.second, query[idx]);
      }
      else
      {
        rng.first = fm_index.lf_mapping (rng.first, query[idx]);
        rng.second = fm_index.lf_mapping (rng.second, query[idx]);
      }

      if (rng.first >= rng.second)
      {
        if (idx < query.size() - para_pack.min_prefix_len)
        {
          if constexpr (IsRC)
              real_pos = rc_fm_index.sbwt_range_to_seq_idx(previous_rng, query.size() - idx - 1);
          else
              real_pos = fm_index.sbwt_range_to_seq_idx(previous_rng, query.size() - idx - 1);
          
          for (auto&& pos : real_pos)
          {
            seed_mismatch_record.emplace_back(std::make_tuple(pos, start_pos, mismatch_char));
            tail_pos_record.push_back(idx);
            last_match_pos_record.push_back(previous_rng);
          }
        }
        break;
      }

      if (idx == 0)
      { 
          if constexpr (IsRC)
              real_pos = rc_fm_index.sbwt_range_to_seq_idx(rng, query.size());
          else
              real_pos = fm_index.sbwt_range_to_seq_idx(rng, query.size());

          for (auto&& pos : real_pos)
          {
            seed_mismatch_record.emplace_back(std::make_tuple(pos, start_pos, mismatch_char));
            tail_pos_record.emplace_back(-1);
            last_match_pos_record.push_back(rng);
          }

          break;
      }
      --idx;
    }
  } 
  
  template <bool IsRC>
  void check_tail_mismatch(
    const SEQ& query,
    const FMIdxRange& begin_rng,
    int tail_begin_pos,
    std::vector<Mismatchtype>& tail_mismatch_record,
    std::vector<FMIdxRange>& last_search_range,
    std::vector<int>& tail_pos_record
  ) const
  {
    tail_pos_record.clear();
    if (tail_begin_pos < 2)   // or tail_begin_pos == -1 is contained
        return;

    int next_pos;
    auto origin_pos = tail_begin_pos;
    FMIdxRange rng;
    FMIdxRange previous_rng;
    char c;
    for (auto i = 0; i < 4; ++i)
    {
      c = ch_set[i];
      rng = begin_rng;
      next_pos = origin_pos - 1;

      if (query[origin_pos] == c) 
          continue;
      
      if constexpr (IsRC)
      {
        rng.first = rc_fm_index.lf_mapping (rng.first, c);
        rng.second = rc_fm_index.lf_mapping (rng.second, c);
      }
      else
      {
        rng.first = fm_index.lf_mapping (rng.first, c);
        rng.second = fm_index.lf_mapping (rng.second, c); 
      }

      if (rng.first >= rng.second) 
          continue;
      
      while (1)
      {
        previous_rng = rng;
        
        if constexpr (IsRC)
        {
          rng.first = rc_fm_index.lf_mapping (rng.first, query[next_pos]);
          rng.second = rc_fm_index.lf_mapping (rng.second, query[next_pos]);
        }
        else
        {
          rng.first = fm_index.lf_mapping (rng.first, query[next_pos]);
          rng.second = fm_index.lf_mapping (rng.second, query[next_pos]); 
        }
        
        if (rng.first >= rng.second)
        {
          if (next_pos == origin_pos - 1)
              break;
          
          tail_pos_record.push_back(next_pos);
          last_search_range.push_back(previous_rng);
          tail_mismatch_record.emplace_back(std::make_pair(origin_pos, c));

          break;
        }

        if (next_pos == 0)
        {
          tail_pos_record.push_back(-1);
          last_search_range.push_back(rng);
          tail_mismatch_record.emplace_back(std::make_pair(origin_pos, c));
          break;
        }
        --next_pos;
      }
    }
  }  

  template <bool IsRC = false>
  void record_res (
    const Fastq& fq,
    const FMIdxRange& rng,
    int tail_pos,
    std::vector<AlignedReads>& res,
    const Mismatchtype& tail_mismatch = std::make_pair(-1, ' ')
  ) const
  {
    std::uint16_t flag;
    auto genome_len = fm_index.seg_info.back().first + fm_index.seg_info.back().second;
    std::vector<IndexType> real_pos;
    std::vector<std::string> chrs;
    IndexType prefix_match_len = fq.seq.size() - tail_pos - 1;
    
    if constexpr (IsRC)
    {
      flag = 0;
      real_pos = rc_fm_index.sbwt_range_to_seq_idx(rng, prefix_match_len);
      if (para_pack.min_multi != 0 && real_pos.size() > para_pack.min_multi) 
          return;
      for (auto& i : real_pos)
      {
        i = genome_len - i - prefix_match_len + 1;
        chrs.emplace_back(find_chr(i));
      }
    }
    else
    {
      flag = 16;
      real_pos = fm_index.sbwt_range_to_seq_idx(rng, prefix_match_len);
      if (para_pack.min_multi != 0 && real_pos.size() > para_pack.min_multi)
          return;
      for (auto& i : real_pos)
      {
        i += 1;
        chrs.emplace_back(find_chr(i));
      }
    }
   
    if (real_pos.empty())
        return;
     
    bool has_tail_mismatch;
    
    if (!(tail_mismatch.first == -1 && tail_mismatch.second == ' '))
      has_tail_mismatch = true;
    else
      has_tail_mismatch = false;

    res.emplace_back( AlignedReads{
      true,
      fq,
      flag,
      tail_pos,
      rng,
      real_pos,
      chrs,
      rng.second - rng.first,
      prefix_match_len,
      false,
      has_tail_mismatch,
      std::make_pair(-1, ' '),
      tail_mismatch
    });
  }

  template <bool IsRC = false>
  void record_res_seed_mm (
    const Fastq& fq,
    const FMIdxRange& rng,
    int tail_pos,
    std::vector<AlignedReads>& res,
    const SeedMMType& seed_mismatch,
    const Mismatchtype& tail_mismatch = std::make_pair(-1, ' ')
  ) const
  {
    bool has_tail_mismatch;
    if (!(tail_mismatch.first == -1 && tail_mismatch.second == ' '))
      has_tail_mismatch = true;
    else
      has_tail_mismatch = false;
    
    std::uint16_t flag;
    auto genome_len = fm_index.seg_info.back().first + fm_index.seg_info.back().second;
    std::vector<IndexType> real_pos;
    std::vector<std::string> chrs;
    IndexType prefix_match_len = fq.seq.size() - tail_pos - 1;
    
    if constexpr (IsRC)
    {
      flag = 0;
      if (has_tail_mismatch)
      {
        real_pos = rc_fm_index.sbwt_range_to_seq_idx(rng, prefix_match_len);
        if (para_pack.min_multi != 0 && real_pos.size() > para_pack.min_multi) 
          return;
      }
      else
      {
        real_pos.push_back(std::get<0>(seed_mismatch));
      }
      
      for (auto& i : real_pos)
      {
        i = genome_len - i - prefix_match_len + 1;
        chrs.emplace_back(find_chr(i));
      }
    }
    else
    {
      flag = 16;
      if (has_tail_mismatch)
      {
        real_pos = fm_index.sbwt_range_to_seq_idx(rng, prefix_match_len);
        if (para_pack.min_multi != 0 && real_pos.size() > para_pack.min_multi)
          return;
      }
      else
      {
        real_pos.push_back(std::get<0>(seed_mismatch));
      }
      
      for (auto& i : real_pos)
      {
        i += 1;
        chrs.emplace_back(find_chr(i));
      }
    }
   
    if (real_pos.empty())
        return;

    res.emplace_back( AlignedReads{
      true,
      fq,
      flag,
      tail_pos,
      rng,
      real_pos,
      chrs,
      rng.second - rng.first,
      prefix_match_len,
      true,
      has_tail_mismatch,
      std::make_pair(std::get<1>(seed_mismatch), std::get<2>(seed_mismatch)),
      tail_mismatch
    });
  }

  auto find_chr(IndexType& real_pos) const
  {
    auto comp(
      [](const std::pair<IndexType, IndexType>& a, const IndexType& b)
      {
        return (a.first + a.second) <= b;
      });
    auto it(std::lower_bound(fm_index.seg_info.cbegin(), fm_index.seg_info.cend(), real_pos, comp));
    real_pos -= (*it).first;
    return fm_index.chr_names[std::distance(fm_index.seg_info.cbegin(), it)];
  }
};
}
