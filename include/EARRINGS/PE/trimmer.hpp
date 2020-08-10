#pragma once

#define FLANKING_BASE 16

#define FLANKING_BYTE (FLANKING_BASE / 4)
#define MAX_BASE 256
#define ALIGN_BYTE (ALIGN_BASE / 4)
#define ALIGN_BASE 64

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <boost/mpl/string.hpp>
#include <simdpp/simd.h>
#include <algorithm>

using SINGLE = boost::mpl::string<'SING', 'LE'>;
using PAIRED = boost::mpl::string<'PAIR', 'ED'>;

namespace EARRINGS
{

template <typename TRIMMER_TYPE, typename FORMAT>
class Trimmer;

template <typename FORMAT>
class Trimmer <PAIRED, FORMAT>
{
  public:
    using SIMD_Vector = simdpp::uint8<ALIGN_BYTE>;

  private:
    static constexpr auto ones_in_byte = []() constexpr
    {
	std::array<uint8_t, 256> ret{};
	
	for (size_t i(0); i < ret.size(); i++)
	{
	    uint8_t i_copy(i);

	    do
	    {
		if ((i_copy & 3) == 3)
		    ret[i]++;
		i_copy = i_copy >> 2;
	    }
	    while (i_copy);
	}

	return ret;
    }();
    static SIMD_Vector extract_mask;
    static std::vector<SIMD_Vector> erase_mask;

  public:
    struct TraitParm
    {
        float match_rate = 0.7, seq_cmp_rate = 0.9, adapter_cmp_rate = 0.8;
        size_t adapter_min_len = 5;

        TraitParm() = default;
        TraitParm (float match, float seq_cmp, float adap_cmp)
            : match_rate(match)
            , seq_cmp_rate(seq_cmp)
            , adapter_cmp_rate(adap_cmp)
        {
        }
    };

    TraitParm trait_parm;

    inline void cut_off_longer_seq(FORMAT& fq1, FORMAT& fq2) const
    {
        if (fq1.seq.size() > fq2.seq.size())
            fq1.trim(fq2.seq.size());
        else if (fq1.seq.size() < fq2.seq.size())
            fq2.trim(fq1.seq.size());
    }
         
    // find reverse complement position
    // RC check
    template<bool Check_Adapter=false, class SEQ>
    size_t find_rc_pos(FORMAT& fq1, FORMAT& fq2, std::vector<SEQ>& adapters) const
    {
        std::vector<size_t> possible_pos1, possible_pos2, 
            intersect_pos;

        auto rc1 = get_rc(fq1.seq, FLANKING_BASE);
        auto rc2 = get_rc(fq2.seq, FLANKING_BASE);
        
        find_possible_pos(possible_pos1, fq1.seq, rc2);
        find_possible_pos(possible_pos2, fq2.seq, rc1);
    
        for (auto it1(possible_pos1.begin()), 
            it2(possible_pos2.begin());
            it1 != possible_pos1.end() && 
            it2 != possible_pos2.end(); )
        {
            if (*it1 > *it2)
                it2++;
            else if (*it1 < *it2)
                it1++;
            else
            {
                intersect_pos.emplace_back(*it1);
                it1++;
                it2++;
            }
        }

        size_t trim_pos(fq1.seq.size());
        size_t idx(0);

        // std::cerr << "intersec_pos: " << intersect_pos.size() << "\n";
        size_t min_len1(0), min_len2(0);
        SEQ subseq1, subseq2, r_seq;
        for (auto r_it(intersect_pos.rbegin()); 
            r_it != intersect_pos.rend(); 
            r_it++)
        {
            idx = *r_it;

            r_seq = get_rc(fq2.seq, idx);
            if ((float)cal_match_num(fq1.seq, r_seq, idx) / idx
                    >= trait_parm.seq_cmp_rate)
            {
                if constexpr (Check_Adapter)
                {
                    if (idx == fq1.seq.size())
                    {
                        continue;
                    }
                    
                    subseq1 = SEQ(fq1.seq.begin() + idx, fq1.seq.end());
                    subseq2 = SEQ(fq2.seq.begin() + idx, fq2.seq.end());
                    min_len1 = std::min(subseq1.size(), adapters[0].size());
                    min_len2 = std::min(subseq2.size(), adapters[1].size());
                    if ((float)cal_match_num(subseq1, adapters[0], min_len1) 
                        / min_len1 >= trait_parm.adapter_cmp_rate ||
                        (float)cal_match_num(subseq2, adapters[1], min_len2)
                        / min_len2 >= trait_parm.adapter_cmp_rate)
                    {
                        return idx;
                    }
                }
                else
                {
                    return idx;
                }
                
            }
        }
        
        possible_pos1.clear();
        possible_pos2.clear();
        
        // if no adapters are detected in the sequence,
        // check if adapters are located at the head of the reads
        if (Check_Adapter && intersect_pos.size() == 0)
        {
            subseq1 = SEQ(fq1.seq.begin(), fq1.seq.begin() + FLANKING_BASE * 2);
            size_t pos = find_adapter_pos(subseq1, adapters[0]);
            if (pos != -1)
            {
                return pos;
            }
            subseq2 = SEQ(fq2.seq.begin(), fq2.seq.begin() + FLANKING_BASE * 2);
            pos = find_adapter_pos(subseq2, adapters[1]);
            if (pos != -1)
            {
                return pos;
            }
        }
        else;  // do nothing
        

        return trim_pos;
    }


  private:
    template <typename SEQ>
    SEQ get_rc(const SEQ& seq, const size_t size) const
    {
        SEQ reverse_seq(size);

        for (size_t i(0); i < reverse_seq.size(); i++)
            reverse_seq[i] = seq[reverse_seq.size() - 1 - i];
        
        reverse_seq.flip();
        return reverse_seq;
    }

    template<typename SEQ>
    size_t find_adapter_pos(SEQ &seq, SEQ &seq_adapt) const
    {
        size_t min_len = std::min(seq.size(), (size_t)(FLANKING_BASE * 2));
        SEQ subseq(seq.begin(), seq.begin() + min_len);
        SIMD_Vector adapter(simdpp::load(seq_adapt[0].get_seg()))
                  , seq_buf(simdpp::load(subseq[0].get_seg()));
        
        size_t match_num;
        for (size_t i(0); i < FLANKING_BASE; ++i)
        {
            min_len = std::min({subseq.size() - i, (size_t)FLANKING_BASE, seq_adapt.size()});
            if (min_len < trait_parm.adapter_min_len)
                break;
            auto &mask = erase_mask[min_len];

            match_num = fixed_len_cal_match_num<ALIGN_BYTE>(
                        ~((seq_buf & mask) ^ (adapter | ~mask)), 
                        std::make_index_sequence<ALIGN_BYTE>{}  
                    );

            if ((float)match_num / min_len > trait_parm.adapter_cmp_rate) 
                return i;

            seq_buf = (seq_buf >> 2) | 
                simdpp::move16_l<1>(
                (seq_buf & extract_mask) << 6
                );
        }

        return -1;
    }

    template <typename SEQ>
    void find_possible_pos(std::vector<size_t>& possible_pos, 
	    const SEQ& seq, const SEQ& seq_adapt) const
    {
        SIMD_Vector adapter(simdpp::load(seq_adapt[0].get_seg()));
        SIMD_Vector seq_buf, next_seq_buf;
        size_t remaining_base;
        size_t min_len = std::min({seq.size(), (size_t)FLANKING_BASE, seq_adapt.size()});

        // if sequence size is smaller than the align base
        if (seq.size() < ALIGN_BASE)
        {
            seq_buf = simdpp::load(seq[0].get_seg());
            for (auto i(0); i < FLANKING_BASE; ++i)
            {
                min_len = std::min({seq.size() - i, (size_t)FLANKING_BASE, seq_adapt.size()});
                if (min_len < trait_parm.adapter_min_len) break;
                auto& mask = erase_mask[min_len];
                auto match_num = fixed_len_cal_match_num<ALIGN_BYTE>(
                    ~((seq_buf & mask) ^ (adapter | ~mask)), 
                    std::make_index_sequence<FLANKING_BASE>{}  
                    );
                if ((float)match_num / min_len > trait_parm.match_rate)
                {
                    possible_pos.emplace_back(i);
                }

                seq_buf = (seq_buf >> 2) | 
                    simdpp::move16_l<1>(
                    (seq_buf & extract_mask) << 6
                    );
            }
        }
        
        auto& mask = erase_mask[min_len];
        for (size_t i(0); i < seq.size() / ALIGN_BASE; i++)
        {
            seq_buf = simdpp::load(
                seq[i * ALIGN_BASE].get_seg()
            );
            next_seq_buf = simdpp::load(
                seq[(i + 1) * ALIGN_BASE].get_seg()
            );

            if (i == seq.size() / ALIGN_BASE - 1)
            remaining_base = seq.size() - (i + 1) * ALIGN_BASE 
                + 1 + (ALIGN_BASE - FLANKING_BASE);
            else
            remaining_base = ALIGN_BASE;

            for (size_t j(0); j < remaining_base; j++)
            {
                auto match_num = fixed_len_cal_match_num<ALIGN_BYTE>(
                    ~((seq_buf & mask) ^ (adapter | ~mask)), 
                    std::make_index_sequence<FLANKING_BYTE>{}  
                    );
                if ((float)match_num / min_len >= 
                    trait_parm.match_rate)
                {
                    possible_pos.emplace_back(
                        i * ALIGN_BASE + j + FLANKING_BASE
                    );
                }

                seq_buf = (seq_buf >> 2) | 
                    simdpp::move16_l<1>(
                    (seq_buf & extract_mask) << 6
                    );
                seq_buf = simdpp::insert<ALIGN_BYTE - 1>(
                    seq_buf, 
                    simdpp::extract<ALIGN_BYTE - 1>(seq_buf) | 
                    (simdpp::extract<0>(next_seq_buf) & 3) << 6
                    );
                next_seq_buf = (next_seq_buf >> 2) | 
                    simdpp::move16_l<1>(
                    (next_seq_buf & extract_mask) << 6
                    );
            }
        }
    }

    template <uint32_t I, size_t... IDX>
    uint32_t fixed_len_cal_match_num(const simdpp::uint8<I>& cmp, 
	    std::index_sequence<IDX...>) const
    {
	uint32_t count(0);

	((count += ones_in_byte[simdpp::extract<IDX>(cmp)]), ...);

	return count;
    }

    template <typename SEQ>
    uint32_t cal_match_num(const SEQ& seq1, const SEQ& seq2, size_t idx) const
    {
        if (seq1.empty() || seq2.empty())
        {
            return 0;
        }
        
		uint32_t match_num(0);
		SIMD_Vector cmp_buf1, cmp_buf2;

		for (size_t i(0); i < idx / ALIGN_BASE; i++)
		{
			cmp_buf1 = simdpp::load(seq1[i * ALIGN_BASE].get_seg());
            cmp_buf2 = simdpp::load(seq2[i * ALIGN_BASE].get_seg());

            match_num += fixed_len_cal_match_num<ALIGN_BYTE>(
                ~(cmp_buf1 ^ cmp_buf2),
                std::make_index_sequence<ALIGN_BYTE>{}
                );
		}

		auto& mask(
			erase_mask[idx - idx / ALIGN_BASE * ALIGN_BASE]
			);

        // idx-1
		cmp_buf1 = simdpp::load(
			seq1[idx / ALIGN_BASE * ALIGN_BASE].get_seg()
			);

        cmp_buf2 = simdpp::load(
            seq2[idx / ALIGN_BASE * ALIGN_BASE].get_seg()
        );
        match_num += fixed_len_cal_match_num<ALIGN_BYTE>(
            ~((cmp_buf1 & mask) ^ (cmp_buf2 | ~mask)), 
            std::make_index_sequence<ALIGN_BYTE>{}
        );
        
        return match_num;
    }

};

template <typename FORMAT>
simdpp::uint8<ALIGN_BYTE> Trimmer<PAIRED, FORMAT>::extract_mask 
    = simdpp::make_uint(3);

template <typename FORMAT>
std::vector<simdpp::uint8<ALIGN_BYTE>> 
    Trimmer<PAIRED, FORMAT>::erase_mask = [](){
	std::vector<simdpp::uint8<ALIGN_BYTE>> v(ALIGN_BASE);
	std::vector<uint64_t> mask_val(ALIGN_BYTE, 0xFF);

	for (auto i(ALIGN_BASE - 1); i >= 0; i--)
	{
	    mask_val[i * 2 / 8] >>= 2;
	    v[i] = simdpp::make_uint(
		    mask_val[0], mask_val[1], mask_val[2], 
		    mask_val[3], mask_val[4], mask_val[5], 
		    mask_val[6], mask_val[7], mask_val[8], 
		    mask_val[9], mask_val[10], mask_val[11], 
		    mask_val[12], mask_val[13], mask_val[14], 
		    mask_val[15]
		);
	}
	
	return v;
    }();
}