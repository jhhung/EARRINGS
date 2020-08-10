#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <numeric>
#include <vector>
//#include "/home/noreason/lab_project/PEAT/include/libsimdpp/simdpp/simd.h"

namespace biovoltron::string_sorter
{

template <typename SuffixArrayType = std::vector<uint32_t>>
class BucketSort
{
	//using SIMD_Seq = simdpp::uint8<16>;
	//using SIMD_SA = simdpp::uint32<4>;
	using IndexType = typename SuffixArrayType::value_type;
	using BucketType = std::vector<std::vector<IndexType>>;

	BucketType bucket_;

  public:
	explicit BucketSort(size_t bucket_size)
		: bucket_{
				std::vector<IndexType>(), 
				std::vector<IndexType>(bucket_size), 
				std::vector<IndexType>(bucket_size), 
				std::vector<IndexType>(bucket_size), 
				std::vector<IndexType>(bucket_size)
			}
	{
	}

	template <typename SEQ, typename Mapper>
	explicit BucketSort(
		size_t bucket_size, 
		SEQ&& seq, 
		SuffixArrayType& sa, 
		Mapper&& char_to_order, 
		IndexType prefix_len = 0
	)
		: bucket_{
				std::vector<IndexType>(), 
				std::vector<IndexType>(bucket_size), 
				std::vector<IndexType>(bucket_size), 
				std::vector<IndexType>(bucket_size), 
				std::vector<IndexType>(bucket_size)
			}
	{
		sort(std::forward<SEQ>(seq), sa, char_to_order, prefix_len);
	}

	template <typename SEQ, typename Mapper>
	void sort(
		SEQ&& seq, 
		SuffixArrayType& sa, 
		Mapper&& char_to_order, 
		IndexType prefix_len = 0)
	{
		if (sa.size() == 0)
		{
			sa.resize(seq.size());
			std::iota(sa.begin(), sa.end(), 0);
		}

		bucket_sort(
			seq, sa.begin(), sa.end(), char_to_order, 
			0, prefix_len);
	}

	template <typename SEQ, typename Mapper>
	void bucket_sort(
		SEQ&& seq, 
		typename SuffixArrayType::iterator begin, 
		typename SuffixArrayType::iterator end, 
		Mapper&& char_to_order, 
		IndexType offset, 
		IndexType prefix_len)
	{
		if (std::distance(begin, end) <= 1)
			return;

		std::array<IndexType, 4+1> count;
		
		count.fill(0);
		if (prefix_len == 0)
			prefix_len = seq.size();
		bucket_[0].resize(
			std::min(
				(size_t)prefix_len, 
				bucket_[1].size()
			));

		for (auto it(begin); it != end; it++)
		{
			IndexType tmp(*it);

			if (tmp + offset >= seq.size())
			{
				bucket_[0][count[0]++] = tmp;
				continue;
			}
			
			IndexType order(char_to_order[seq[tmp + offset]] + 1);

			bucket_[order][count[order]++] = tmp;
		}

		auto sa_it(begin);
		for (size_t i(0); i < count.size(); i++)
		{
			auto tmp(count[i]);

			if (tmp != 0)
			{
				std::copy(
					bucket_[i].begin(), 
					bucket_[i].begin() + tmp, 
					sa_it
				);
			}
			sa_it += tmp;
		}

		if (offset == prefix_len)
			return;

		sa_it = begin;
		for (size_t i(0); i < count.size(); i++)
		{
			auto tmp(count[i]);

			if (tmp != 0)
			{
				bucket_sort(
					seq, sa_it, sa_it + tmp, char_to_order, 
					offset + 1, prefix_len
				);
			}
			sa_it += tmp;
		}
	}

	inline size_t get_bucket_size()
	{
		return bucket_[1].size();
	}
};

}
