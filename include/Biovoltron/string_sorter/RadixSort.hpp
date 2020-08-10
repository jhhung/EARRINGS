#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <memory>
#include <vector>
#include <Biovoltron/string_sorter/utility.hpp>
//#include "/home/noreason/lab_project/PEAT/include/libsimdpp/simdpp/simd.h"

namespace biovoltron::string_sorter
{

template <typename SuffixArrayType = std::vector<uint32_t>>
class RadixSort
{
	//using SIMD_Seq = simdpp::uint8<16>;
	//using SIMD_SA = simdpp::uint32<4>;
	using IndexType = typename SuffixArrayType::value_type;
	using BucketType = std::vector<std::vector<IndexType>>;

	BucketType bucket_;

  public:
	explicit RadixSort(size_t bucket_size)
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
	explicit RadixSort(
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

		radix_sort(
			std::forward<SEQ>(seq), sa.begin(), sa.end(), 
			std::forward<Mapper>(char_to_order), 0, prefix_len
		);
	}

	template <typename SEQ, typename Mapper>
	void radix_sort(
		SEQ&& seq, 
		typename SuffixArrayType::iterator begin, 
		typename SuffixArrayType::iterator end, 
		Mapper&& char_to_order, 
		IndexType offset, 
		IndexType prefix_len)
	{
		if (std::distance(begin, end) <= 1)
			return;
		
		if (std::distance(begin, end) < 8)
		{
			std::sort(begin, end, 
				[&seq, offset, prefix_len](IndexType a, IndexType b)
				{
					return suffix_less(seq, a, b, prefix_len, offset);
				}
			);

			return;
		}

		if (prefix_len == 0)
			prefix_len = seq.size();

		bucket_[0].resize(
			std::min(
				(size_t)prefix_len, 
				bucket_[1].size()
			));

		radix_sort_impl(
			std::forward<SEQ>(seq), begin, end, 
			std::forward<Mapper>(char_to_order), offset, prefix_len
		);
	}

	inline void set_bucket_size(size_t s)
	{
		for (size_t i(1); i < bucket_.size(); i++)
			bucket_[i].resize(s);
	}

	inline size_t get_bucket_size()
	{
		return bucket_[1].size();
	}

  private:
	
	template <typename SEQ, typename Mapper>
	void radix_sort_impl(
		SEQ&& seq, 
		typename SuffixArrayType::iterator begin, 
		typename SuffixArrayType::iterator end, 
		Mapper&& char_to_order, 
		IndexType offset, 
		IndexType prefix_len)
	{
		if (std::distance(begin, end) <= 1)
			return;
		
		if (std::distance(begin, end) < 8)
		{
			std::sort(begin, end, 
				[&seq, offset, prefix_len](IndexType a, IndexType b)
				{
					return suffix_less(seq, a, b, prefix_len, offset);
				}
			);

			return;
		}

		std::array<IndexType, 4+1> count;
		
		count.fill(0);

		IndexType tmp, order;
		for (auto it(begin); it != end; it++)
		{
			tmp = *it;

			if (tmp + offset >= seq.size())
			{
				bucket_[0][count[0]++] = tmp;
				continue;
			}
			
			order = char_to_order[seq[tmp + offset]] + 1;
			bucket_[order][count[order]++] = tmp;
		}

		auto sa_it(begin);
		// for (size_t i(0); i < count.size(); i++)
		// {
		// 	tmp = count[i];

		// 	if (tmp != 0)
		// 	{
		// 		std::copy(
		// 			bucket_[i].begin(), 
		// 			bucket_[i].begin() + tmp, 
		// 			sa_it
		// 		);
		// 	}
		// 	sa_it += tmp;
		// }

		auto* p(std::addressof(*begin));
		tmp = count[0];
		if (tmp != 0)
		{
			std::memcpy(p, bucket_[0].data(), tmp*4);
		}
		p += tmp;
		tmp = count[1];
		if (tmp != 0)
		{
			std::memcpy(p, bucket_[1].data(), tmp*4);
		}
		p += tmp;
		tmp = count[2];
		if (tmp != 0)
		{
			std::memcpy(p, bucket_[2].data(), tmp*4);
		}
		p += tmp;
		tmp = count[3];
		if (tmp != 0)
		{
			std::memcpy(p, bucket_[3].data(), tmp*4);
		}
		p += tmp;
		tmp = count[4];
		if (tmp != 0)
		{
			std::memcpy(p, bucket_[4].data(), tmp*4);
		}
		p += tmp;

		if (offset + 1 == prefix_len)
			return;

		sa_it = begin;
		// for (size_t i(0); i < count.size(); i++)
		// {
		// 	tmp = count[i];

		// 	if (tmp != 0)
		// 	{
		// 		radix_sort_impl(
		// 			seq, sa_it, sa_it + tmp, char_to_order, 
		// 			offset + 1, prefix_len
		// 		);
		// 	}
		// 	sa_it += tmp;
		// }

		tmp = count[0];
		if (tmp != 0)
		{
			radix_sort_impl(
				seq, sa_it, sa_it + tmp, char_to_order, 
				offset + 1, prefix_len
			);
		}
		sa_it += tmp;
		tmp = count[1];
		if (tmp != 0)
		{
			radix_sort_impl(
				seq, sa_it, sa_it + tmp, char_to_order, 
				offset + 1, prefix_len
			);
		}
		sa_it += tmp;
		tmp = count[2];
		if (tmp != 0)
		{
			radix_sort_impl(
				seq, sa_it, sa_it + tmp, char_to_order, 
				offset + 1, prefix_len
			);
		}
		sa_it += tmp;
		tmp = count[3];
		if (tmp != 0)
		{
			radix_sort_impl(
				seq, sa_it, sa_it + tmp, char_to_order, 
				offset + 1, prefix_len
			);
		}
		sa_it += tmp;
		tmp = count[4];
		if (tmp != 0)
		{
			radix_sort_impl(
				seq, sa_it, sa_it + tmp, char_to_order, 
				offset + 1, prefix_len
			);
		}
		sa_it += tmp;
	}
};

}
