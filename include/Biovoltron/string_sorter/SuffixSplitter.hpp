#pragma once

#include <Biovoltron/string_sorter/utility.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>
#include <type_traits>
//#include "/home/noreason/lab_project/PEAT/include/libsimdpp/simdpp/simd.h"

namespace biovoltron::string_sorter
{

template <
	typename LinearSorter, 
	typename SuffixArrayType = std::vector<uint32_t>
>
class SuffixSplitter
{
	using IndexType = typename SuffixArrayType::value_type;
	//using SIMD_Seq = simdpp::uint8<16>;
	//using SIMD_SA = simdpp::uint32<4>;

	inline static std::mt19937 gen_{std::random_device{}()};
	uint32_t arch_count_;
	IndexType splitter_idx_;

	static constexpr auto non_zero_pos_ = []() constexpr
    {
		std::array<uint8_t, 256> ret{};
		
		for (size_t i(0); i < ret.size(); i++)
		{
			uint8_t i_copy(i);

			for (uint8_t j(0); j <= 4; j++)
			{
				if (i_copy == 0)
				{
					ret[i] = 4 - j;
					break;
				}

				i_copy = i_copy >> 2;
			}
		}

		return ret;
    }();

  public:
	LinearSorter sorter;
	SuffixArrayType splitters;
	IndexType partition_threshold;

	explicit SuffixSplitter(
		IndexType file_size_threshold, 
		IndexType bucket_size_threshold
	)
		: arch_count_(0)
		, splitter_idx_(0)
		, partition_threshold(file_size_threshold)
		, sorter(bucket_size_threshold)
	{
	}

	template <
		typename SEQ, 
		typename = std::enable_if_t<
				!std::is_integral_v<SEQ> && 
				!std::is_convertible_v<SEQ, SuffixSplitter>
			>
	>
	explicit SuffixSplitter(
		SEQ&& seq, 
		IndexType file_size_threshold, 
		IndexType bucket_size_threshold
	)
		: arch_count_(0)
		, splitter_idx_(0)
		, partition_threshold(file_size_threshold)
		, sorter(bucket_size_threshold)
	{
		sample(std::forward<SEQ>(seq), SuffixArrayType{});
	}

	template <
		typename SEQ, 
		typename = std::enable_if_t<
				!std::is_integral_v<SEQ> && 
				!std::is_convertible_v<SEQ, SuffixSplitter>
			>
	>
	explicit SuffixSplitter(
		SEQ&& seq, 
		const SuffixArrayType& sa, 
		IndexType file_size_threshold, 
		IndexType bucket_size_threshold
	)
		: arch_count_(0)
		, splitter_idx_(0)
		, partition_threshold(file_size_threshold)
		, sorter(bucket_size_threshold)
	{
		sample(std::forward<SEQ>(seq), sa);
	}

	template <typename SEQ, typename Mapper>
	void sample(
		SEQ&& seq, 
		Mapper&& char_to_order, 
		IndexType prefix_len = 0)
	{
		IndexType group_num(seq.size() / partition_threshold);

		splitter_idx_ = 0;
		splitters.reserve(group_num);
		
		if (group_num != 0)
		{
			std::uniform_int_distribution<IndexType> 
				uid(0, seq.size() - 1);
			std::vector<IndexType> tmp_v(group_num * 4);
			for (auto& i : tmp_v)
				i = uid(gen_);
			sorter.sort(seq, tmp_v, char_to_order, prefix_len);
			for (IndexType i(1); i <= group_num; i++)
				splitters.emplace_back(
					tmp_v[tmp_v.size() * i / (group_num + 1)]
				);
		}
	}

	template <typename SEQ>
	void split(
		SEQ&& seq, 
		IndexType prefix_len = 0)
	{
		SuffixArrayType found;

		if (splitters.size() == 0)
		{
			auto out(make_archive(0));

			found.resize(seq.size());
			std::iota(found.begin(), found.end(), 0);
			*(out.second.back()) & found;
			out.first.back()->close();
			return;
		}

		found.reserve(partition_threshold);
		
		while (splitter_idx_ <= splitters.size())
		{
			split_groups(seq, prefix_len, found);
			found.clear();
		}

		return;
	}

	template <typename SEQ>
	void split_groups(
		SEQ&& seq, 
		IndexType prefix_len, 
		SuffixArrayType& found
	)
	{
		SuffixArrayType empty{};
		if (splitter_idx_ > 0 && splitter_idx_ < splitters.size())
		{
			for (IndexType i(0); i < seq.size(); i++)
			{
				if (suffix_less(seq, i, splitters[splitter_idx_], 
						prefix_len) && 
					suffix_less(seq, splitters[splitter_idx_ - 1], 
						i, prefix_len))
				{
					found.emplace_back(i);
					if (if_larger_than_threshold(
							seq, prefix_len, empty, found, i + 1))
					{
						auto out(make_archive(arch_count_++));

						*(out.second.back()) & found;
						out.first.back()->close();
						return;
					}
				}
			}
		}
		else if (splitter_idx_ == 0)
		{
			for (IndexType i(0); i < seq.size(); i++)
			{
				if (suffix_less(seq, i, splitters[splitter_idx_], 
						prefix_len))
				{
					found.emplace_back(i);
					if (if_larger_than_threshold(
							seq, prefix_len, empty, found, i + 1))
					{
						auto out(make_archive(arch_count_++));

						*(out.second.back()) & found;
						out.first.back()->close();
						return;
					}
				}	
			}
		}
		else if (splitter_idx_ == splitters.size())
		{
			for (IndexType i(0); i < seq.size(); i++)
			{
				if (suffix_less(seq, splitters[splitter_idx_ - 1], 
						i, prefix_len))
				{
					found.emplace_back(i);
					if (if_larger_than_threshold(
							seq, prefix_len, empty, found, i + 1))
					{
						auto out(make_archive(arch_count_++));

						*(out.second.back()) & found;
						out.first.back()->close();
						return;
					}
				}
			}
		}
		else
		{
			return;
		}

		for (auto it(splitters.cbegin() + splitter_idx_ + 1); 
			it < splitters.cend(); 
			it++)
		{
			if (std::find(found.cbegin(), found.cend(), *it) 
				!= found.cend())
				splitters.erase(it--);
		}
		
		if (splitter_idx_ < splitters.size())
			found.emplace_back(splitters[splitter_idx_]);

		auto out(make_archive(arch_count_++));

		*(out.second.back()) & found;
		out.first.back()->close();
		splitter_idx_++;
		return;
	}

	/*
	template <typename SEQ>
	static bool suffix_less(
		SEQ&& seq, 
		IndexType a, 
		IndexType b, 
		IndexType prefix_len = 0)
	{
		if (a == b)
			return false;
		
		SIMD_Seq simd_seq1, simd_seq2;
		
		for (IndexType i(0), 
			times(std::min(
					prefix_len - 1, 
					(IndexType)seq.size() - 1 - std::max(a, b)
				) / 64); 
			i <= simd_times;
			i++)
		{
			simd_seq1 = simdpp::load(seq[a + i * 64].get_seg());
			simd_seq2 = simdpp::load(seq[b + i * 64].get_seg());
			simd_seq1 = simd_seq1 ^ simd_seq2;

			size_t pos;
			if ((pos = find_mismatch_pos(
					simd_seq1, std::make_index_sequence<16>{}
				)) != 64)
			{
				if (seq[a + i*64 + pos] < seq[b + i*64 + pos])
					return true;
				else if (seq[a + i*64 + pos] > seq[b + i*64 + pos])
					return false;
			}
		}

		if (a > b)
			return true;
		else
			return false;
	}

	template <typename SIMD, size_t... IDX>
	static inline size_t find_mismatch_pos(
		SIMD&& v, 
		std::index_sequence<IDX...>)
	{
		size_t result, count(0), tmp;

		(((tmp = non_zero_pos_[simdpp::extract<IDX>(v)]) == 4 ?
			count += 4 : result = count + tmp
		), ...);

		return result;
	}
	*/

	inline void reset(IndexType threshold)
	{
		splitter_idx_ = 0;
		arch_count_ = 0;
		splitters.clear();
		partition_threshold = threshold;
	}

  private:

	template <typename SEQ>
	inline bool if_larger_than_threshold(
		SEQ&& seq, 
		IndexType prefix_len, 
		SuffixArrayType& found, 
		SuffixArrayType& new_found, 
		IndexType idx)
	{
		if (new_found.size() == partition_threshold)
		{
			found.clear();
			found.shrink_to_fit();

			auto less_than = 
			[&seq, &new_found, prefix_len]
			(IndexType a, IndexType b)
			{
				return suffix_less(
						seq, new_found[a], new_found[b], prefix_len);
			};

			std::vector<IndexType> candidates(3), tmp(9);
			
			for (size_t i(1); i <= tmp.size(); i++)
				tmp[i] = new_found.size() * i / (9 + 1);
			std::sort(tmp.begin(), tmp.begin() + 3, less_than);
			std::sort(tmp.begin() + 3, tmp.begin() + 6, less_than);
			std::sort(tmp.begin() + 6, tmp.end(), less_than);
			candidates[0] = tmp[1];
			candidates[1] = tmp[4];
			candidates[2] = tmp[7];
			std::sort(
				candidates.begin(), 
				candidates.end(), 
				less_than
			);

			IndexType new_splitter(new_found[candidates[1]]);
			auto it = std::find(
				splitters.cbegin(), splitters.cend(), new_splitter);
			
			if (it != splitters.cend() && *it == new_splitter)
				splitters.erase(it);
			
			if (splitter_idx_ == splitters.size())
				splitters.insert(splitters.cend(), new_splitter);
			else
				splitters.insert(
					splitters.cbegin() + splitter_idx_, new_splitter);

			new_found = resplit(
					std::forward<SEQ>(seq), 
					prefix_len, 
					new_found, 
					idx
				);

			return true;
		}
		return false;
	}

	template <typename SEQ>
	SuffixArrayType resplit(
		SEQ&& seq, 
		IndexType prefix_len, 
		SuffixArrayType& found, 
		IndexType idx)
	{
		SuffixArrayType new_found;

		new_found.reserve(partition_threshold);
		
		if (splitter_idx_ > 0 && splitter_idx_ < splitters.size())
		{
			for (auto i : found)
			{
				if (suffix_less(seq, i, splitters[splitter_idx_], 
						prefix_len) && 
					suffix_less(seq, splitters[splitter_idx_ - 1], 
						i, prefix_len))
				{
					new_found.emplace_back(i);
					if (if_larger_than_threshold(
							seq, prefix_len, found, new_found, idx)
					)
						return new_found;
				}
			}

			for (IndexType i(idx); i < seq.size(); i++)
			{
				if (suffix_less(seq, i, splitters[splitter_idx_], 
						prefix_len) && 
					suffix_less(seq, splitters[splitter_idx_ - 1], 
						i, prefix_len))
				{
					new_found.emplace_back(i);
					if (if_larger_than_threshold(
							seq, prefix_len, found, new_found, i + 1)
					)
						return new_found;
				}
			}
		}
		else if (splitter_idx_ == 0)
		{
			for (auto i : found)
			{
				if (suffix_less(seq, i, splitters[splitter_idx_], 
						prefix_len))
				{
					new_found.emplace_back(i);
					if (if_larger_than_threshold(
							seq, prefix_len, found, new_found, idx)
					)
						return new_found;
				}
			}

			for (IndexType i(idx); i < seq.size(); i++)
			{
				if (suffix_less(seq, i, splitters[splitter_idx_], 
						prefix_len))
				{
					new_found.emplace_back(i);
					if (if_larger_than_threshold(
							seq, prefix_len, found, new_found, i + 1)
					)
						return new_found;
				}
			}
		}
		else if (splitter_idx_ == splitters.size())
		{
			for (auto i : found)
			{
				if (suffix_less(seq, splitters[splitter_idx_ - 1], 
						i, prefix_len))
				{
					new_found.emplace_back(i);
					if (if_larger_than_threshold(
							seq, prefix_len, found, new_found, idx)
					)
						return new_found;
				}
			}

			for (IndexType i(idx); i < seq.size(); i++)
			{
				if (suffix_less(seq, splitters[splitter_idx_ - 1], 
						i, prefix_len))
				{
					new_found.emplace_back(i);
					if (if_larger_than_threshold(
							seq, prefix_len, found, new_found, i + 1)
					)
						return new_found;
				}
			}
		}
		else
		{
			return new_found;
		}

		for (auto it(splitters.cbegin() + splitter_idx_ + 1); 
			it < splitters.cend(); 
			it++)
		{
			if (std::find(new_found.cbegin(), new_found.cend(), *it) 
				!= new_found.cend())
				splitters.erase(it--);
		}
		
		if (splitter_idx_ < splitters.size())
			new_found.emplace_back(splitters[splitter_idx_]);
		splitter_idx_++;

		return new_found;
	}
};

}
