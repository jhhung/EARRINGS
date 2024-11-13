#pragma once

#include <OldBiovoltron/string_sorter/SuffixSplitter.hpp>

#include <boost/archive/binary_iarchive.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace biovoltron::string_sorter
{

template <typename SuffixArrayType, typename Sorter>
class SBWT
{
	using IndexType = typename SuffixArrayType::value_type;
	uint32_t sort_count_;
	std::vector<IndexType> read_times_of_groups_;

  public:
	SuffixSplitter<Sorter> splitter;

	explicit SBWT(
		IndexType file_size_threshold, 
		IndexType bucket_size_threshold
	)
		: sort_count_(0)
		, splitter(file_size_threshold, bucket_size_threshold)
	{
	}

	template <typename SEQ, typename Mapper>
	SuffixArrayType sort_some(
		SEQ&& seq, 
		Mapper&& char_to_order, 
		IndexType prefix_len = 0)
	{
		if (sort_count_ == 0)
		{
			splitter.sample(seq, char_to_order, prefix_len);
			splitter.split(seq, prefix_len);
		}

		if (prefix_len == 0)
			prefix_len = seq.size();

		SuffixArrayType v, buf, tmp;
		std::ifstream ifs(
			"group_" + std::to_string(sort_count_), 
			std::ios::binary | std::ios::in);
		boost::archive::binary_iarchive i_arch(ifs);
		
		i_arch & v;
		split_sort(seq, v, v.begin(), v.end(), 
			char_to_order, 0, prefix_len);
		sort_count_++;

		return v;
	}

	inline void reset(IndexType threshold)
	{
		sort_count_ = 0;
		splitter.reset(threshold);
	}

  private:
	
	template <typename SEQ, typename Mapper>
	void split_sort(
		SEQ&& seq, 
		SuffixArrayType& v, 
		typename SuffixArrayType::iterator begin, 
		typename SuffixArrayType::iterator end, 
		Mapper&& char_to_order, 
		IndexType offset, 
		IndexType prefix_len)
	{
		if (begin + 1 >= end || offset == prefix_len)
			return;
		
		size_t len(std::distance(begin, end));

		if (len < splitter.sorter.get_bucket_size())
		{
			splitter.sorter.radix_sort(
				seq, begin, end, char_to_order, 0, prefix_len
			);
			return;
		}

		auto digit_comp = 
			[&seq, offset](IndexType a, IndexType b)
			{
				auto offset_a(a + offset), offset_b(b + offset);
				if (offset_a < seq.size())
				{
					if (offset_b < seq.size())
						return sign(seq[offset_a] - seq[offset_b]);
					else
						return 1;
				}
				else
				{
					if (offset_b < seq.size())
						return -1;
					else
						return 0;
				}
			};
		
		auto less_than = 
			[&seq, begin, prefix_len](IndexType a, IndexType b)
			{
				return suffix_less(
						seq, *(begin + a), *(begin + b), prefix_len);
			};

		std::vector<IndexType> candidates(3), tmp(9);

		for (size_t i(0); i < tmp.size(); i++)
			tmp[i] = (len - 1) * i / 8;
		std::sort(tmp.begin(), tmp.begin() + 3, less_than);
		std::sort(tmp.begin() + 3, tmp.begin() + 6, less_than);
		std::sort(tmp.begin() + 6, tmp.end(), less_than);
		candidates[0] = tmp[1];
		candidates[1] = tmp[4];
		candidates[2] = tmp[7];
		std::sort(candidates.begin(), candidates.end(), less_than);
		std::iter_swap(begin, begin + candidates[1]);

		auto it(begin + 1), l_it(begin + 1), r_it(end);
		while (it != r_it)
		{
			auto c = digit_comp(*it, *begin);
			switch (c)
			{
				case -1:
					std::iter_swap(it++, l_it++);
					break;
				case 1:
					std::iter_swap(it, --r_it);
					break;
				case 0:
					it++;
					break;
			}
		}
		std::iter_swap(begin, --l_it);
		
		split_sort(
			seq, v, begin, l_it, char_to_order, 
			offset, prefix_len);
		split_sort(
			seq, v, l_it, r_it, char_to_order, 
			offset + 1, prefix_len);
		split_sort(
			seq, v, r_it, end, char_to_order, 
			offset, prefix_len);
	}
};

}
