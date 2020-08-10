#pragma once

#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <cstdint>
#include <fstream>

template <typename T> 
inline static int32_t sign(T a) 
{
	return (a > T(0)) - (a < T(0));
}

template <typename SEQ, typename IndexType>
bool suffix_less(
	SEQ&& seq, 
	IndexType a, 
	IndexType b, 
	IndexType prefix_len, 
	IndexType offset = 0)
{
	if (a == b)
		return false;
	
	for (IndexType i(offset), 
		times(
			std::min(
				prefix_len - 1, 
				(IndexType)seq.size() - 1 - std::max(a, b)
			)); 
		i <= times;
		i++)
	{
		auto sign_val(sign(seq[b + i] - seq[a + i]));

		switch(sign_val)
		{
			case 1:
				return true;
			case -1:
				return false;
		}
	}

	return (a > b ? true : false);
}

inline auto make_archive(uint32_t count)
{
	std::pair<
		std::vector<std::unique_ptr<std::ofstream>>, 
		std::vector<
			std::unique_ptr<
				boost::archive::binary_oarchive
			>>
	> out;
	
	out.first.emplace_back(
		std::make_unique<std::ofstream>(
			"group_" + std::to_string(count), 
	 		std::ios::binary | std::ios::out
		));
	out.second.emplace_back(
		std::make_unique<boost::archive::binary_oarchive>(
			*(out.first.back())
		));
	
	return out;
}
