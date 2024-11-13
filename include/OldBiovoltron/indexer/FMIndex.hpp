#pragma once

#include <OldBiovoltron/format/fastq.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

static std::vector<unsigned long>
    fast_rand_seed = {123456789, 362436069, 521288629};

namespace biovoltron::indexer
{

template <
	typename SEQ, 
	typename SuffixArrayType, 
	typename Sorter, 
	uint32_t LogInterval, 
	uint32_t CharTypeNum, 
	uint32_t LookupStrLen, 
	uint32_t PrefixLen, 
	bool IsSBWT = false
>
class FMIndex
{
  public:
  
  static constexpr auto lookup_len = LookupStrLen;
  
  using BaseType = typename SEQ::value_type;
	using SizeType = typename SEQ::size_type;
	using SeqType = SEQ;
	using IndexType = typename SuffixArrayType::value_type;
	using FASTQ = biovoltron::format::FASTQ<SEQ>;
	using SamType = std::tuple<
			std::string, //QNAME
			uint16_t, //SAM_FLAG, //FLAG
			std::string, //RNAME
			uint32_t, //POS
			uint8_t , //MAPQ
			std::string, //CIGAR
			std::string, //RNEXT
			uint32_t, //PNEXT
			int32_t, //TLEN
			std::string, //SEQ
			std::string, //QUAL
			uint32_t
		>;

  private:
  
  std::array<IndexType, 256> char_to_order_;
	std::array<IndexType, CharTypeNum> order_to_char_;

  public:
	
  static constexpr uint32_t 
		interval{(uint32_t)(1 << LogInterval)};
	SEQ bwt_seq;
	IndexType seq_end_pos;
  std::vector<std::string> chr_names;
	std::vector<IndexType> c_table, lookup_table;
	std::vector<std::pair<IndexType, IndexType>> 
		seg_info, n_table, loc_table, lookup_exception;
    std::vector<std::vector<IndexType>> occ_table;
	std::vector<bool> is_sampled_table;

	explicit FMIndex(
		const std::array<IndexType, 256>& c_to_o, 
		const std::array<IndexType, CharTypeNum>& o_to_c
	)
		: seq_end_pos(0)
		, char_to_order_(c_to_o)
		, order_to_char_(o_to_c)
		, c_table(CharTypeNum, 0)
		, occ_table(CharTypeNum)
		, lookup_table(std::pow(CharTypeNum, LookupStrLen), (IndexType)-1)
	{
		if (PrefixLen != 0 && PrefixLen < LookupStrLen)
			throw std::runtime_error(
				"ERROR: FMIndex() fail, because PrefixLen is less "
				"than LookupStrLen. This will cause incorrect "
				"results.\n"
			);
	}

	FMIndex(const FMIndex&) = delete;

	FMIndex(FMIndex&&) = default;

	FMIndex& operator=(const FMIndex&) = delete;

	FMIndex& operator=(FMIndex&& fm)
	{
		this->char_to_order_ = std::move(fm.char_to_order_);
		this->order_to_char_ = std::move(fm.order_to_char_);
		this->bwt_seq = std::move(fm.bwt_seq);
		this->seq_end_pos = std::move(fm.seq_end_pos);
		this->c_table = std::move(fm.c_table);
		this->loc_table = std::move(fm.loc_table);
		this->lookup_table = std::move(fm.lookup_table);
		this->lookup_exception = std::move(fm.lookup_exception);
		this->occ_table = std::move(fm.occ_table);
		this->seg_info = std::move(fm.seg_info);
		this->n_table = std::move(fm.n_table);
	}

	template <
		typename IStream, 
		typename = std::enable_if_t<
			std::is_base_of_v<std::istream, std::decay_t<IStream>>
		>>
	void build(
		IStream&& is,
		Sorter& sorter)
	{
		SEQ seq;
		IndexType seg_pos(0);
        
		for (std::string buf; std::getline(is, buf);)
		{
			if (buf.front() != '>')
				throw std::runtime_error(
					"ERROR: Input file format error\n"
				);

			std::getline(is, buf);
			seg_info.emplace_back(seg_pos, buf.size());
			seg_pos += buf.size();
			make_seq(seq, buf);
		}

		build(seq, sorter);
	}

	/***********************************************************
	 * If users use this function to build FMIndex, they should 
	 * maintain n_table, seg_info by themselves.
	 **********************************************************/
	void build(
		const SEQ& seq, 
		Sorter& sorter)
	{
		//std::chrono::duration<double, std::milli> dur, dur2;
		std::chrono::duration<double, std::ratio<60>> dur, dur2;
		auto clock(std::chrono::steady_clock::now()), 
			clock2(std::chrono::steady_clock::now());

		bwt_seq.reserve(seq.size() + 1);
        for (auto& i : occ_table)
            i.reserve(
				bwt_seq.capacity() >> LogInterval << LogInterval
			);
        loc_table.reserve(
			bwt_seq.capacity() >> LogInterval << LogInterval
		);
		is_sampled_table.reserve(bwt_seq.capacity());

		handle_dollar_sign(seq);
		if constexpr (IsSBWT)
		{
			IndexType cumulative_idx(1);
			SuffixArrayType group;
			
			for (size_t i(0); 
				i < sorter.splitter.splitters.size() + 1; 
				i++)
			{
				group = sorter.sort_some(
						seq, char_to_order_, PrefixLen
					);
				dur += std::chrono::steady_clock::now() - clock;
				clock2 = std::chrono::steady_clock::now();
				build_impl(seq, group, cumulative_idx);
				dur2 += std::chrono::steady_clock::now() - clock2;
				clock = std::chrono::steady_clock::now();
			}
			std::cout << "sort: "  << dur.count() << " min\n";
		}
		else
		{
			IndexType cumulative_idx(1);
			SuffixArrayType group;

			sorter.sort(seq, group, char_to_order_, PrefixLen);
			build_impl(seq, group, cumulative_idx);
		}
		
		clock2 = std::chrono::steady_clock::now();
		c_table[0]++;

        for (auto it(c_table.begin()); it != c_table.end() - 1; it++)
            *(it + 1) += *it;
		for (auto rit(c_table.rbegin()); 
			rit != c_table.rend() - 1; 
			rit++)
			*rit = *(rit + 1);
		c_table.front() = 1;

		if (lookup_table.back() == (IndexType)-1)
			lookup_table.back() = bwt_seq.size();
		for (auto rit(lookup_table.rbegin()); 
			rit != lookup_table.rend() - 1; 
			rit++)
		{
			if (*(rit + 1) == (IndexType)-1)
				*(rit + 1) = *rit;
		}
		dur2 += std::chrono::steady_clock::now() - clock2;
		std::cout << "build: " << dur2.count() << " min\n";
	}

	template <typename IStream, typename OStream>
	void map(
		IStream&& is, 
		OStream&& os, 
		IndexType max_candidate_num = 0)
	{
		std::chrono::duration<double, std::milli> dur;
		auto clock(std::chrono::steady_clock::now());
		if (!is.is_open())
			throw std::runtime_error(
					"ERROR: istream passed into map() is not opened\n"
				);
		if constexpr (
			!std::is_same_v<
				std::decay_t<OStream>, 
				std::ostream
			>)
		{
			if (!os.is_open())
				throw std::runtime_error(
						"ERROR: ostream passed into map() "
						"is not opened\n"
					);
		}

		FASTQ fq;
		while (is >> fq)
		{
			std::vector<IndexType> seq_indice;
			std::pair<IndexType, std::pair<IndexType, IndexType>> 
				result;
			auto& range(result.second);

			if constexpr (IsSBWT)
			{
				result = sbwt_exact_match(fq.seq, max_candidate_num);

				if (result.first == -1)
				{
					seq_indice = sbwt_range_to_seq_idx(
							range, fq.seq.size()
						);
				}
			}
			else
			{
				result = exact_match(fq.seq, max_candidate_num);

				if (result.first == -1)
				{
					for (IndexType i(range.first); 
						i < range.second; 
						i++)
					{
						seq_indice.emplace_back(
							bwt_idx_to_seq_idx(i)
						);
					}
				}
			}
			os << range.first << '~' << range.second << ":\n";
			std::copy(seq_indice.cbegin(), seq_indice.cend(), 
				std::ostream_iterator<IndexType>(os, " "));
			os << '\n';
		}
		dur = std::chrono::steady_clock::now() - clock;
		std::cout << "map: " << dur.count() << "ms\n";
	}

	template <typename SeqType>
	auto exact_match(
		SeqType&& query, 
		IndexType max_candidate_num)
	{
		if (c_table.front() == 0)
			throw std::runtime_error(
				"ERROR: Should not do exact_match before building "
				"fm-index\n"
			);

		if (PrefixLen - (1 << LogInterval) < query.size())
			throw std::runtime_error(
				"ERROR: Length of query too large for "
				"PrefixLen of FMIndex"
			);
		
		std::pair<IndexType, std::pair<IndexType, IndexType>> ret;
		auto& range(ret.second);
		auto rit(query.crbegin());

		if (query.size() >= LookupStrLen)
		{
			lookup(query, range);
			if (range.first < range.second)
				rit += LookupStrLen;
      else
      {
        range.first = 0;
        range.second = bwt_seq.size();
      }
		}
		else
		{
			range.first = 0;
			range.second = bwt_seq.size();
		}
		

		for (auto rend(query.crend()); rit != rend; rit++)
		{
			range.first = lf_mapping(range.first, *rit);
			range.second = lf_mapping(range.second, *rit);
			if (range.first >= range.second)
			{
				ret.first = std::distance(query.crbegin(), rit);
				return ret;
			}
		}

		if (max_candidate_num != 0 && 
			range.first + max_candidate_num < range.second)
		{
			ret.first = -2;
			return ret;
		}
		
		ret.first = -1;
		return ret;
	}

	template <typename SeqType>
	auto exact_match(
		SeqType&& query, 
		IndexType max_candidate_num, 
		std::vector<std::pair<IndexType, IndexType>>& range_records) const
	{
		if (c_table.front() == 0)
			throw std::runtime_error(
				"ERROR: Should not do exact_match before building "
				"fm-index\n"
			);

		if (PrefixLen - (1 << LogInterval) < query.size())
			throw std::runtime_error(
				"ERROR: Length of query too large for "
				"PrefixLen of FMIndex"
			);
		
		std::pair<IndexType, std::pair<IndexType, IndexType>> ret;
		auto& range(ret.second);
		auto rit(query.crbegin());

		range_records.reserve(query.size());

		if (query.size() >= LookupStrLen)
		{
			lookup(query, range);
			if (range.first < range.second)
			{
				rit += LookupStrLen;
				range_records.emplace_back(range);
			}
		  else
      {
        range.first = 0;
        range.second = bwt_seq.size();
      }
    }
		else
		{
			range.first = 0;
			range.second = bwt_seq.size();
		}
		

		for (auto rend(query.crend()); rit != rend; rit++)
		{
			range.first = lf_mapping(range.first, *rit);
			range.second = lf_mapping(range.second, *rit);
			if (range.first >= range.second)
			{
				ret.first = std::distance(query.crbegin(), rit);
				range_records.shrink_to_fit();
				return ret;
			}
			range_records.emplace_back(range);
		}

		if (max_candidate_num != 0 && 
			range.first + max_candidate_num < range.second)
		{
			ret.first = -2;
			range_records.shrink_to_fit();
			return ret;
		}
		
		ret.first = -1;
		range_records.shrink_to_fit();
		return ret;
	}
	
  template <typename SeqType>
	auto exact_match_by_base(
		SeqType&& query, 
		IndexType max_candidate_num, 
		std::vector<std::pair<IndexType, IndexType>>& range_records) const
	{
		if (c_table.front() == 0)
			throw std::runtime_error(
				"ERROR: Should not do exact_match before building "
				"fm-index\n"
			);

		if (PrefixLen - (1 << LogInterval) < query.size())
			throw std::runtime_error(
				"ERROR: Length of query too large for "
				"PrefixLen of FMIndex"
			);
		
		std::pair<IndexType, std::pair<IndexType, IndexType>> ret;
		auto& range(ret.second);
		auto rit(query.crbegin());
    
		range_records.reserve(query.size());
    
    range.first = 0;
    range.second = bwt_seq.size();
		range_records.emplace_back(range);

		for (auto rend(query.crend()); rit != rend; rit++)
		{
			range.first = lf_mapping(range.first, *rit);
			range.second = lf_mapping(range.second, *rit);
			if (range.first >= range.second)
			{
				ret.first = std::distance(query.crbegin(), rit);
				range_records.shrink_to_fit();
				return ret;
			}
			range_records.emplace_back(range);
		}

		if (max_candidate_num != 0 && 
			range.first + max_candidate_num < range.second)
		{
			ret.first = -2;
			range_records.shrink_to_fit();
			return ret;
		}
		
		ret.first = -1;
		range_records.shrink_to_fit();
		return ret;
	}

	template <typename SeqType>
	auto sbwt_exact_match(
		SeqType&& query, 
		IndexType max_candidate_num)
	{
		if (c_table.front() == 0)
			throw std::runtime_error(
				"ERROR: Should not do sbwt_exact_match "
				"before building fm-index\n"
			);
		
		SEQ sbwt_query;

		make_rc_query(sbwt_query, query);
		
		return exact_match(sbwt_query, max_candidate_num);
	}

	template <typename SeqType>
	auto sbwt_exact_match(
		SeqType&& query, 
		IndexType max_candidate_num, 
		std::vector<std::pair<IndexType, IndexType>>& range_records) const
	{
		if (c_table.front() == 0)
			throw std::runtime_error(
				"ERROR: Should not do sbwt_exact_match "
				"before building fm-index\n"
			);
		
		SEQ sbwt_query;
		
		make_rc_query(sbwt_query, query);
		
		return exact_match(sbwt_query, max_candidate_num, 
				range_records
			);
	}
	
  template <typename SeqType>
	auto sbwt_exact_match_by_base(
		SeqType&& query, 
		IndexType max_candidate_num, 
		std::vector<std::pair<IndexType, IndexType>>& range_records) const
	{
		if (c_table.front() == 0)
			throw std::runtime_error(
				"ERROR: Should not do sbwt_exact_match "
				"before building fm-index\n"
			);
		
		SEQ sbwt_query;
		
		make_rc_query(sbwt_query, query);
		
		return exact_match_by_base(sbwt_query, max_candidate_num, 
				range_records
			);
	}

	inline IndexType lf_mapping(
		IndexType idx, 
		typename SEQ::value_type ch) const
	{
		return c_table[char_to_order_[ch]] + 
			occ_tracking(idx, ch, 
				idx >> LogInterval << LogInterval
			);
	}

	inline IndexType bwt_idx_to_seq_idx(IndexType idx)
	{
		IndexType count(0);
		for (; !is_sampled_table[idx]; count++)
			idx = lf_mapping(idx, bwt_seq[idx]);

		return std::lower_bound(
				loc_table.cbegin(), 
				loc_table.cend(), 
				idx, 
				[](const std::pair<IndexType, IndexType>& a, 
					IndexType b)
				{
					return a.first < b;
				}
			)->second + count;
	}

	std::vector<IndexType> sbwt_range_to_seq_idx(
		const std::pair<IndexType, IndexType>& range, 
		const IndexType& query_len) const
	{
		std::vector<IndexType> results;
		
		if (range.first >= range.second)
			return results;

		results.reserve(range.second - range.first);
		//results.resize(range.second - range.first);
		IndexType candidate_num = range.second - range.first;
    sbwt_range_to_seq_idx_impl(
			results, candidate_num, range, 0, query_len
		);
		return results;
	}

	void save(
		const std::string& filename, 
		bool is_binary_archive = true)
	{
		if (is_binary_archive)
		{
			std::ofstream ofs_b(
				filename, std::ios::out | std::ios::binary
			);

			if (!ofs_b.good())
				throw std::runtime_error(
					"ERROR: open file in save() fail\n"
				);
			
			boost::archive::binary_oarchive arch(ofs_b);
			arch << LogInterval << CharTypeNum << LookupStrLen << 
				PrefixLen << char_to_order_ << order_to_char_ << 
				bwt_seq << seq_end_pos << c_table << loc_table << 
				lookup_table << lookup_exception << occ_table << 
				seg_info << n_table << chr_names;
		}
		else
		{
			std::ofstream ofs_t(filename, std::ios::out);

			if (!ofs_t.good())
				throw std::runtime_error(
					"ERROR: open file in save() fail\n"
				);

			boost::archive::text_oarchive arch(ofs_t);
			arch << LogInterval << CharTypeNum << LookupStrLen << 
				PrefixLen << char_to_order_ << order_to_char_ << 
				bwt_seq << seq_end_pos << c_table << loc_table << 
				lookup_table << lookup_exception << occ_table << 
				seg_info << n_table << chr_names;
		}
	}

	void load(
		const std::string& filename, 
		bool is_binary_archive = true)
	{
		uint32_t log_interval, char_type_num, lookup_str_len, 
			prefix_len;
		if (is_binary_archive)
		{
			std::ifstream ifs_b(
				filename, std::ios::in | std::ios::binary
			);

			if (!ifs_b.good())
				throw std::runtime_error(
					"ERROR: open file in load() fail\n"
				);

			boost::archive::binary_iarchive arch(ifs_b);
			arch >> log_interval >> char_type_num >> 
				lookup_str_len >> prefix_len;

			if (log_interval != LogInterval || 
				char_type_num != CharTypeNum || 
				lookup_str_len != LookupStrLen || 
				prefix_len != PrefixLen)
			{
				throw std::runtime_error(
					"ERROR: Object traits are different from traits "
					"in file.\n"
				);
			}

			arch >> char_to_order_ >> order_to_char_ >> bwt_seq >> 
				seq_end_pos >> c_table >> loc_table >> 
				lookup_table >> lookup_exception >> occ_table >> 
				seg_info >> n_table >> chr_names;
		}
		else
		{
			std::ifstream ifs_t(filename, std::ios::in);

			if (!ifs_t.good())
				throw std::runtime_error(
					"ERROR: open file in load() fail\n"
				);

			boost::archive::text_iarchive arch(ifs_t);
			arch >> log_interval >> char_type_num >> 
				lookup_str_len >> prefix_len;

			if (log_interval != LogInterval || 
				char_type_num != CharTypeNum || 
				lookup_str_len != LookupStrLen || 
				prefix_len != PrefixLen)
			{
				throw std::runtime_error(
					"ERROR: Object traits are different from traits "
					"in file.\n"
				);
			}

			arch >> char_to_order_ >> order_to_char_ >> bwt_seq >> 
				seq_end_pos >> c_table >> loc_table >> 
				lookup_table >> lookup_exception >> occ_table >> 
				seg_info >> n_table >> chr_names;
		}
	}

  private:
	
	template <typename SeqType>
	void make_seq(SEQ& seq, SeqType&& buf)
	{
		seq.reserve(seq.capacity() + buf.size());
		for (auto it(buf.cbegin()); it != buf.cend(); it++)
		{
			if (*it == 'N' || *it == 'n')
			{
				n_table.emplace_back(
					std::distance(buf.cbegin(), it), 1
				);
				seq.push_back(
					order_to_char_[fast_rand() % CharTypeNum]
				);

				auto& n_count(n_table.back().second);
				for (it++; it != buf.cend(); it++)
				{
					switch (*it)
					{
						case 'N':
						case 'n':
							n_count++;
							seq.push_back(
								order_to_char_[
									fast_rand() % CharTypeNum
								]
							);
							break;
						default:
							it = buf.cend() - 1;
					}
				}
			}

			seq.push_back(*it);
		}
	}

	template <typename SeqType>
	void make_rc_query(SEQ& seq, SeqType&& buf) const
	{
		seq.reserve(buf.size());
		for (auto rit(buf.crbegin()); rit != buf.crend(); rit++)
		{
			switch(*rit)
			{
				case 'A':
					seq.push_back('T');
					break;
				case 'T':
					seq.push_back('A');
					break;
				case 'C':
					seq.push_back('G');
					break;
				case 'G':
					seq.push_back('C');
					break;
				default:
					throw std::runtime_error(
						"query contains invalid character"
					);
			}
		}
	}

	inline void handle_dollar_sign(const SEQ& seq)
	{
		for (IndexType i(0); i < c_table.size(); i++)
			occ_table[i].emplace_back(c_table[i]);
		
		if ((seq.size() & interval - 1) == 0)
			loc_table.emplace_back(0, seq.size());
		
		bwt_seq.push_back(seq.back());
		c_table[char_to_order_[bwt_seq.back()]]++;
		is_sampled_table.emplace_back(true);
	}
	
	inline void build_impl(
		const SEQ& seq, 
		const SuffixArrayType& group, 
		IndexType& cumulative_idx)
	{
		IndexType lookup_count(0);
		  /*constexpr size_t mask(0b11), 
			  bases_per_byte(SEQ::bases_per_word / sizeof(size_t)), 
			  log_bases_per_byte(std::log2(bases_per_byte)), 
			  bits_per_base(CHAR_BIT / bases_per_byte), 
			  log_bits_per_base(std::log2(bits_per_base));
      */
    size_t shift_base, reverse_lookup;

		for (auto it(group.cbegin());
			it < group.cend();
			it++, cumulative_idx++)
		{
			if ((cumulative_idx & interval - 1) == 0)
				for (IndexType i(0); i < c_table.size(); i++)
					occ_table[i].emplace_back(c_table[i]);

			switch (*it & interval - 1)
			{
				case 0:
					loc_table.emplace_back(cumulative_idx, *it);
					is_sampled_table.emplace_back(true);
					break;
				default:
					is_sampled_table.emplace_back(false);
			}
			
			if (*it + LookupStrLen <= seq.size())
			{
				lookup_count = 0;
				if constexpr (
					std::is_same_v<
						SEQ, 
						biovoltron::vector<typename SEQ::value_type>
					>)
				{
		      constexpr size_t mask(0b11), 
			      bases_per_byte(SEQ::bases_per_word / sizeof(size_t)), 
			      log_bases_per_byte(std::log2(bases_per_byte)), 
			      bits_per_base(CHAR_BIT / bases_per_byte), 
			      log_bits_per_base(std::log2(bits_per_base));
					
          shift_base = *it & SEQ::bases_per_word - 1; 
					reverse_lookup = 
						*(size_t*)(
							(uint8_t*)seq[*it].get_seg() + 
							(shift_base >> log_bases_per_byte)
						) >> 
						((shift_base & mask) << log_bits_per_base);

					for (size_t j(0); j < LookupStrLen; j++)
					{
						lookup_count |= 
							(reverse_lookup & mask) << 
							((LookupStrLen - 1 - j) << log_bits_per_base);
						reverse_lookup >>= 2;
					}
				}
				else
				{
					for (auto j(*it), end(*it + LookupStrLen); 
						j < end; 
						j++)
					{
						lookup_count += char_to_order_[seq[j]] *
							std::pow(CharTypeNum, end - j - 1);
					}
				}
				

				if (lookup_table[lookup_count] == (IndexType)-1)
					lookup_table[lookup_count] = cumulative_idx;
			}
			else
			{
				if (lookup_exception.size() == 0 || 
					lookup_exception.back().first != lookup_count)
				{
					lookup_exception.emplace_back(
						lookup_count, cumulative_idx
					);
				}
			}

			if (*it != 0)
			{
				bwt_seq.push_back(seq[*it - 1]);
				c_table[char_to_order_[bwt_seq.back()]]++;
			}
			else
			{
				bwt_seq.push_back(order_to_char_.front());
				seq_end_pos = cumulative_idx;
			}
		}
	}

	IndexType occ_tracking(
		IndexType idx, 
		typename SEQ::value_type ch, 
		IndexType checkpoint) const
	{
		IndexType occ;

		if (checkpoint == bwt_seq.size())
			checkpoint -= interval;

		if ((idx & interval - 1) >> (LogInterval - 1) == 0 || 
			checkpoint + interval >= bwt_seq.size())
		{
			occ = occ_table[char_to_order_[ch]]
				[checkpoint >> LogInterval];
			for (IndexType i(checkpoint); i < idx; i++)
				if (bwt_seq[i] == ch)
					occ++;
      if (ch == 'A' && checkpoint <= seq_end_pos && idx > seq_end_pos)
        occ--;
		}
		else
		{
			checkpoint += interval;
			occ = occ_table[char_to_order_[ch]]
				[checkpoint >> LogInterval];
			for (IndexType i(checkpoint - 1); i >= idx; i--)
			{
				if (bwt_seq[i] == ch)
					occ--;
			}
      if (ch == 'A' && checkpoint > seq_end_pos && idx <= seq_end_pos)
        occ++;
		}

		return occ;
	}

	template <typename SeqType>
	void lookup(
		SeqType&& query, 
		std::pair<IndexType, IndexType>& range) const
	{
		uint64_t lookup_count(0);

		if constexpr (
			std::is_same_v<
				SEQ, 
				biovoltron::vector<typename SEQ::value_type>
			>)
		{
			constexpr size_t mask(0b11), 
				bases_per_byte(SEQ::bases_per_word / sizeof(size_t)), 
				log_bases_per_byte(std::log2(bases_per_byte)), 
				log_bits_per_base(
					std::log2(CHAR_BIT / bases_per_byte)
				);
			size_t idx(query.size() - LookupStrLen), 
				shift_base(idx & SEQ::bases_per_word - 1), 
				reverse_lookup(
					*(size_t*)(
						(uint8_t*)query[idx].get_seg() + 
						(shift_base >> log_bases_per_byte)
					) >> 
					((shift_base & mask) << log_bits_per_base)
				);

			for (size_t j(0); j < LookupStrLen; j++)
			{
				lookup_count |= 
					(reverse_lookup & mask) << 
					((LookupStrLen - 1 - j) << log_bits_per_base);
				reverse_lookup >>= 2;
			}
		}
		else
		{
			for (auto i(query.size() - LookupStrLen); 
				i < query.size(); 
				i++)
			{
				lookup_count += char_to_order_[query[i]] *
					std::pow(CharTypeNum, query.size() - i - 1);
			}
		}
		
		auto it = std::lower_bound(
			lookup_exception.cbegin(), 
			lookup_exception.cend(), 
			lookup_count, 
			[](const std::pair<IndexType, IndexType>& a, 
				IndexType b)
			{
				return a.first < b;
			}
		);

		range.first = lookup_table[lookup_count];
		if (it != lookup_exception.cend() && 
			it->first == lookup_count
		)
			range.second = it->second;
		else
			range.second = lookup_table[lookup_count + 1];
	}

	void sbwt_range_to_seq_idx_impl(
		std::vector<IndexType>& results, 
		IndexType& candidate_num,
    const std::pair<IndexType, IndexType>& range, 
    IndexType depth, 
		const IndexType& query_len) const
	{
    //if (results.size() == results.capacity())
		if (candidate_num == 0)  
      return;
    if (depth == interval)
    {
        candidate_num = 0;
        return;
    }

		for (IndexType i(range.first); i < range.second; i++)
    {
			auto it(std::lower_bound(
				loc_table.cbegin(), 
				loc_table.cend(), 
				i, 
				[](const std::pair<IndexType, IndexType>& a, 
					IndexType b)
				{
					return a.first < b;
				}
			));

			if (it->first == i)
			{
				record_if_valid(
					results, it->second + depth, query_len, candidate_num 
				);
        if (candidate_num ==0)
            return;
			}
		}

		std::pair<IndexType, IndexType> range_copy;
		for (auto ch : order_to_char_)
		{
			range_copy.first = lf_mapping(range.first, ch);
			range_copy.second = lf_mapping(range.second, ch);

      IndexType num(range_copy.second - range_copy.first);
			if ((int32_t)num <= 0)
				continue;
			
			sbwt_range_to_seq_idx_impl(
				results, num, range_copy, depth + 1, query_len
			);
		}
	}

	inline void record_if_valid( 
		std::vector<IndexType>& results, 
		IndexType pos, 
		IndexType query_len,
    IndexType& candidate_num
) const
	{
		auto comp(
			[](IndexType a, 
				const std::pair<IndexType, IndexType>& b)
			{
				return a < b.first;
			});
    auto comp_lower(
       [](const std::pair<IndexType, IndexType>& a, 
         const IndexType& b)
      {
        return (a.first + a.second) <= b;
      });
		auto it(std::upper_bound(
				seg_info.cbegin(), seg_info.cend(), pos, comp
			));
		auto n_it(std::lower_bound(
        n_table.cbegin(), n_table.cend(), pos, comp_lower
      ));

		if ((it == std::upper_bound(
					seg_info.cbegin(), seg_info.cend(), 
					pos + query_len - 1, comp
				)
			) && 
			(std::upper_bound(
				n_table.cbegin(), n_table.cend(), pos, comp
			) == 
			std::upper_bound(
				n_table.cbegin(), n_table.cend(), 
				pos + query_len - 1, comp
			)) &&
      (n_it == n_table.cend() ||
       pos + query_len - 1 < (*n_it).first
      ))
		{
      //std::cout << "a";
      results.emplace_back(pos);
		}
    //std::cout << "b\n";
    candidate_num--; 
    
	}

	auto to_sam(
		const FASTQ& fq, 
		const std::vector<std::pair<IndexType, IndexType>>& results, 
		bool is_rc = false)
	{		
		std::vector<SamType> ret;
		ret.reserve(results.size());

		if (is_rc)
		{
			for (const auto& result : results)
			{
				ret.emplace_back(
					fq.name, 
					16, 
					"chr" + std::to_string(result.first), 
					result.second + 1, 
					255, 
					std::to_string(fq.seq.size()) + 'M', 
					"*", 
					0, 
					0, 
					fq.seq,
					"*", 
					results.size()
				);
			}
		}
		else
		{
			for (const auto& result : results)
			{
				ret.emplace_back(
					fq.name, 
					0, 
					"chr" + std::to_string(result.first), 
					result.second + 1, 
					255, 
					std::to_string(fq.seq.size()) + 'M', 
					"*", 
					0, 
					0, 
					fq.seq,
					"*", 
					results.size()
				);
			}
		}

		return ret;
	}

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
};

}
