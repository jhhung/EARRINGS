/**
 * @file vcf.hpp
 * @brief a parser of vcf format file
 *
 * @author JHH corp
 */

#pragma once

#include <exception>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <Biovoltron/base_vector.hpp>

namespace biovoltron::format{

/**
 * @class FASTQException
 * @brief An exception class which inherit std::runtime_error
 *
 * This class will be thrown if actions which belongs to VCF have 
 * exceptions.
 */
class FASTQException : public std::runtime_error
{
  public:
    /**
     * @brief Message assign constructor
     *
     * @param m Message string wanted to be assigned
     */
    FASTQException (const std::string& msg)
	: runtime_error(std::string("FASTQException " + msg))
    {
    }
};

/**
 * @brief A static vector which stores 3 fast random seeds for fast 
 * random function of FASTQ class
 */
static std::vector<unsigned long> 
    fast_rand_seed = {123456789, 362436069, 521288629};

/**
 * @class FASTQ
 * @brief A FASTQ object which stores a line of fastq format, that is 
 * an entry of fastq data.
 * 
 * Member variable:<br>
 *	State state_<br>
 *	std::string name<br>
 *	NTable n_base_info_table<br>
 *	Sequence seq<br>
 *	std::string seq_qual<br>
 *
 * This class can store an entry of fastq data, and also provides 
 * entry parsing, access, modify and entry processing functions, 
 * like substr(), trim(), get_antisense().<br>
 * And also provides sequence compression by using 
 * biovoltron::base_vector, so that the size of sequence can be 
 * decressed to 1/4.
 *
 * @tparam Sequence Type to store sequences, std::string is defaulted 
 * . If need compression, can use biovoltron::base_vector instead.
 * @tparam QualType A type which can store the quality line of a 
 * fastq entry, std::string is default. If user want to use their 
 * type, that type should work like std::string, that means it 
 * should at least support "clear()", "reserve()", "size()", 
 * "push_back()", "append()", "substr()", "resize()" "and operator[]"
 * . For example, biovoltron::vector is a proper type.
 *
 */
template <typename Sequence = std::string, typename QualType = std::string>
class FASTQ
{
    /**
     * @brief Simulate a table by using 
     * vector<pair<long int, long int>>
     *
     * NTable stores start point and length of every continueous 
     * n-base sequences
     */
    using NTable = std::vector<std::pair<uint64_t, uint64_t>>;

    /**
     * @enum state
     * @brief An enum class which indicates different state during 
     * parsing fastq files.
     *
     * This enum contains name, seq, plus, qual, end. The last one 
     * "end" is used for distinguish the end of a round. And others 
     * indicate different categories of lines in fastq format.
     */
    enum class State {name, seq, plus, qual, end};

    ///An State variable which is initialized to State::name
    State state_ = State::name;

  public:
    
    using SEQ = Sequence;
    ///A string which store the name line of a fastq entry
    std::string name;

    /**
     * A NTable which store the start point and length of all 
     * continuous n-base sub-sequence
     */
    NTable n_base_info_table;

    ///A Sequence which store sequence, Sequence is template parameter
    Sequence seq;

    ///A QualType variable which store the quality line of a fastq entry
    QualType seq_qual;

    /**
     * @brief >> operator which get entry from an istream
     * 
     * @param is An istream where to get entry
     * @param fq A FASTQ object to store gotten entry
     * 
     * @return Is identical to parameter is
     *
     * Just call FASTQ::get_obj() to help parse given istream, and 
     * store them in given FASTQ object.<br>
     * Should have the some property to get_obj function.
     *
     * @sa FASTQ::get_obj()
     */
    friend std::istream& operator>>(std::istream& is, FASTQ& fq)
    {
	return get_obj(is, fq);
    }

    /**
     * @brief << operator which write entry to an ostream
     * 
     * @param os An ostream where to write entry
     * @param fq A FASTQ object which is going to write to os
     * 
     * @return Is identical to parameter os
     *
     * Just call FASTQ::to_string() to help change data of given FASTQ  
     * object to string then output to given ostream, and store them 
     * in given FASTQ object.<br>
     * Should have the some property to FASTQ::to_string()
     *
     * @sa FASTQ::to_string()
     */
    friend std::ostream& operator<<(std::ostream& os, const FASTQ& fq)
    {
	os << fq.to_string();
	return os;
    }

    /**
     * @brief Parse fastq entry from specific istream
     *
     * @param is An istream contains entry data
     * @param fq A FASTQ object used to store parsed entry data
     * @return An istream identical to parameter is
     *
     * Use std::getline to get lines from istream and store in fq. 
     * Because there are 4 types of line for each fastq entry, use 
     * for loop to iterate state_ and use switch to compare states. 
     * If the format is correct, then just store data to respective 
     * variable. Otherwise, throw FASTQException with the reason to 
     * warn user.<br>
     * Notice that, if there are n-bases in sequence. we store start 
     * point and length of every continueous n-base sub-sequence in 
     * n_base_info_table, and use fast_rand to generate a basic base 
     * used to replace a n-base. To do this, because we want to use
     * 2 bits to compress the sequence, so we need to use a basic
     * base to replace n-base.
     *
     * Time complexity: O(n)<br>
     *	    \e n: the length of sequence<br>
     *
     * @sa FASTQ::fast_rand()
     */
    static std::istream& get_obj(std::istream& is, FASTQ& fq)
    {
	std::string buf;

	for (fq.state_ = State::name; 
	    fq.state_ != State::end; 
	    fq.state_ = (State)((size_t)fq.state_ + 1))
	{
	    switch (fq.state_)
	    {
		case State::name:
		    std::getline(is, fq.name, '\n');
		    if (!is.good())
			return is;

		    if (*fq.name.begin() == '@')
			fq.name.erase(fq.name.begin());
		    else
			throw FASTQException(
			    "ERROR: get_obj(): Can't find \'@\' from "
			    "input when state is equal to name\n"
			);

		    break;
		case State::seq:
		    fq.seq.clear();
		    std::getline(is, buf, '\n');

		    if (!is.good())
			throw FASTQException(
			    "ERROR: get_obj(): seq field of an entry "
			    "is disappear\n"
			);
		    fq.seq.reserve(buf.size());

		    for (size_t i(0); i < buf.size(); i++)
		    {
			switch (buf[i])
			{
			    case 'A':
			    case 'a':
				fq.seq.push_back('A');
				break;
			    case 'C':
			    case 'c':
				fq.seq.push_back('C');
				break;
			    case 'G':
			    case 'g':
				fq.seq.push_back('G');
				break;
			    case 'T':
			    case 't':
				fq.seq.push_back('T');
				break;
			    case 'N':
			    case 'n':
				fq.n_base_info_table.emplace_back(
				    i, 0
				);
				i--;
				while (buf[i + 1] == 'N' ||
					buf[i + 1] == 'n')
				{
				    switch (fast_rand() % 4)
				    {
					case 0:
					    fq.seq.push_back('A');
					    break;
					case 1:
					    fq.seq.push_back('C');
					    break;
					case 2:
					    fq.seq.push_back('G');
					    break;
					case 3:
					    fq.seq.push_back('T');
					    break;
				    }
				    fq.n_base_info_table.back().second++;
				    i++;
				}

				break;
			    default:
				throw FASTQException(
				    "ERROR: get_obj(): invalid input "
				    "seq charactor\n"
				);
			}
		    }

		    break;
		case State::plus:
		    std::getline(is, buf, '\n');

		    if (buf[0] != '+')
			throw FASTQException(
			    "ERROR: get_obj(): There is not \'+\' "
			    "after seq line\n"
			);
		    else
			if (buf.size() > 1 && fq.name != buf.substr(1))
			    throw FASTQException(
				"ERROR: get_obj(): There is string"
				"after \'+\' but not equal to string "
				"after \'@\'\n"
			    );

		    break;
		case State::qual:
			fq.seq_qual.clear();

		    if (!is.good())
			throw FASTQException(
			    "ERROR: get_obj(): seq_qual of an entry "
			    "is disappear\n"
			);

		    if constexpr (std::is_same_v<
			    QualType, 
			    biovoltron::vector<biovoltron::char_type>
			>)
		    {
				std::getline(is, buf);
				fq.seq_qual.reserve(buf.size());

				for (auto i : buf)
					fq.seq_qual.push_back(i);
		    }
		    else
				std::getline(is, fq.seq_qual);

		    for (unsigned i(0); i < fq.seq_qual.size(); i++)
				if (fq.seq_qual[i] > '~' || fq.seq_qual[i] < '!')
					throw FASTQException(
					"ERROR: get_obj(): wrong charactor in "
					"quality string\n"
					);
			
		    break;
		default:
		    break;
	    }
	}

	return is;
    }

    /**
     * @brief Parse fastq entry from specific container which stores 
	 * 4 lines of fastq entry in 4 strings continuous
     *
     * @param is An container iterator points to the first line(aka 
	 * name line of a fastq entry)
     * @return The fastq entry required
     *
     * This function copy 4 lines of string from iterator to 
	 * iterator + 3.<br>
     * Because there are 4 types of line for each fastq entry, use 
     * for loop to iterate state_ and use switch to compare states. 
     * If the format is correct, then just store data to respective 
     * variable. Otherwise, throw FASTQException with the reason to 
     * warn user.<br>
     * Notice that, if there are n-bases in sequence. we store start 
     * point and length of every continueous n-base sub-sequence in 
     * n_base_info_table, and use fast_rand to generate a basic base 
     * used to replace a n-base. To do this, because we want to use
     * 2 bits to compress the sequence, so we need to use a basic
     * base to replace n-base.
     *
     * Time complexity: O(n)<br>
     *	    \e n: the length of sequence<br>
     *
     * @sa FASTQ::fast_rand()
     */
	template <typename Iterator>
	static FASTQ parse_obj(Iterator it)
	{
		FASTQ fq;
		for (fq.state_ = State::name; 
			fq.state_ != State::end; 
			fq.state_ = (State)((size_t)fq.state_ + 1), it++)
		{
			switch (fq.state_)
			{
				case State::name:
					fq.name = std::move(*it);

					if (fq.name.size() > 0 && fq.name.front() == '@')
						fq.name.erase(fq.name.begin());
					else
						throw FASTQException(
						"ERROR: get_obj(): format of name field is "
						"invalid\n"
						);

					break;
				case State::seq:
					fq.seq.reserve(it->size());

					for (size_t i(0); i < it->size(); i++)
					{
						switch (it->at(i))
						{
						case 'A':
						case 'a':
							fq.seq.push_back('A');
							break;
						case 'C':
						case 'c':
							fq.seq.push_back('C');
							break;
						case 'G':
						case 'g':
							fq.seq.push_back('G');
							break;
						case 'T':
						case 't':
							fq.seq.push_back('T');
							break;
						case 'N':
						case 'n':
							fq.n_base_info_table.emplace_back(i, 0);
							for (; i < it->size() &&
									(it->at(i) == 'N' ||
									it->at(i) == 'n'); 
								i++)
							{
							switch (fast_rand() % 4)
							{
								case 0:
								fq.seq.push_back('A');
								break;
								case 1:
								fq.seq.push_back('C');
								break;
								case 2:
								fq.seq.push_back('G');
								break;
								case 3:
								fq.seq.push_back('T');
								break;
							}
							fq.n_base_info_table.back().second++;
							}
							i--;

							break;
						default:
							throw FASTQException(
							"ERROR: get_obj(): invalid input "
							"seq charactor\n"
							);
						}
					}

					break;
				case State::plus:
					if (it->size() == 0 || it->at(0) != '+')
						throw FASTQException(
						"ERROR: get_obj(): There is not \'+\' "
						"after seq line\n"
						);
					else
						if (it->size() > 1 && fq.name != it->substr(1))
						throw FASTQException(
							"ERROR: get_obj(): There is string"
							"after \'+\' but not equal to string "
							"after \'@\'\n"
						);

					break;
				case State::qual:
					if constexpr (std::is_same_v<
							QualType, 
							biovoltron::vector<biovoltron::char_type>
						>)
					{
						fq.seq_qual.reserve(it->size());
						fq.seq_qual.assign(it->cbegin(), it->cend());
					}
					else
						fq.seq_qual = std::move(*it);

					if (fq.seq_qual.size() != fq.seq.size())
						throw FASTQException(
						"ERROR: get_obj(): length of seq_qual field "
						"is different to length of seq field\n"
						);


					for (auto i : fq.seq_qual)
						if (i > '~' || i < '!')
						{
							throw FASTQException(
							"ERROR: get_obj(): wrong charactor in "
							"quality string\n"
							);
						}
			}
		}

		return fq;
    }



    /**
     * @brief Can get data of this FASTQ object
     *
     * @return A string of data in fastq format
     *
     * Concatenate string of 4 lines by order, and if there are data 
     * in n_base_info_table, change those sub-sequence back to 
     * n-bases.
     *
     * Time Complexity: O(n)<br>
     *	    \e n is the length of sequence
     */
    std::string to_string() const
    {
	std::string buf("@");

	buf.reserve(seq.size() * 3);
	buf.append(name);
	buf.append("\n");

	size_t seq_pos(buf.size());

	for (const auto i : seq)
	{
	    switch (i)
	    {
		case 'A':
		    buf.append("A");
		    break;
		case 'C':
		    buf.append("C");
		    break;
		case 'G':
		    buf.append("G");
		    break;
		case 'T':
		    buf.append("T");
		    break;
		default:
		    throw FASTQException(
			"ERROR: to_string(): invalid character in seq\n"
		    );
	    }
	}

	for (auto& i : n_base_info_table)
	    for (unsigned j(0); j < i.second; j++)
		buf[seq_pos + i.first + j] = 'N';

	buf.append("\n+\n");
	buf.append(seq_qual);

	return buf;
    }

    /**
     * @brief Convert vector of FASTQ to string then dump to an 
     * ostream
     *
     * @param os An ostream to dump fastq format data
     * @param v_fastq A vector<FASTQ>, which we want to dump them to 
     * fastq format
     *
     * Convert vector<FASTQ> to string by using range-for to append 
     * each fastq string which generated by FASTQ::to_string() to a 
     * buffer. And write the buffer to os at the end.
     *
     * Time complexity: O(m * n)<br>
     *	    \e m: size of v_fastq
     *	    \e n: average length of seq of FASTQ objects
     *
     * @sa FASTQ::to_string()
     */
    static void dump(std::ostream& os, 
	    const std::vector<FASTQ>& v_fastq)
    {
	std::string buf("");

	for (const auto& i : v_fastq)
	{
	    buf.append(i.to_string());
	    buf.append("\n");
	}

	buf.pop_back();
	os << buf;
    }

    /**
     * @brief return a FASTQ object with the same name but have 
     * specified sub-sequence and sub-sequence_quality
     *
     * @param pos The start point of specified sub-sequence
     * @param count The length of specified sub-sequence, if (pos + 
     * count) is greater than size of sequence, count will change to 
     * (size of sequence - pos) automatically
     * @return A FASTQ object with the same name like original one, 
     * but has sub-sequence and sub-sequence_quality
     *
     * Generate a FASTQ object with the same name like origin, but
     * sub-sequence and sub-sequence_quality are the specific part
     * of origin.<br>
     * Also modify n_base_info_table if there are continueous 
     * n-base sub-sequence are cut off.
     *
     * Time complexity: O(n)<br>
     *	    \e n: size of n_base_info_table
     */
    FASTQ substr(size_t pos, size_t count = -1) const
    {
	if (pos >= seq.size())
	    throw FASTQException(
		"ERROR: substr(): out_of_range_exception\n"
	    );

	FASTQ tmp;
	const auto& fq_n(n_base_info_table);
	size_t n_end, begin, end;

	count = std::min(count, seq.size() - pos);

	tmp.name = name;
	tmp.seq = {seq.begin() + pos, seq.begin() + pos + count};
	tmp.seq_qual = seq_qual.substr(pos, count);
	tmp.n_base_info_table.reserve(fq_n.size());

	for (size_t i(0); i < fq_n.size(); i++)
	{
	    n_end = fq_n[i].first + fq_n[i].second;

	    if (pos < fq_n[i].first)	begin = fq_n[i].first - pos;
	    else if (pos >= n_end)	continue;
	    else			begin = 0;

	    if (pos + count < fq_n[i].first)	break;
	    else if (pos + count >= n_end)	end = n_end - pos;
	    else				end = count;

	    if (begin < end)
		tmp.n_base_info_table.emplace_back(begin, end - begin);
	}

	return tmp;
    }

    /**
     * @brief cut off sequence and sequence_quality after pos
     *
     * @param pos The start point of the part you want to cut off
     *
     * Cut off tail-end of sequence and sequence_quality, only 
     * reserve those parts before pos. Also modified 
     * n_base_info_table to make sure its make sence.
     *
     * Time complexity: O(n)<br>
     *	    \e n: size of n_base_info_table
     *
     * @sa FASTQ::substr()
     */
    void trim(size_t pos)
    {
	if (pos > seq.size())
	    throw FASTQException(
		"ERROR: trim(): out_of_range_exception\n"
	    );

	auto& fq_n(n_base_info_table);
	size_t n_end;

	seq.resize(pos);
	seq_qual.resize(pos);

	for (size_t i(0); i < fq_n.size(); i++)
	{
	    n_end = fq_n[i].first + fq_n[i].second;
	    
	    if (n_end <= pos)
		continue;
	    else if (n_end > pos && fq_n[i].first < pos)
		fq_n[i].second = pos - fq_n[i].first;
	    else
	    {
		fq_n.resize(i);
		break;
	    }
	}
    }

    /**
     * @brief return reverse complement of sequence, only enable if 
     * template parameter "Sequence" is biovoltron::base_vector
     *
     * @return A sequence which is reverse complement of origin, and 
     * the type is biovoltron::base_vector
     *
     * Generate a biovoltron::base_vector object with r_iterator to 
     * get reverse sequence then use the function flip defined in 
     * biovoltron::base_vector to implement every base.
     *
     * Time Complexity: O(n)<br>
     *	    \e n: length of sequence
     *
     * @sa base_vector
     */
    template <typename T = Sequence>
    typename std::enable_if<
	std::is_same<
	    T, typename biovoltron::Sequence
	>::value, Sequence
    >::type get_antisense() const noexcept
    {
	Sequence s(seq.rbegin(), seq.rend());
	s.flip();

	return s;
    }

    /**
     * @brief return reverse complement of sequence, only enable if 
     * template parameter "Sequence" isn't biovoltron::base_vector
     *
     * @return A sequence which is reverse complement of origin, and 
     * the type is biovoltron::base_vector
     *
     * Generate a Sequence object with r_iterator to get reverse 
     * sequence then flip every base to implement every base.
     *
     * Time Complexity: O(n)<br>
     *	    \e n: length of sequence
     */
    template <typename T = Sequence>
    typename std::enable_if<
	!std::is_same<
	    T, typename biovoltron::Sequence
	>::value, Sequence
    >::type get_antisense() const noexcept
    {
	Sequence s(seq.rbegin(), seq.rend());

	for (auto& i : s)
	{
	    switch(i)
	    {
		case 'A':
		    i = 'T';
		    break;
		case 'C':
		    i = 'G';
		    break;
		case 'G':
		    i = 'C';
		    break;
		case 'T':
		    i = 'A';
		    break;
	    }
	}

	return s;
    }

    /**
     * @brief generate a random number quickly
     *
     * @return A random number
     *
     * Generate a random number by doing bit-wise operation to three
     * seeds
     */
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

    template <typename Iterator>
    void to_iterator(Iterator it)
    {
        (*it) = "@" + std::move(name);
        (*(it + 3)) = std::move(seq_qual);

        for( std::size_t i = (*(it + 1)).size(); i > seq.size(); --i )
        {
            (*(it + 1)).pop_back();
        }
    }
};

}
