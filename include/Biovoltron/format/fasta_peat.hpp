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
 * @class fastaException
 * @brief An exception class which inherit std::runtime_error
 *
 * This class will be thrown if actions which belongs to VCF have 
 * exceptions.
 */
class fastaException : public std::runtime_error
{
  public:
    /**
     * @brief Message assign constructor
     *
     * @param m Message string wanted to be assigned
     */
    fastaException (const std::string& msg)
	: runtime_error(std::string("fastaException " + msg))
    {
    }
};

/**
 * @class FASTA_PE
 * @brief A FASTA_PE object which stores a line of fasta format, that is 
 * an entry of fasta data.
 * 
 * Member variable:<br>
 *	State state_<br>
 *	std::string name<br>
 *	NTable n_base_info_table<br>
 *	Sequence seq<br>
 *
 * This class can store an entry of fasta data, and also provides 
 * entry parsing, access, modify and entry processing functions, 
 * like substr(), trim(), get_antisence().<br>
 * And also provides sequence compression by using 
 * biovoltron::base_vector, so that the size of sequence can be 
 * decressed to 1/4.
 */
template <typename Sequence = std::string>
class FASTA_PE
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
     * parsing fasta files.
     *
     * This enum contains name, seq, plus, qual, end. The last one 
     * "end" is used for distinguish the end of a round. And others 
     * indicate different categories of lines in fasta format.
     */
    enum class State {name, seq, end};

    ///An State variable which is initialized to State::name
    State state_ = State::name;

  public:
    ///A string which store the name line of a fasta entry
    std::string name;

    /**
     * A NTable which store the start point and length of all 
     * continuous n-base sub-sequence
     */
    NTable n_base_info_table;

    ///A Sequence which store sequence, Sequence is template parameter
    Sequence seq;

    /**
     * @brief >> operator which get entry from an istream
     * 
     * @param is An istream where to get entry
     * @param fa A FASTA_PE object to store gotten entry
     * 
     * @return Is identical to parameter is
     *
     * Just call FASTA_PE::get_obj() to help parse given istream, and 
     * store them in given FASTA_PE object.<br>
     * Should have the some property to get_obj function.
     *
     * @sa FASTA_PE::get_obj()
     */
    friend std::istream& operator>>(std::istream& is, FASTA_PE& fa)
    {
	    return get_obj(is, fa);
    }

    /**
     * @brief << operator which write entry to an ostream
     * 
     * @param os An ostream where to write entry
     * @param fa A FASTA_PE object which is going to write to os
     * 
     * @return Is identical to parameter os
     *
     * Just call FASTA_PE::to_string() to help change data of given FASTA_PE  
     * object to string then output to given ostream, and store them 
     * in given FASTA_PE object.<br>
     * Should have the some property to FASTA_PE::to_string()
     *
     * @sa FASTA_PE::to_string()
     */
    friend std::ostream& operator<<(std::ostream& os, const FASTA_PE& fa)
    {
	    os << fa.to_string();
	    return os;
    }

    /**
     * @brief Parse fasta entry from specific istream
     *
     * @param is An istream contains entry data
     * @param fa A FASTA_PE object used to store parsed entry data
     * @return An istream identical to parameter is
     *
     * Use std::getline to get lines from istream and store in fa. 
     * Because there are 4 types of line for each fasta entry, use 
     * for loop to iterate state_ and use switch to compare states. 
     * If the format is correct, then just store data to respective 
     * variable. Otherwise, throw fastaException with the reason to 
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
     * @sa FASTA_PE::fast_rand()
     */
    static std::istream& get_obj(std::istream& is, FASTA_PE& fa)
    {
	std::string buf;

	for (fa.state_ = State::name; 
	    fa.state_ != State::end; 
	    fa.state_ = (State)((size_t)fa.state_ + 1))
	{
	    switch (fa.state_)
	    {
		case State::name:
		    if (!is.good())
			return is;

		    std::getline(is, fa.name, '\n');

		    if (fa.name.size() > 0 && fa.name.front() == '>')
			fa.name.erase(fa.name.begin());
		    else
			throw fastaException(
			    "ERROR: get_obj(): Can't find \'>\' from "
			    "input when state is equal to name\n"
			);

		    break;
		case State::seq:
		    fa.seq.clear();

		    if (!is.good())
			throw fastaException(
			    "ERROR: get_obj(): seq field of an entry "
			    "is disappear\n"
			);

		    std::getline(is, buf, '\n');
		    fa.seq.reserve(buf.size());

		    for (size_t i(0); i < buf.size(); i++)
		    {
			switch (buf.at(i))
			{
			    case 'A':
			    case 'a':
				fa.seq.push_back('A');
				break;
			    case 'C':
			    case 'c':
				fa.seq.push_back('C');
				break;
			    case 'G':
			    case 'g':
				fa.seq.push_back('G');
				break;
			    case 'T':
			    case 't':
				fa.seq.push_back('T');
				break;
			    case 'N':
			    case 'n':
				fa.n_base_info_table.emplace_back(i, 0);
				for (; i < buf.size() &&
						(buf.at(i) == 'N' ||
						buf.at(i) == 'n'); 
					i++)
				{
				switch (fast_rand() % 4)
				{
					case 0:
					fa.seq.push_back('A');
					break;
					case 1:
					fa.seq.push_back('C');
					break;
					case 2:
					fa.seq.push_back('G');
					break;
					case 3:
					fa.seq.push_back('T');
					break;
				}
				fa.n_base_info_table.back().second++;
				}
				i--;

				break;
			    default:
				throw fastaException(
				    "ERROR: get_obj(): invalid input "
				    "seq charactor\n"
				);
			}
		    }
		    break;

		default:
		    break;
	    }
	}

	return is;
    }

    /**
     * @brief Parse fasta entry from specific container which stores 
	 * 4 lines of fasta entry in 4 strings continuous
     *
     * @param is An container iterator points to the first line(aka 
	 * name line of a fasta entry)
     * @return The fasta entry required
     *
     * This function copy 4 lines of string from iterator to 
	 * iterator + 3.<br>
     * Because there are 4 types of line for each fasta entry, use 
     * for loop to iterate state_ and use switch to compare states. 
     * If the format is correct, then just store data to respective 
     * variable. Otherwise, throw fastaException with the reason to 
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
     * @sa FASTA_PE::fast_rand()
     */
	template <typename Iterator>
	static FASTA_PE parse_obj(Iterator it)
	{
		FASTA_PE fa;
		for (fa.state_ = State::name; 
			fa.state_ != State::end; 
			fa.state_ = (State)((size_t)fa.state_ + 1), it++)
		{
			switch (fa.state_)
			{
				case State::name:
					fa.name = std::move(*it);

					if (fa.name.size() > 0 && fa.name.front() == '>')
						fa.name.erase(fa.name.begin());
					else
						throw fastaException(
						"ERROR: get_obj(): format of name field is "
						"invalid\n"
						);

					break;
				case State::seq:
					fa.seq.reserve(it->size());

					for (size_t i(0); i < it->size(); i++)
					{
						switch (it->at(i))
						{
						case 'A':
						case 'a':
							fa.seq.push_back('A');
							break;
						case 'C':
						case 'c':
							fa.seq.push_back('C');
							break;
						case 'G':
						case 'g':
							fa.seq.push_back('G');
							break;
						case 'T':
						case 't':
							fa.seq.push_back('T');
							break;
						case 'N':
						case 'n':
							fa.n_base_info_table.emplace_back(i, 0);
							for (; i < it->size() &&
									(it->at(i) == 'N' ||
									it->at(i) == 'n'); 
								i++)
							{
							switch (fast_rand() % 4)
							{
								case 0:
								fa.seq.push_back('A');
								break;
								case 1:
								fa.seq.push_back('C');
								break;
								case 2:
								fa.seq.push_back('G');
								break;
								case 3:
								fa.seq.push_back('T');
								break;
							}
							fa.n_base_info_table.back().second++;
							}
							i--;

							break;
						default:
							throw fastaException(
							"ERROR: get_obj(): invalid input "
							"seq charactor\n"
							);
						}
					}

					break;
			}
		}

		return fa;
    }

    /**
     * @brief Can get data of this FASTA_PE object
     *
     * @return A string of data in fasta format
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
        std::string buf(">");

        buf.reserve(seq.size() * 2);
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
                throw fastaException(
                "ERROR: to_string(): invalid character in seq\n"
                );
            }
        }

        for (auto& i : n_base_info_table)
            for (unsigned j(0); j < i.second; j++)
                buf.at(seq_pos + i.first + j) = 'N';

        return buf;
    }

    template <typename Iterator>
    void to_iterator(Iterator it)
    {
        (*it) = ">" + std::move(name);
    }

    /**
     * @brief Convert vector of FASTA_PE to string then dump to an 
     * ostream
     *
     * @param os An ostream to dump fasta format data
     * @param v_fasta A vector<FASTA_PE>, which we want to dump them to 
     * fasta format
     *
     * Convert vector<FASTA_PE> to string by using range-for to append 
     * each fasta string which generated by FASTA_PE::to_string() to a 
     * buffer. And write the buffer to os at the end.
     *
     * Time complexity: O(m * n)<br>
     *	    \e m: size of v_fasta
     *	    \e n: average length of seq of FASTA_PE objects
     *
     * @sa FASTA_PE::to_string()
     */
    static void dump(std::ostream& os, 
	    const std::vector<FASTA_PE>& v_fasta)
    {
	std::string buf("");

	for (const auto& i : v_fasta)
	{
	    buf.append(i.to_string());
	    buf.append("\n");
	}

	buf.pop_back();
	os << buf;
    }

    /**
     * @brief return a FASTA_PE object with the same name but have 
     * specified sub-sequence and sub-sequence_quality
     *
     * @param pos The start point of specified sub-sequence
     * @param count The length of specified sub-sequence, if (pos + 
     * count) is greater than size of sequence, count will change to 
     * (size of sequence - pos) automatically
     * @return A FASTA_PE object with the same name like original one, 
     * but has sub-sequence and sub-sequence_quality
     *
     * Generate a FASTA_PE object with the same name like origin, but
     * sub-sequence and sub-sequence_quality are the specific part
     * of origin.<br>
     * Also modify n_base_info_table if there are continueous 
     * n-base sub-sequence are cut off.
     *
     * Time complexity: O(n)<br>
     *	    \e n: size of n_base_info_table
     */
    FASTA_PE substr(size_t pos, size_t count = -1) const
    {
	if (pos >= seq.size())
	    throw fastaException(
		"ERROR: substr(): out_of_range_exception\n"
	    );

	FASTA_PE tmp;
	const auto& fa_n(n_base_info_table);
	size_t n_end, begin, end;

	count = std::min(count, seq.size() - pos);

	tmp.name = name;
	tmp.seq = {seq.begin() + pos, seq.begin() + pos + count};
	tmp.n_base_info_table.reserve(fa_n.size());

	for (size_t i(0); i < fa_n.size(); i++)
	{
	    n_end = fa_n.at(i).first + fa_n.at(i).second;

	    if (pos < fa_n.at(i).first)	begin = fa_n.at(i).first - pos;
	    else if (pos >= n_end)	continue;
	    else			begin = 0;

	    if (pos + count < fa_n.at(i).first)	break;
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
     * @sa FASTA_PE::substr()
     */
    void trim(size_t pos)
    {
	if (pos > seq.size())
	    throw fastaException(
		"ERROR: trim(): out_of_range_exception\n"
	    );

	auto& fa_n(n_base_info_table);
	size_t n_end;

	seq.resize(pos);

	for (size_t i(0); i < fa_n.size(); i++)
	{
	    n_end = fa_n.at(i).first + fa_n.at(i).second;
	    
	    if (n_end <= pos)
		continue;
	    else if (n_end > pos && fa_n.at(i).first < pos)
		fa_n.at(i).second = pos - fa_n.at(i).first;
	    else
	    {
		fa_n.resize(i);
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
    >::type get_antisence() const noexcept
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
    >::type get_antisence() const noexcept
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
	fast_rand_seed.at(0) ^= fast_rand_seed.at(0) << 16;
	fast_rand_seed.at(0) ^= fast_rand_seed.at(0) >> 5;
	fast_rand_seed.at(0) ^= fast_rand_seed.at(0) << 1;

	t = fast_rand_seed.at(0);
	fast_rand_seed.at(0) = fast_rand_seed.at(1);
	fast_rand_seed.at(1) = fast_rand_seed.at(2);
	fast_rand_seed.at(2) = t ^ fast_rand_seed.at(0) ^ fast_rand_seed.at(1);

	return fast_rand_seed.at(2);
    }
};

}
