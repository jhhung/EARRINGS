/**
 * @file newlineSeperatedFasta.hpp
 * @brief a parser of newlineSeperatedFasta format file, 
 * only support std::string seq type, and each reads has
 * only one line.
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

namespace biovoltron::format{

/**
 * @class FASTAException
 * @brief An exception class which inherit std::runtime_error
 *
 * This class will be thrown if actions which belongs to Fasta have 
 * exceptions.
 */
class FASTAException : public std::runtime_error
{
  public:
    /**
     * @brief Message assign constructor
     *
     * @param m Message string wanted to be assigned
     */
    FASTAException (const std::string& msg)
	: runtime_error(std::string("FASTAException " + msg))
    {
    }
};

/**
 * @brief A static vector which stores 3 fast random seeds for fast 
 * random function of FASTQ class
 */
static std::vector<unsigned long> 
    fast_rand_seed_a = {123456789, 362436069, 521288629};

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
template <typename Sequence = std::string>
class FASTA
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
    enum class State {name, seq, end};

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
    friend std::istream& operator>>(std::istream& is, FASTA& fa)
    {
      return get_obj(is, fa);
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
    friend std::ostream& operator<<(std::ostream& os, const FASTA& fa)
    {
      os << fa.to_string();
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
    static std::istream& get_obj(std::istream& is, FASTA& fa)
    {
      std::string buf;

      for (fa.state_ = State::name; 
          fa.state_ != State::end; 
          fa.state_ = (State)((size_t)fa.state_ + 1))
      {
          switch (fa.state_)
          {
            case State::name:
              std::getline(is, fa.name, '\n');
              if (!is.good())
                return is;

              if (*fa.name.begin() == '>')
                fa.name.erase(fa.name.begin());
              else
                throw FASTAException(
                  "ERROR: get_obj(): Can't find \'>\' from "
                  "input when state is equal to name\n"
                );

              break;
            case State::seq:
              fa.seq.clear();
              std::getline(is, buf, '\n');

              if (!is.good())
                throw FASTAException(
                  "ERROR: get_obj(): seq field of an entry "
                  "is disappear\n"
                );
              fa.seq.reserve(buf.size());

              for (size_t i(0); i < buf.size(); i++)
              {
                switch (buf[i])
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
                    fa.n_base_info_table.emplace_back(
                      i, 0
                    );
                    i--;
                    while (buf[i + 1] == 'N' ||
                      buf[i + 1] == 'n')
                    {
                      switch (fast_rand_a() % 4)
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
                      i++;
                    }

                    break;
                  default:
                    throw FASTAException(
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
      std::string buf(">");

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
              throw FASTAException(
                "ERROR: to_string(): invalid character in seq\n"
              );
          }
      }

      for (auto& i : n_base_info_table)
        for (unsigned j(0); j < i.second; j++)
          buf[seq_pos + i.first + j] = 'N';

      buf.append("\n");
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
	    const std::vector<FASTA>& v_fasta)
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
    FASTA substr(size_t pos, size_t count = -1) const
    {
	    if (pos >= seq.size())
	      throw FASTAException(
		      "ERROR: substr(): out_of_range_exception\n"
	      );

	    FASTA tmp;
	    const auto& fa_n(n_base_info_table);
	    size_t n_end, begin, end;

	    count = std::min(count, seq.size() - pos);

	    tmp.name = name;
	    tmp.seq = {seq.begin() + pos, seq.begin() + pos + count};
	    tmp.n_base_info_table.reserve(fa_n.size());

	    for (size_t i(0); i < fa_n.size(); i++)
	    {
        n_end = fa_n[i].first + fa_n[i].second;

        if (pos < fa_n[i].first)	begin = fa_n[i].first - pos;
        else if (pos >= n_end)	continue;
        else			begin = 0;

        if (pos + count < fa_n[i].first)	break;
        else if (pos + count >= n_end)	end = n_end - pos;
        else				end = count;

        if (begin < end)
          tmp.n_base_info_table.emplace_back(begin, end - begin);
	    }

	    return tmp;
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
	    !std::is_same<T, typename biovoltron::Sequence>::value
    , Sequence
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
    static unsigned long fast_rand_a()
    {
      unsigned long t;
      fast_rand_seed_a[0] ^= fast_rand_seed_a[0] << 16;
      fast_rand_seed_a[0] ^= fast_rand_seed_a[0] >> 5;
      fast_rand_seed_a[0] ^= fast_rand_seed_a[0] << 1;

      t = fast_rand_seed_a[0];
      fast_rand_seed_a[0] = fast_rand_seed_a[1];
      fast_rand_seed_a[1] = fast_rand_seed_a[2];
      fast_rand_seed_a[2] = t ^ fast_rand_seed_a[0] ^ fast_rand_seed_a[1];

      return fast_rand_seed_a[2];
    }
};

}
