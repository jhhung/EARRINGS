#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>
#include <random>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <OldBiovoltron/base_vector.hpp>
#include <EARRINGS/fasta.hpp>

using namespace biovoltron;

namespace EARRINGS {
  class FastqException : public std::runtime_error {
  public:
    FastqException(const std::string &msg)
            : runtime_error(std::string("FastqException " + msg)) {}
  };

  enum class FastqState { name, seq, plus, qual, end };

  FastqState& operator++(FastqState& state) {
    state = static_cast<FastqState>(static_cast<size_t>(state) + 1);
    return state;
  }

  template <bool Encoded = false, typename Sequence = void>
  class Fastq : public biovoltron::FastqRecord<Encoded> {
  protected:
    using NTable = std::vector<std::pair<uint64_t, uint64_t>>;
    FastqState state_ = FastqState::name;

    inline static char normalize_base(char base, const std::string& function_name) {
      if (!Codec::is_valid(base)) {
        throw FastqException("Error: " + function_name + "(): invalid input seq character\n");
      }
      return Codec::to_char(Codec::to_int(base));
    }

    static int fast_rand() noexcept {
      static std::mt19937 rng;
      std::uniform_int_distribution<int> dist(0, 3);
      return dist(rng);
    }

  public:
    using biovoltron::FastqRecord<Encoded>::START_SYMBOL;
    using biovoltron::FastqRecord<Encoded>::DELIM;
    using biovoltron::FastqRecord<Encoded>::name;
    std::conditional_t<std::is_void_v<Sequence>, decltype(biovoltron::FastqRecord<Encoded>::seq), Sequence> seq;
    using biovoltron::FastqRecord<Encoded>::qual;
    NTable n_base_info_table;

    Fastq(const EARRINGS::Fasta<Encoded, Sequence>& fa) {
      name = fa.name;
      if (name.front() == Fasta<>::START_SYMBOL) name[0] = START_SYMBOL;
      seq = fa.seq;
      qual = std::string(fa.seq.size(), 'I');
    }

    friend std::istream& operator>>(std::istream& is, Fastq& fq) {
      return get_obj(is, fq);
    }

    friend std::ostream& operator<<(std::ostream& os, const Fastq& fq) {
      os << fq.to_string();
      return os;
    }

    static std::istream& get_obj(std::istream& is, Fastq& fq) {
      std::string buf;

      for (fq.state_ = FastqState::name; fq.state_ != FastqState::end; ++fq.state_) {
        switch (fq.state_) {
          case FastqState::name:
            std::getline(is, fq.name, '\n');
            if (!is.good())
              return is;

            if (*fq.name.begin() == START_SYMBOL)
              fq.name.erase(fq.name.begin());
            else
              throw FastqException(std::string("ERROR: get_obj(): Can't find '") + START_SYMBOL +
                                   "' from input when state is equal to name\n");

            break;
          case FastqState::seq:
            fq.seq.clear();
            std::getline(is, buf, '\n');

            if (!is.good())
              throw FastqException("ERROR: get_obj(): seq field of an entry is disappear\n");

            fq.seq.reserve(buf.size());

            for (size_t i = 0; i < buf.size(); ++i) {
              char base = Fastq::normalize_base(buf.at(i), __func__);
              if (base == 'N') {
                fq.n_base_info_table.emplace_back(i, 0);
                for (; i < buf.size() && (buf.at(i) == 'N' || buf.at(i) == 'n'); ++i) {
                  fq.seq.push_back(Codec::to_char(Fastq::fast_rand()));
                  fq.n_base_info_table.back().second++;
                }
                i--;
              } else {
                fq.seq.push_back(base);
              }
            }
            break;
          case FastqState::plus:
            std::getline(is, buf, '\n');

            if (buf[0] != DELIM)
              throw FastqException(std::string("ERROR: get_obj(): There is not '") + DELIM + "' after seq line\n");
            else
            if (buf.size() > 1 && fq.name != buf.substr(1))
              throw FastqException(
                      std::string("ERROR: get_obj(): There is string after '")
                      + DELIM + "' but not equal to string after " + START_SYMBOL + "'\n"
              );

            break;
          case FastqState::qual:
            fq.qual.clear();

            if (!is.good())
              throw FastqException("ERROR: get_obj(): qual of an entry is disappear\n");

            /*if constexpr (std::is_same_v<QualType, biovoltron::vector<biovoltron::char_type>>)
            {
              std::getline(is, buf);
              fq.qual.reserve(buf.size());

              for (auto i : buf)
                fq.qual.push_back(i);
            }
            else*/
            std::getline(is, fq.qual);

            for (unsigned i(0); i < fq.qual.size(); ++i)
              if (fq.qual[i] > '~' || fq.qual[i] < '!')
                throw FastqException("ERROR: get_obj(): wrong charactor in quality string\n");

            break;
          default:
            break;
        }
      }

      return is;
    }

    template <typename Iterator>
    static Fastq parse_obj(Iterator it) {
      Fastq fq;
      for (fq.state_ = FastqState::name; fq.state_ != FastqState::end; ++fq.state_, ++it) {
        switch (fq.state_) {
          case FastqState::name:
            fq.name = std::move(*it);

            if (fq.name.size() > 0 && fq.name.front() == START_SYMBOL)
              fq.name.erase(fq.name.begin());
            else
              throw FastqException("ERROR: get_obj(): format of name field is invalid\n");

            break;
          case FastqState::seq:
            fq.seq.reserve(it->size());

            for (size_t i = 0; i < it->size(); ++i) {
              char base = Fastq::normalize_base(it->at(i), __func__);
              if (base == 'N') {
                fq.n_base_info_table.emplace_back(i, 0);
                for (; i < it->size() && (it->at(i) == 'N' || it->at(i) == 'n'); ++i) {
                  fq.seq.push_back(Codec::to_char(Fastq::fast_rand()));
                  fq.n_base_info_table.back().second++;
                }
                i--;
              } else {
                fq.seq.push_back(base);
              }
            }
            break;
          case FastqState::plus:
            if (it->size() == 0 || it->at(0) != DELIM)
              throw FastqException(
                      std::string("ERROR: get_obj(): There is not '") +
                      DELIM + "' after seq line\n"
              );
            else
            if (it->size() > 1 && fq.name != it->substr(1))
              throw FastqException(
                      std::string("ERROR: get_obj(): There is string after '") +
                      DELIM + "' but not equal to string after '" + START_SYMBOL + "'\n"
              );

            break;
          case FastqState::qual:
            /*if constexpr (std::is_same_v<
                          QualType,
                          biovoltron::vector<biovoltron::char_type>
                          >)
            {
              fq.qual.reserve(it->size());
              fq.qual.assign(it->cbegin(), it->cend());
            }
            else*/
            fq.qual = std::move(*it);

            if (fq.qual.size() != fq.seq.size())
              throw FastqException(
                      "ERROR: get_obj(): length of qual field "
                      "is different to length of seq field\n"
              );


            for (auto i : fq.qual)
              if (i > '~' || i < '!') {
                throw FastqException("ERROR: get_obj(): wrong charactor in quality string\n");
              }

            break;
          default:
            break;
        }
      }

      return fq;
    }

    std::string to_string() const
    {
      std::string buf(START_SYMBOL);

      buf.reserve(seq.size() * 3);
      buf.append(name);
      buf.append("\n");

      size_t seq_pos(buf.size());

      for (const auto i : seq) {
        if (Codec::is_valid(i)) {
          buf.push_back(i);
        } else {
          throw FastqException("ERROR: to_string(): invalid character in seq\n");
        }
      }

      for (auto& i : n_base_info_table)
        for (unsigned j(0); j < i.second; ++j)
          buf[seq_pos + i.first + j] = 'N';

      buf.append("\n+\n");
      buf.append(qual);

      return buf;
    }

    template <typename Iterator>
    void to_iterator(Iterator it) {
      (*it) = START_SYMBOL + std::move(name);
      (*(it + 3)) = std::move(qual);

      for( std::size_t i = (*(it + 1)).size(); i > seq.size(); --i) {
        (*(it + 1)).pop_back();
      }
    }

    static void dump(std::ostream& os, const std::vector<Fastq>& v_fastq) {
      std::string buf("");

      for (const auto& i : v_fastq) {
        buf.append(i.to_string());
        buf.append("\n");
      }

      buf.pop_back();
      os << buf;
    }

    Fastq substr(size_t pos, size_t count = -1) const
    {
      if (pos >= seq.size())
        throw FastqException("ERROR: substr(): out_of_range_exception\n");

      Fastq tmp;
      const auto& fq_n = n_base_info_table;
      size_t n_end, begin, end;

      count = std::min(count, seq.size() - pos);

      tmp.name = name;
      tmp.seq = {seq.begin() + pos, seq.begin() + pos + count};
      tmp.qual = qual.substr(pos, count);
      tmp.n_base_info_table.reserve(fq_n.size());

      for (size_t i = 0; i < fq_n.size(); ++i) {
        n_end = fq_n[i].first + fq_n[i].second;

        if (pos < fq_n[i].first) begin = fq_n[i].first - pos;
        else if (pos >= n_end) continue;
        else begin = 0;

        if (pos + count < fq_n[i].first) break;
        else if (pos + count >= n_end) end = n_end - pos;
        else end = count;

        if (begin < end)
          tmp.n_base_info_table.emplace_back(begin, end - begin);
      }

      return tmp;
    }

    void trim(size_t pos)
    {
      if (pos > seq.size())
        throw FastqException("ERROR: trim(): out_of_range_exception\n");

      auto& fq_n = n_base_info_table;
      size_t n_end;

      seq.resize(pos);
      qual.resize(pos);

      for (size_t i = 0; i < fq_n.size(); ++i) {
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

    template <typename T = Sequence>
    T get_antisense() const noexcept
    {
      T s(seq.rbegin(), seq.rend());

      if constexpr (std::is_same<T, biovoltron::Sequence>::value) {
        s.flip();
      } else {
        for (auto& i : s) {
          switch (i) {
            case 'A': i = 'T'; break;
            case 'C': i = 'G'; break;
            case 'G': i = 'C'; break;
            case 'T': i = 'A'; break;
          }
        }
      }

      return s;
    }
  };
}
