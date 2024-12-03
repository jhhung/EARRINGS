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

using namespace biovoltron;

namespace EARRINGS {
  class FastaException : public std::runtime_error {
  public:
    FastaException(const std::string &msg)
            : runtime_error(std::string("FastaException " + msg)) {}
  };

  enum class FastaState { name, seq, end };

  FastaState& operator++(FastaState& state) {
    return state = static_cast<FastaState>(static_cast<size_t>(state) + 1);
  }

  template <bool Encoded = false, typename Sequence = void>
  class Fasta : public biovoltron::FastaRecord<Encoded> {
  protected:
    using NTable = std::vector<std::pair<uint64_t, uint64_t>>;
    FastaState state_ = FastaState::name;

    inline static char normalize_base(char base, const std::string& function_name) {
      if (!Codec::is_valid(base)) {
        throw FastaException("Error: " + function_name + "(): invalid input seq character\n");
      }
      return Codec::to_char(Codec::to_int(base));
    }

    static int fast_rand() noexcept {
      static std::mt19937 rng(std::random_device{}());
      static std::uniform_int_distribution<int> dist(0, 3);
      return dist(rng);
    }

  public:
    using biovoltron::FastaRecord<Encoded>::START_SYMBOL;
    using biovoltron::FastaRecord<Encoded>::name;
    std::conditional_t<std::is_void_v<Sequence>, decltype(biovoltron::FastaRecord<Encoded>::seq), Sequence> seq;
    NTable n_base_info_table;

    friend std::istream& operator>>(std::istream& is, Fasta& fa) {
      return get_obj(is, fa);
    }

    friend std::ostream& operator<<(std::ostream& os, const Fasta& fa) {
      os << fa.to_string();
      return os;
    }

    static std::istream& get_obj(std::istream& is, Fasta& fa) {
      std::string buf;

      for (fa.state_ = FastaState::name; fa.state_ != FastaState::end; ++fa.state_) {
        switch (fa.state_) {
          case FastaState::name:
            if (!is.good())
              return is;

            std::getline(is, fa.name, '\n');

            if (fa.name.size() > 0 && fa.name.front() == Fasta::START_SYMBOL)
              fa.name.erase(fa.name.begin());
            else
              throw FastaException(std::string("ERROR: get_obj(): Can't find '") +
                                   Fasta::START_SYMBOL + "' from input when state is equal to name\n");

            break;
          case FastaState::seq:
            fa.seq.clear();

            if (!is.good())
              throw FastaException("ERROR: get_obj(): seq field of an entry is disappear\n");

            std::getline(is, buf, '\n');
            fa.seq.reserve(buf.size());

            for (size_t i = 0; i < buf.size(); ++i) {
              char base = Fasta::normalize_base(buf.at(i), __func__);
              if (base == 'N') {
                fa.n_base_info_table.emplace_back(i, 0);
                for (; i < buf.size() && (buf.at(i) == 'N' || buf.at(i) == 'n'); ++i) {
                  fa.seq.push_back(Codec::to_char(Fasta::fast_rand()));
                  fa.n_base_info_table.back().second++;
                }
                i--;
              } else {
                fa.seq.push_back(base);
              }
            }
            break;
          default:
            break;
        }
      }
      return is;
    }

    template <typename Iterator>
    static Fasta parse_obj(Iterator it) {
      Fasta fa;

      for (fa.state_ = FastaState::name; fa.state_ != FastaState::end; ++fa.state_, ++it) {
        switch (fa.state_) {
          case FastaState::name:
            fa.name = std::move(*it);

            if (fa.name.size() > 0 && fa.name.front() == Fasta::START_SYMBOL)
              fa.name.erase(fa.name.begin());
            else
              throw FastaException("ERROR: parse_obj(): format of name field is invalid\n");
            break;
          case FastaState::seq:
            fa.seq.reserve(it->size());

            for (size_t i = 0; i < it->size(); ++i) {
              char base = Fasta::normalize_base(it->at(i), __func__);
              if (base == 'N') {
                fa.n_base_info_table.emplace_back(i, 0);
                for (; i < it->size() && (it->at(i) == 'N' || it->at(i) == 'n'); ++i) {
                  fa.seq.push_back(Codec::to_char(Fasta::fast_rand()));
                  fa.n_base_info_table.back().second++;
                }
                i--;
              } else {
                fa.seq.push_back(base);
              }
            }
            break;
          default:
            break;
        }
      }
      return fa;
    }

    std::string to_string() const noexcept {
      std::string buf(Fasta::START_SYMBOL);

      buf.reserve(seq.size() * 2);
      buf.append(name);
      buf.append("\n");

      size_t seq_pos = buf.size();

      for (const auto i : seq) {
        if (Codec::is_valid(i)) {
          buf.push_back(i);
        } else {
          throw FastaException("ERROR: to_string(): invalid character in seq\n");
        }
      }

      for (auto& i : n_base_info_table)
        for (unsigned j = 0; j < i.second; ++j)
          buf.at(seq_pos + i.first + j) = 'N';

      return buf;
    }

    template <typename Iterator>
    void to_iterator(Iterator it) const {
      *it = Fasta::START_SYMBOL + name;
    }

    static void dump(std::ostream& os, const std::vector<Fasta>& v_fasta) {
      std::string buf;

      for (const auto& i : v_fasta) {
        buf.append(i.to_string());
        buf.append("\n");
      }

      if (!buf.empty())
        buf.pop_back();

      os << buf;
    }

    Fasta substr(size_t pos, size_t count = std::string::npos) const {
      if (pos >= seq.size())
        throw FastaException("ERROR: substr(): out_of_range_exception\n");

      Fasta tmp;
      const auto& fa_n = n_base_info_table;
      size_t n_end, begin, end;

      count = std::min(count, seq.size() - pos);

      tmp.name = name;
      tmp.seq = {seq.begin() + pos, seq.begin() + pos + count};
      tmp.n_base_info_table.reserve(fa_n.size());

      for (size_t i = 0; i < fa_n.size(); ++i) {
        n_end = fa_n.at(i).first + fa_n.at(i).second;

        if (pos < fa_n.at(i).first) begin = fa_n.at(i).first - pos;
        else if (pos >= n_end) continue;
        else begin = 0;

        if (pos + count < fa_n.at(i).first) break;
        else if (pos + count >= n_end) end = n_end - pos;
        else end = count;

        if (begin < end)
          tmp.n_base_info_table.emplace_back(begin, end - begin);
      }

      return tmp;
    }

    void trim(size_t pos) {
      if (pos > seq.size())
        throw FastaException("ERROR: trim(): out_of_range_exception\n");

      auto& fa_n = n_base_info_table;
      size_t n_end;

      seq.resize(pos);

      for (size_t i = 0; i < fa_n.size(); ++i) {
        n_end = fa_n.at(i).first + fa_n.at(i).second;

        if (n_end <= pos)
          continue;
        else if (n_end > pos && fa_n.at(i).first < pos)
          fa_n.at(i).second = pos - fa_n.at(i).first;
        else {
          fa_n.resize(i);
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
