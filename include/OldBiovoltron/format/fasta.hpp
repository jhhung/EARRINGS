/**
 * @file fasta.hpp
 * @brief a constructor and parser of fasta file
 *
 * @author JHH corp
 */

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <type_traits>

namespace biovoltron
{
namespace format
{
    /**
     * @class Fasta
     * @brief A class which stores only one chromosome inside, and class SeqVec is 
     * nested inside
     *
     * Member variable:<br>
     * size_t chunk_size_<br>
     * SeqVec seq_<br>
     * std::string name_<br>
     *
     * This class store one chromosome, which can be compressed. And also 
     * provide function of parsing and storing.
     */
    template <class seq_type>
    class Fasta
    {
      private:
        ///the size of each element in vector
        size_t chunk_size_;
      public:
        /**
         * @class SeqVec
         * @brief A class that support [], which is ambigous in std::vector
         *
         * Member varible:<br>
         * private:
         * std::vector<seq_type> sv
         * size_t chunk_size__
         * 
         * This class is an alternative way to support [], which can find the
         * nucleotide in chrmosome based on the index. std::vector<seq_type> is 
         * a two dimension vector, so it is impossible to find the nucletoide by
         * index.
         */
        class SeqVec
        {
          private:
            /// sequence with chunk size, 
            std::vector<seq_type> sv;

            /// chunk size, which is power of 2
            size_t chunk_size__;
          public:

            /// default constructor
            SeqVec ()
                :sv(std::vector<seq_type>(0))
                ,chunk_size__(4096)
            {}

            /** @brief constructor.
             *
             * input chunk_size that will pass by 
             * constructor in Fasta. All constructor should pass in
             * the chunk_size
             */
            SeqVec (size_t chunk_size)
                :sv(std::vector<seq_type>(0))
                ,chunk_size__(chunk_size)
            {}

            /** @brief constructor.
             *
             * Input std::vector<seq_type>, which is the 
             * same type as member variable of SeqVec
             */
            SeqVec (const std::vector<seq_type>& v)
                :sv(v)
                ,chunk_size__(v[0].size())
            {}
            
            /// inititialize with std::string. It will divide the string
            /// into std::vector<seq_type>
            SeqVec (const std::string& s, size_t chunk_size)
                : chunk_size__(chunk_size)
            {
                set_vec(s);
            }

            /// copy constructor
            SeqVec (const SeqVec& s)
                :sv(s.sv)
                ,chunk_size__(s.chunk_size__)
            {}

            /// copy assignment
            SeqVec& operator=(const SeqVec& s)
            {
                sv = s.sv;
                chunk_size__ = s.chunk_size__;
                return *this;
            }

            /// move constructor
            SeqVec (std::vector<seq_type>&& v, size_t chunk_size)
                :sv(std::move(v))
                ,chunk_size__(chunk_size)
            {}

            /// move constructor
            SeqVec (SeqVec&& s)
                :sv(std::move(s.sv))
                ,chunk_size__(s.chunk_size__)
            {}

            /// move assignment
            SeqVec& operator=(SeqVec&& s)
            {
                sv = std::move(s.sv);
                chunk_size__ = s.chunk_size__;
                return *this;
            }
            
            /// show how many base inside
            inline size_t size() const
            {
                size_t sum = 0;
                for(size_t i = 0; i != sv.size(); ++i)
                    sum += sv[i].size();
                return sum;
            }

            inline size_t vec_size() const
            {
                return sv.size();
            }

            /** @brief get "chunk"th chunk
             *
             *  it is the original function of [] in vector
             */
            inline seq_type& get_chunk(size_t chunk)
            {
                return sv[chunk];
            }

            /// same as emplace_back of std::vector
            template<class... Args>
            inline void emplace_back(Args&&... args)
            {
                sv.emplace_back(std::forward<Args>(args)...);
            }

            /** @brief get "num"th nucleotide
             *
             * it is the main function of this class 
             * which can find the nucleotide by index. 
             */
            char operator[](size_t num) const
            {
                return sv[num / (chunk_size__)][num & ((chunk_size__) - 1)];
            }

            /** @brief let a string to be seq_type with chunk
             *
             * divide a string into std::vector<seq_type>
             */
            template <class T>
            inline void set_vec(const T& str)
            {
                for(size_t i = 0, cnt = 0; i < str.size(); i += sv[cnt++].size())
                {
                    if(i + chunk_size__ < str.size())
                        sv.emplace_back(str.substr(i, chunk_size__));
                    else
                        sv.emplace_back(str.substr(i, str.size() - i));
                }
            }

            /// complement single base
            char complement(char ch)
            {
                switch(ch)
                {
                    case('A'):
                        return 'T';
                    case('T'):
                        return 'A';
                    case('C'):
                        return 'G';
                    case('G'):
                        return 'C';
                }
            }

            /** @fn reverse complement
             *  reverse and complement of one sequence
             */
            void reverse_complement()
            {

                // reverse
                for(int i = 0, j = sv[sv.size()- 1].size() - 1, cnt_i = 0, cnt_j = sv.size() - 1
						, cnt = 0; cnt != this -> size() / 2; ++i, --j, ++cnt)
                {
                    if(i == sv[cnt_i].size())
                    {
                        i = 0;
                        ++cnt_i;
                    }
                    if(j == -1)
                        j = sv[--cnt_j].size() - 1;

                    char tmp = sv[cnt_i][i];
                    sv[cnt_i][i] = sv[cnt_j][j];
                    sv[cnt_j][j] = tmp;
                }

                // complement
                for(size_t i = 0; i != sv.size(); ++i)
                    if constexpr(std::is_same<std::string, seq_type>::value)
                    {
                        char r;
                        std::replace_if(sv[i].begin(), sv[i].end(), [&](char c)
                                { 
                                    r = complement(c); 
                                    return (c == 'A') || (c == 'C') || (c == 'G') || (c == 'T');
                                }, r);
                    }
                    else
                        sv[i].flip();
            }
        };

        /// A sequence of one chromosome that is divided into many chunk.
        SeqVec seq_;
        /// the name of the chromosome
        std::string name_;

        /** @brief the buffer of istream.
         *
         * If there are still imformation after the end of a chromosome, store
         * it and read it in the next chrmosome
         */
        static std::string seq_buffer_;

        /// default constructor. The default chunk_size_ is 4096    
        Fasta ()
            :name_ ("")
            ,seq_ (SeqVec(chunk_size_))
            ,chunk_size_(4096)
        {}

        /// constructor. Input name, seq and chunk_size
        Fasta (const std::string& n, const std::string& s
                , size_t chunk_size = 2)
            :name_ (n)
            ,seq_ (s, chunk_size_)
            ,chunk_size_(std::pow(2, chunk_size))
        {}
        
        /// constructor. Input name, seq and chunk_size
        Fasta (const char* n, const char* s, size_t chunk_size = 2)
            :name_ (n)
            ,seq_ (s, chunk_size_)
            ,chunk_size_(std::pow(2, chunk_size))
        {}

        ///copy constructor
        Fasta (const Fasta& k)
            :name_(k.name_)
            ,seq_(k.seq_)
            ,chunk_size_(k.get_chunk_size())
        {}

        ///copy assignment
        Fasta& operator= (const Fasta& k)
        {
            name_ = k.name_;
            seq_ = k.seq_;
            chunk_size_ = k.get_chunk_size();
            return *this;
        }

        ///move constructor
        Fasta (Fasta&& k)
            :name_(std::move(k.name_))
            ,seq_(std::move(k.seq_))
            ,chunk_size_(k.get_chunk_size())
        {}
            

        ///move assignment
        Fasta& operator= (Fasta&& k)
        {
            name_ = std::move(k.name);
            seq_ = std::move(k.seq);
            chunk_size_ = k.get_chunk_size();

            return *this;
        }

        /// get the private member-chunk_size_
        inline size_t get_chunk_size() const
        {
            return chunk_size_;
        }

        /** @brief load one chromosome from file
         *
         * input an istream, Fasta (and buffer) and it will 
         * load a chromosome from file into Fasta. If there
         * is nothing to read in fileafter get_obj, it will 
         * return an eof.
         */
        static inline std::istream& get_obj(std::istream& ifile
                , Fasta<seq_type>& fa, std::string& seq_buffer = seq_buffer_)
        {
            size_t buffer_size = std::pow(2, 12);

            bool is_name = true;
            
            // if read the end of one chromosome, read_end = true
            bool read_end = false;
            
            // the count of chunk
            size_t cnt = 0;
            size_t compress_chunk_size = fa.get_chunk_size() << 2;

            fa.seq_.emplace_back();
            while(!read_end)
            {
                std::string buffer;

                // load to buffer
                if(seq_buffer.empty())
                {
                    buffer = std::string(buffer_size, ' ');
                    ifile.read(&buffer[0], buffer_size);
                }

                // if it is not empty in buffer, load it and clear
                else
                {
                    buffer = std::string(seq_buffer);
                    seq_buffer.clear();
                }

                uint32_t beg = 0;

                for(size_t i = 0; i != buffer.size(); ++i)
                {
                    if(is_name)
                    {
                        switch(buffer[i])
                        {
                            case('\n'):
                                is_name = false;
                                beg = i + 1;
                                break;

                            case('\r'):    case('>'):
                                break;
                            
                            default:
                                fa.name_ += buffer[i];
                                break;
                        }
                    }

                    else
                    {
                        // end of buffer or end of one chromose
                        if((i == buffer.size() - 1) || (buffer[i] == '>') 
                                || (buffer[i] == ' '))
                        {
                            std::string tmp;
                            std::copy_if(buffer.begin() + beg, 
                                buffer.begin() + i + 1, std::back_inserter(tmp),
                                [](char x){ return ((x != '\r') && (x != '\n') 
                                && (x != '>') && (x != ' ')); });
                            
                            if(fa.seq_.get_chunk(cnt).size() + tmp.size() >= compress_chunk_size)
                            {
                                size_t tail_size = compress_chunk_size - fa.seq_.get_chunk(cnt).size();
                                fa.seq_.get_chunk(cnt++) += tmp.substr(beg, tail_size);
                                fa.seq_.emplace_back(tmp.substr(tail_size, tmp.size() - tail_size));
                            }

                            else
                                fa.seq_.get_chunk(cnt) += tmp;

                            // end of one chromosome, if it is at the end of buffer, buffer will be read next round
                            if(buffer[i] == '>' || buffer[i] == ' ')
                            {
                                is_name = true;
                                seq_buffer = std::string(buffer, i + 1, buffer.size() - i);
                                read_end = true;
                                break;
                            }
                        }

                        continue;
                    }
                }    
            }

            return ifile;
        }

        /// write the fasta class information back to file
        static inline void store(std::ostream& os, Fasta<seq_type>& fa)
        {
            os << fa.name_ << '\r' << '\n';
            for(size_t i = 0; i != fa.seq_.vec_size(); ++i)
                os << std::string(fa.seq_.get_chunk(i).begin(), fa.seq_.get_chunk(i).end());
        }

        /// write the vector<fasta> information back to file
        static inline void store(std::ostream& os, std::vector<Fasta<seq_type>>& vec_fa)
        {
            for(size_t i = 0; i != vec_fa.size(); ++i)
            {
                os << vec_fa[i].name_ << '\r' << '\n';
                
                for(size_t j = 0; j != vec_fa[i].seq_.vec_size(); ++j)
                    os << std::string(vec_fa[i].seq_.get_chunk(j).begin()
                            , vec_fa[i].seq_.get_chunk(j).end());

                os << '\r' << '\n';
            }
        }

    /// same as get_obj
    friend std::istream& operator>> (std::istream& is, Fasta<seq_type>& fa)
    {
        get_obj(is, fa, seq_buffer_);
    }

    /// same as store
    friend std::ostream& operator<< (std::ostream& os, Fasta<seq_type>& fa)
    {
        store(os, fa);
    }

    /// same as store
    friend std::ostream& operator<< (std::ostream& os
            , std::vector<Fasta<seq_type>>& vec_fa)
    {
        store(os, vec_fa);
    }
    };

    template <class seq_type>
    std::string Fasta<seq_type>::seq_buffer_ = std::string();
}
}

