/**
 *  @file sam.hpp
 *  @brief The program parsing SAM format file
 *  @author JHHlab corp
 */
#pragma once

#include <Biovoltron/format/sam/header.hpp>

namespace biovoltron::format {
    
    class BAM;
    
    /**
     * @brief The main structure of SAM parser.
     * 
     * A SAM class contains a shared header and a 
     * record of data.There are serveral function
     * provided for user to manipulate or to 
     * check the data.
     */
    class SAM
    {
    friend class BAM;
      public:
        using Header = sam::Header;
        /**
         * @brief The normal constructor of SAM. 
         * @param header The reference of the header.
         *
         * Constructing a SAM class via existent header.
         * All the SAM class share the same header.
         */
        SAM(Header& header)
        : header_       (header)
        , data_members_ ( sam::MemberType {
                            std::string()
                          , 0
                          , std::string()
                          , 0
                          , 0
                          , std::string()
                          , std::string()
                          , 0
                          , 0
                          , std::string()
                          , std::string()
                          , std::vector<sam::OptionalFieldType>()
                        })
        {
            std::get<sam::MEMBER_INDEX::RNAME>(data_members_).reserve(23);
            std::get<sam::MEMBER_INDEX::CIGAR>(data_members_).reserve(6);
            std::get<sam::MEMBER_INDEX::SEQ>(data_members_).reserve(35);
            std::get<sam::MEMBER_INDEX::QUAL>(data_members_).reserve(35);
            std::get<sam::MEMBER_INDEX::OPTIONAL_FIELDS>(data_members_)
                .reserve(6);
        };

        /**
         * @brief Convert constructor from BAM to SAM. 
         * @param bam Source BAM object to construct SAM.
         *
         * Constructing a SAM class from BAM object. 
         * Header in this object will reference to bam's header.
         */
        SAM(const BAM& bam);

        /// Function template used to get the specific member of SAM.
        template <std::size_t N>
        const auto& get_member() const
        {
            return std::get<N>(data_members_);
        };

        /// Function template used to set the specific member of SAM.
        template <std::size_t N>
        void set_member
            (const std::tuple_element_t<N, sam::MemberType>& rhs)
        {
            std::get<N>(data_members_) = rhs;
        }; 

        /// Return the header of SAM class.
        const Header& get_header() const
        {
            return header_;
        };

        /// The public interface to user.
        std::string to_string() const
        {
            return to_string<0>();
        }

        /**
         * @brief Get one data from istream.
         * @param is Istream reference containing data.
         * @param obj Target SAM class to be parsed to. 
         * @return Return istream reference for continuous usage.
         *
         * Static member function used to get one data from istream
         * and parse it to SAM class.
         */
        static std::istream& get_obj(std::istream& is, SAM& obj)
        {
            std::string str;
            getline(is, str);
            if (str.size() != 0)
            {
                sam::VecStr vec;
                vec.reserve(26);
                char* c_str = new char[str.size() + 1];
                char* savestr = nullptr;
                std::memcpy(c_str, str.c_str(), str.size() + 1);
                char* split = strtok_r(c_str, "\t ", &savestr);
                while (split != nullptr)
                {
                    vec.emplace_back(split);
                    split = strtok_r(nullptr, "\t ", &savestr);
                }
                parse_one_line_data<0>(vec, obj);
            }
            return is;
        };

        /**
         * @brief Output the whole information of the SAM file.
         * @param os Ostream reference to output data.
         * @param obj Vector of SAM class.
         *
         * Dump the whole data,including one header and several data.
         */
        static void dump(std::ostream& os, std::vector<SAM>& obj)
        {
            os << obj[0].header_;
            for(std::size_t i(0) ; i<obj.size() ; ++i)
            {
                os << obj[i];
                if(i != obj.size()-1)
                    os << '\n';
            }
        };

        /**
         * @brief Overloaded >> operator.
         * @param in Istream reference containing data.
         * @param rhs Target SAM class to be parsed to.
         * @return Istream reference for continuous usage.
         *
         * This overloaded operator is used to get data from istream
         * to the sam class.
         */
        friend std::istream& operator>> (std::istream& in, SAM& rhs)
        {
            get_obj(in, rhs);
            return in;
        };

        /**
         * @brief Overloaded << operator.
         * @param out Ostream reference to output data.
         * @param rhs Target SAM class to be printed.
         * @return Ostream reference for continuous usage.
         *
         * This overloaded operator is used to output data from
         * the SAM class.
         */
        friend std::ostream& operator<< 
            (std::ostream& out, const SAM& rhs)
        {
            out << rhs.to_string<0>();
            return out;
        };

      private:
        Header& header_;                ///< Reference to the header 
                                        ///< which this alignment 
                                        ///< belongs to.
        sam::MemberType data_members_;  ///< A Tuple which 
                                        ///< contains alignment 
                                        ///< information
    
        /// Function template used to recursively 
        /// print the content of SAM.
        template<std::size_t N>
        void print_all()
        {
            if constexpr( N == sam::OPTIONAL_FIELDS )
            {
                std::cout << opt_to_string() << '\n';
                return;
            }
            else
            {
                std::cout << std::get<N>(data_members_) << " ";
                print_all<N+1>();
            }
        }

        /// Function template used to convert optional field 
        /// in data into string.
        std::string opt_to_string() const
        {
            std::string ret;
            for( auto& item : std::get<sam::OPTIONAL_FIELDS>(data_members_) )
            {
                std::array<char,2>tag = std::get<sam::OPTIONAL_FIELD_INDEX::TAG>(item);
                ret.append( 1, tag[0] );
                ret.append( 1, tag[1] );
                ret.append(":");
                ret.append( 1, std::get<sam::OPTIONAL_FIELD_INDEX::VALUE_TYPE>(item));
                ret.append(":");
                ret.append( std::get<sam::OPTIONAL_FIELD_INDEX::VALUE>(item));
                ret.append("\t");
            }
            if( ret.size() != 0 )
                ret.erase(ret.size()-1);
            return ret;
        }

        /**
         * @brief Function template used to parse one line data 
                  to target object.
         * @param vec Vecstor of string contsaining one line data.
         * @param obj Target SAM class to be parsed to.
         * 
         * It's a private function to recursively parse data to a sam 
         * class. It invokes another private function(parse_optional).
         * The reason to declare is as static is that it will be 
         * called by get_obj(...) which is a static member function. 
         */
        template<std::size_t N>
        static void parse_one_line_data(sam::VecStr& vec, SAM& obj)
        {
            if constexpr( N == sam::OPTIONAL_FIELDS )
                parse_optional(vec, obj);
            else
            {
                if constexpr( std::is_integral<
                                std::tuple_element_t<N,
                                    decltype(obj.data_members_)>
                                >::value )
                {
                    std::int32_t res = 0;
                    if(vec[N].size() != 0){
                        auto i = vec[N].begin();
                        boost::spirit::qi::parse(i, vec[N].end(), 
                            boost::spirit::int_, res);
                    }
                    std::get<N>(obj.data_members_) = res;
                }
                else
                    std::get<N>(obj.data_members_) = vec[N];
                parse_one_line_data<N+1>(vec, obj);
            }
            return;
        }

        /// Static function used to parse optional field 
        /// from istream to target sam class.
        static void parse_optional(sam::VecStr& vec, SAM& obj)
        {
            std::get<sam::OPTIONAL_FIELDS>(obj.data_members_).clear();
            for(std::size_t i(11) ; i<vec.size() ; ++i)
            {
                std::string tmp = vec[i];
                std::array<char,2> tag;
                tag[0] = tmp[0];
                tag[1] = tmp[1];
                char vt = tmp[3];
                std::get<sam::OPTIONAL_FIELDS>(obj.data_members_).
                    emplace_back(tag, vt, tmp.substr(5));
            }
        }
        
        /**
         * @brief Function template used to convert sam data into string.
         * 
         * Recursively add data to a string a return it.
         * The returned string represents one data in the sam file.
         */
        template<std::size_t N>
        std::string to_string() const
        {
            std::string ret;
            if constexpr( N == sam::OPTIONAL_FIELDS )
            {
               ret.append(opt_to_string());
            }
            else{  
                if constexpr(
                    std::is_same<
                        std::tuple_element_t< N, sam::MemberType>
                      , std::int32_t
                    >::value )
                    ret.append(
                        std::to_string(std::get<N>(data_members_)));
                else
                    ret.append(std::get<N>(data_members_));
                ret.append("\t");
                ret.append(to_string<N+1>());
            }
            if constexpr ( N == sam::QUAL )
            {
                if (ret.back() == '\t')
                    ret.pop_back();
            }
            return ret;
        };
    };
};
