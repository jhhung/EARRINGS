/**
 *  @file header.hpp
 *  @brief The program parsing SAM's header
 *  @author JHHlab corp
 */
#pragma once

#include <map>
#include <tuple>
#include <vector>
#include <string>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <array>
#include <memory>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <type_traits>
#include <assert.h>
#include <zlib.h>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/include/qi_numeric.hpp>

namespace biovoltron::format
{
    class SAM;
};

namespace biovoltron::format::sam
{
    /// The abbreviation of a vector of string.
    using VecStr = std::vector<std::string>;
    
    /// The abbreviation of a pair of string.
    using PairStr = std::pair<std::string, std::string>;
    
    /// Enumerating different types of sorting.
    enum class SORT_TYPE
    {
        UNKNOWN 
      , UNSORTED
      , QUERYNAME
      , COORDINATE
    };

    /// Mapping a string to a SORT_TYPE.
    inline std::map< std::string, SORT_TYPE > str_to_sort = 
    {
        { "unknown",    SORT_TYPE::UNKNOWN      }
      , { "unsorted",   SORT_TYPE::UNSORTED     }
      , { "queryname",  SORT_TYPE::QUERYNAME    }
      , { "coordinate", SORT_TYPE::COORDINATE   }
    };

    /// Enumerating different types of grouping. 
    enum class GROUP_TYPE 
    {
        NONE
      , QUERY
      , REFERENCE
    };

    /// Mapping a string to a GROUP_TYPE.
    inline std::map<std::string, GROUP_TYPE> str_to_group = 
    {
        { "none",     GROUP_TYPE::NONE      }
      , { "query",    GROUP_TYPE::QUERY     }
      , { "reference",GROUP_TYPE::REFERENCE }
    };
    
    /// Enumerating different types of platforms. 
    enum class PLATFORM
    {
        CAPILLARY
      , LS454
      , ILLUMINA
      , SOLID
      , HELICOS
      , IONTORRENT
      , ONT
      , PABIO
    };
    
    /// Mapping a string to a PLATFORM type.
    inline std::map<std::string, PLATFORM> str_to_platform = 
    {
        { "CAPILLARY",   PLATFORM::CAPILLARY    }
      , { "LS454",       PLATFORM::LS454        }
      , { "ILLUMINA",    PLATFORM::ILLUMINA     }
      , { "SOLID",       PLATFORM::SOLID        }
      , { "HELICOS",     PLATFORM::HELICOS      }
      , { "IONTORRENT",  PLATFORM::IONTORRENT   }
      , { "ONT",         PLATFORM::ONT          }
      , { "PABIO",       PLATFORM::PABIO        }
    };
    
    /**
     *  @brief Indices of shift amount of FLAG field
     *
     *  This is a enumerator to record how many bits should
     *  FLAG field in SAM/BAM right shift to get specific FLAG in this 
     *  alignment. After shifting, LSB is the FLAG you desired.<BR>
     *  ex. (sam.getMember<sam::SAM_MEMBER_INDEX::FLAG>() >> 
     *  sam::FLAG::UNMAPPED) % 2 == 0
     */
    enum FLAG
    {
        MULTI_SEG
      , EXACT_MATCH
      , UNMAPPED
      , NEXT_UNMAPPED
      , REVERSED
      , NEXT_REVERSED
      , FIRST_SEG
      , LAST_SEG
      , SECONDARY
      , QT_FAILED
      , DUPLICATED
      , SUPPLEMENTARY
    };

    /// Indices of SQ field in header.
    enum REFERENCE_INDEX
    {
        REFERENCE_NAME
      , REFERENCE_LENGTH
      , ALTERNATE_LOCUS
      , ALTERNATE_REFERENCE_NAME
      , GENOME_ASSEM_ID
      , SPECIES
    };
    
    /// Indices of RG field in header.
    enum READ_GROUP_INDEX
    {
        READ_GROUP_ID
      , BARCODE
      , READ_GROUP_DESCRIPTION
      , FLOW_ORDER
      , KEY_SEQ
      , LIBRARY
      , READ_GROUP_PLATFORM
    };

    /// Indices of PG field in header.
    enum PROGRAM_INDEX
    {
        PROGRAM_ID
      , PROGRAM_NAME
      , COMMAND
      , PROGRAM_DESCRIPTION
      , PROGRAM_VERSION
    }; 

    /// Indices of the whole header.
    enum HEADER_INDEX
    {
        VERSION
      , ALIGNMENT_SORT_ORDER
      , ALIGNMENT_GROUPING
      , REFERENCE
      , READ_GROUP
      , PROGRAM
      , COMMENT
      , PLAIN_TEXT
    }; 

    /// Indices of the optional field in SAM class.
    enum OPTIONAL_FIELD_INDEX
    {
        TAG
      , VALUE_TYPE
      , VALUE
    };

    /// Indices of the member in SAM class.
    enum MEMBER_INDEX
    {
        QNAME
      , FLAG
      , RNAME
      , POS
      , MAPQ
      , CIGAR
      , RNEXT
      , PNEXT
      , TLEN
      , SEQ
      , QUAL
      , OPTIONAL_FIELDS
    };

    /// Data structure of SQ field in header class.
    using ReferenceType = std::tuple<
        std::string
      , std::int32_t
      , std::string
      , std::string
      , std::string
      , std::string
    >;

    /// Data structure of RG field in header class.
    using ReadGroupType = std::tuple<
        std::string
      , std::string
      , std::string
      , std::string
      , std::string
      , std::string
      , PLATFORM
    >;
    
    /// Data structure of PG field in header class.
    using ProgramType = std::tuple<
        std::string
      , std::string
      , std::string
      , std::string
      , std::string
    >;
    
    /// Data structure of the whole header class.
    using HeaderType = std::tuple<
        std::string
      , SORT_TYPE
      , GROUP_TYPE
      , std::vector<ReferenceType>
      , std::vector<ReadGroupType>
      , std::vector<ProgramType>
      , std::vector<std::string>
      , std::string
    >; 

    /// Data structure of the optional field in SAM class.
    using OptionalFieldType = std::tuple<
        std::array<char,2>
      , char
      , std::string
    >;

    /// Data structure of the member data in SAM class.
    using MemberType = std::tuple<
        std::string                     ///< QRNAME
      , std::int32_t                    ///< FLAG   
      , std::string                     ///< RNAME
      , std::int32_t                    ///< POS
      , std::int32_t                    ///< MAPQ
      , std::string                     ///< CIGAR
      , std::string                     ///< RNEXT
      , std::int32_t                    ///< PNEXT
      , std::int32_t                    ///< TLEN
      , std::string                     ///< SEQ
      , std::string                     ///< QUAL
      , std::vector<OptionalFieldType>  ///< OPTIONAL FIELD
      >;
               
    /**
     *  @brief The class of SAM header.
     *  
     *  This class is essential to sam class.
     *  Which means that every sam class must have
     *  a header class.And each sam class will share
     *  the same header.
     */
    class Header
    {
    friend class ::biovoltron::format::SAM;
      public:
        /// Default constructor of header class.
        Header()
        : header_     ( HeaderType {
                        std::string()  
                      , SORT_TYPE::UNKNOWN
                      , GROUP_TYPE::NONE
                      , std::vector<ReferenceType>()
                      , std::vector<ReadGroupType>()
                      , std::vector<ProgramType>()
                      , std::vector<std::string>()
                      , std::string()
                      })
        {
            std::get<HEADER_INDEX::REFERENCE>(header_).reserve(3);
            std::get<HEADER_INDEX::COMMENT>(header_).reserve(1);
            std::get<HEADER_INDEX::PLAIN_TEXT>(header_).reserve(50);
        };

        /** 
         * @brief Istream constructor of header class.
         * @param in Istream containing header information.  
         *
         * Constructing Header class through istream.
         * This constuctor will parse the whole header
         * in a sam file.
         */
        Header(std::istream& in) 
        {
            preparse(in);
        };

        /// Function template used to get a specific field in header.
        template <std::size_t n>
        const auto& get_member() const
        {
            return std::get<n>(header_);
        };

        /** 
         * @brief Overloaded >> operator. 
         * @param in Istream reference constaining header information.
         * @param rhs Target header instance.
         * @return Reference of istream.
         *
         * This overloaded operator is used to parse the
         * information in istream to a header. 
         */
        friend std::istream& operator>> 
            (std::istream& in, Header& rhs)
        {
            rhs.preparse(in);
            return in;
        }

        /**
         * @brief Overloaded << operator.
         * @param out Ostream reference outputing header content. 
         * @param rhs Target header instance.
         * @return Reference of ostream.
         *
         * This overloaded operator is used to output the
         * information of a header.
         */
        friend std::ostream& operator<< 
            (std::ostream& out, const Header& rhs)
        {
            out << rhs.to_string();
            return out;
        }
        
        /// Out put the whole header information.
        std::string to_string() const
        {
            return std::get<PLAIN_TEXT>(header_);
        }

        /**
         * @brief Public interface to pre-parse a header. 
         * @param in Istream reference containing header information.
         *
         * This function is to recognize headers,classify them
         * and invoke corresponding handling function. 
         */
        void preparse(std::istream& in)
        {
            VecStr lines;
            std::string line;
            int c = in.peek();
            while( c == '@')
            {
                getline(in, line);
                lines.emplace_back(std::move(line));
                c = in.peek();
            }
            preparse_impl(lines);
        };

      protected:
        HeaderType header_; ///< This is a protected data member 
                            ///< because bam::Header has exactly 
                            ///< the same header contents

        /**
         * @brief This is a private function to parse header one by one.
         * @param lines Lines which splited from header
         *
         * This function is implementation details for parsing headers. 
         * Set to protected member function is because 
         * bam::Header also need it.
         */
        void preparse_impl(VecStr& lines)
        {
            VecStr vec;
            std::string& text = std::get<PLAIN_TEXT>(header_);
            for(std::size_t i = 0;i < lines.size();++i)
            {
                boost::split(vec,lines[i], boost::is_any_of(" \t"), 
                                boost::token_compress_on);
                if( vec[0] == "@HD" )
                    parse_HD(vec);
                else if( vec[0] == "@SQ" )
                    parse_SQ(vec);
                else if( vec[0] == "@RG" )
                    parse_RG(vec);
                else if( vec[0] == "@PG" )
                    parse_PG(vec);
                else if( vec[0] == "@CO" )
                    parse_CO(vec);
                else
                    ;
                text.append(lines[i]);
                text.append("\n");
            }
        }

      private:
        /// Parsing HD field to header.
        void parse_HD(VecStr& vec)
        {
            for( std::size_t i(1) ; i<vec.size() ; ++i )
            {
                PairStr tmp = split_colon(vec[i]);
                if( tmp.first == "VN" )
                    std::get<VERSION>(header_) = tmp.second;
                
                else if( tmp.first == "SO" )
                    std::get<ALIGNMENT_SORT_ORDER>(header_) = 
                        str_to_sort[tmp.second];
                
                else if( tmp.first == "GO" )
                    std::get<ALIGNMENT_GROUPING>(header_) = 
                        str_to_group[tmp.second];
                
                else
                    ;
            }
        }

        /// Parsing SQ field to header.
        void parse_SQ(VecStr& vec)
        {
            ReferenceType ret;
            for( std::size_t i(1) ; i<vec.size() ; ++i )
            {
                PairStr tmp = split_colon(vec[i]);
                if( tmp.first == "SN" )
                {
                    std::get<REFERENCE_NAME>(ret) = tmp.second;
                }
                else if( tmp.first == "LN" )
                {
                    std::get<REFERENCE_LENGTH>(ret) = 
                        std::stoi(tmp.second);
                }
                else if( tmp.first == "AH" )
                {
                    std::get<ALTERNATE_LOCUS>(ret) = tmp.second;
                }
                else if( tmp.first == "AN" )
                {
                    std::get<ALTERNATE_REFERENCE_NAME>(ret) = tmp.second;
                }
                else if( tmp.first == "AS" )
                {
                    std::get<GENOME_ASSEM_ID>(ret) = tmp.second;
                }
                else if( tmp.first == "SP" )
                {
                    std::get<SPECIES>(ret) = tmp.second;
                }
                else
                    ;
            }
            std::get<REFERENCE>(header_).emplace_back(ret);
        }

        /// Parsing RG field to header.
        void parse_RG(VecStr& vec)
        {
            ReadGroupType ret;
            for( std::size_t i(1) ; i<vec.size() ; ++i )
            {
                PairStr tmp = split_colon(vec[i]);
                if( tmp.first == "ID" )
                {
                    std::get<READ_GROUP_ID>(ret) = tmp.second;
                }
                else if( tmp.first == "BC" )
                {
                    std::get<BARCODE>(ret) = tmp.second;
                }
                else if( tmp.first == "DS" )
                {
                    std::get<READ_GROUP_DESCRIPTION>(ret) = tmp.second;
                }
                else if( tmp.first == "FO" )
                {
                    std::get<FLOW_ORDER>(ret) = tmp.second;
                }
                else if( tmp.first == "KS" )
                {
                    std::get<KEY_SEQ>(ret) = tmp.second;
                }
                else if( tmp.first == "LB" )
                {
                    std::get<LIBRARY>(ret) = tmp.second;
                }
                else if( tmp.first == "PL" )
                {
                    std::get<READ_GROUP_PLATFORM>(ret) = 
                        str_to_platform[tmp.second];
                }
                else
                    ;
            }
            std::get<READ_GROUP>(header_).emplace_back(ret);
        }

        /// Parsing PG field to header.
        void parse_PG(VecStr& vec)
        {
            ProgramType ret;
            for( std::size_t i(1) ; i<vec.size() ; ++i )
            {
                PairStr tmp = split_colon(vec[i]);
                if( tmp.first == "ID" )
                {
                    std::get<PROGRAM_ID>(ret) = tmp.second;
                }
                else if( tmp.first == "PN" )
                {
                    std::get<PROGRAM_NAME>(ret) = tmp.second; 
                }
                else if( tmp.first == "CL" )
                {
                    std::get<COMMAND>(ret) = tmp.second;
                }
                else if( tmp.first == "DS" )
                {
                    std::get<PROGRAM_DESCRIPTION>(ret) = tmp.second;
                }
                else if( tmp.first == "VN" )
                {
                    std::get<PROGRAM_VERSION>(ret) = tmp.second;
                }
                else
                    ;
            }
            std::get<PROGRAM>(header_).emplace_back(ret);
        }

        /// Parsing CO field to header.
        void parse_CO(VecStr& vec)
        {
            std::string ret;
            for( std::size_t i(1) ; i<vec.size() ; ++i )
            {
                ret.append(vec[i]);
                if( i != vec.size()-1 )
                    ret.append("\t");
            }
            std::get<COMMENT>(header_).emplace_back(ret);
        }
        
        /**
         * @brief Function to split string to tokens by colon 
         * @param A reference to string 
         */
        PairStr split_colon(std::string& str)
        {
            PairStr ret;
            std::size_t pos = str.find(":");
            ret.first = str.substr(0, pos);
            ret.second = str.substr(pos + 1);
            return ret;
        }
    };
};
