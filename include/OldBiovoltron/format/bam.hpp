/**
 *  @file bam.hpp
 *  @brief The program parsing BAM format file
 *  @author JHHlab corp
 */
#pragma once

#include <OldBiovoltron/format/sam.hpp>
#include <OldBiovoltron/format/bam/header.hpp>

namespace biovoltron::format {

    class BAM;
    using namespace bam;
    
    /**
     * @brief Class to store BAI file index information
     * @warning Only BAI file read successfully then could active 
     * random access or an error might occur.
     *
     * A class to store indexing information in BAI file to 
     * accomplish fast random access in BAM file.<BR>
     * Class also provides functions to set user's interested
     * reference region to help locate right file offset. 
     */
    class BAI
    {
    friend class BAM;
      public:
        /**
         *  @brief Filename constructor to read BAI file.
         *  @param filename BAI file's filename
         *  @see load()
         * 
         *  Construct object and load index information simultaneously. 
         *  Equlas to call load() to try to load BAI file into object.<BR>
         *  This is the only way to construct BAI file. An empty 
         *  BAI object is meanless.
         */
        BAI(std::string filename)
        : region_           ( RegionType{-1, -1, -1, -1} )
        , current_index_    (             -1             )
        , n_no_coor_        (              0             )
        {
            load(filename);   
        };

        BAI(const BAI& rhs) = default;
        BAI(BAI&& rhs) = default;
        BAI& operator=(const BAI& rhs) = default;
        BAI& operator=(BAI&& rhs) = default;
        ~BAI() = default;

        /// Clear all data in this object
        void reset()
        {
            index_info_.clear();
            overlap_chunks_.clear();
            region_ = RegionType{-1, -1, -1, -1};
            current_index_ = -1;
            n_no_coor_ = 0;
        }

        /**
         *  @brief Load index information from BAI file.
         *  @param filename BAI file's filename
         *  @return Whether information load sucessfully
         */
        bool load(std::string filename)
        {
            reset();
            std::ifstream infile(filename, std::ios::binary);
            if (!infile)
            {
                std::cerr << "ERROR: File not open" << std::endl;
                return false;
            }
            // magic
            char magic[4];
            infile.read(magic, 4);
            if (std::memcmp(magic, BAI_MAGIC, 4) != 0)
            {
                std::cerr << "ERROR: BAI magic string not match";
                std::cerr << std::endl;
                return false;
            }
            // n_ref
            std::int32_t n_ref, n_bin, n_chunk, n_intv;
            std::uint32_t bin_num;
            std::uint64_t chunk_beg, chunk_end, ioffset;
            BinType bin_buf;
            LinearIndexType index_buf;
            std::vector<ChunkType> chunks_buf;
            infile.read(reinterpret_cast<char*>(&n_ref)
                        , sizeof(std::int32_t));
            for (std::size_t i = 0;i < n_ref;++i)
            {
                // n_bin
                infile.read(reinterpret_cast<char*>(&n_bin)
                            , sizeof(std::int32_t));
                for (std::size_t j = 0;j < n_bin;++j)
                {
                    // bin
                    infile.read(reinterpret_cast<char*>(&bin_num)
                                , sizeof(std::uint32_t));
                    // n_chunk
                    infile.read(reinterpret_cast<char*>(&n_chunk)
                                , sizeof(std::int32_t));
                    for (std::size_t k = 0;k < n_chunk;++k)
                    {
                        infile.read(reinterpret_cast<char*>(&chunk_beg)
                                    , sizeof(std::uint64_t));
                        infile.read(reinterpret_cast<char*>(&chunk_end)
                                    , sizeof(std::uint64_t));
                        chunks_buf.emplace_back(chunk_beg, chunk_end);
                    }
                    bin_buf.emplace(bin_num, chunks_buf);
                    chunks_buf.clear();
                }
                // n_intv
                infile.read(reinterpret_cast<char*>(&n_intv)
                            , sizeof(std::int32_t));
                for (std::size_t j = 0;j < n_intv;++j)
                {
                    // ioffset
                    infile.read(reinterpret_cast<char*>(&ioffset)
                                , sizeof(std::uint64_t));
                    index_buf.emplace_back(ioffset);
                }
                index_info_.emplace_back(bin_buf, index_buf);
                index_buf.clear();
                bin_buf.clear();
            }
            infile.read(reinterpret_cast<char*>(&n_no_coor_)
                        , sizeof(std::uint64_t));
            if (infile.peek() != EOF)
            {
                std::cerr << "ERROR: BAI file doesn't read properly";
                std::cerr << std::endl;
                return false;
            }
            return true;
        };

        /**
         *  @brief Test whether index information load successfully.
         *  @return Whether information load successfully
         */
        bool is_useable() const
        {
            return !(index_info_.empty());
        };

        /**
         *  @brief Set interested reference region.
         *  @param r The region is interested
         * 
         *  Set interested reference region 
         *  and seek to nearest block address.<BR>
         *  Will check whether region is valid then use index information 
         *  to seek to a alignment might overlapped with region.
         */
        void set_region(const bam::RegionType& r)
        {
            if (!is_useable())
            {
                std::cerr << "ERROR: Index information is empty.";
                std::cerr << std::endl;
                return;
            }
            std::int64_t left_ref = 
                (std::get<REGION_INDEX::LEFT_REF_ID>(r) <= -1) ? 
                    0 : std::get<REGION_INDEX::LEFT_REF_ID>(r);
            std::int64_t left_pos = 
                (std::get<REGION_INDEX::LEFT_POS>(r) <= -1) ? 
                    0 : std::get<REGION_INDEX::LEFT_POS>(r);
            std::int64_t right_ref = 
                (std::get<REGION_INDEX::RIGHT_REF_ID>(r) <= -1) ? 
                    index_info_.size() - 1 
                    : std::get<REGION_INDEX::RIGHT_REF_ID>(r);
            std::int64_t right_pos = 
                (std::get<REGION_INDEX::RIGHT_POS>(r) <= -1) ? 
                    LONGEST_REFERENCE_LENGTH 
                    : std::get<REGION_INDEX::RIGHT_POS>(r);
            if (left_ref >= index_info_.size())
            {
                std::cerr << "ERROR: Left reference out of range";
                std::cerr << std::endl;
                return;
            }
            else if (right_ref >= index_info_.size())
            {
                std::cerr << "ERROR: Right reference out of range";
                std::cerr << std::endl;
                return;
            }
            else if (left_ref > right_ref || 
                        (left_ref == right_ref && left_pos > right_pos))
            {
                std::cerr << "ERROR: Left position ";
                std::cerr << "is greater than right position";
                std::cerr << std::endl;
                return;
            }

            region_ = r;
            current_index_ = -1;
            std::size_t min_index;  // min_offset's index in lenear_offset
            std::uint64_t min_offset;  // linear offset that include beg
            // buffer region in the middle of processing
            std::int64_t current_left_pos, current_right_pos;
            // overall overlap bins in region
            std::vector<std::uint32_t> overlap_bins;
            BAI::BinType::iterator bins_iter;
            overlap_chunks_.clear();
            for (std::int64_t current_ref = left_ref;
                    current_ref <= right_ref; ++current_ref)
            {
                BAI::BinType& bins = 
                    std::get<BAI::REFERENCE_BIN_INFO_INDEX::BINS>(
                        index_info_[current_ref]);
                std::vector<std::uint64_t>& linear_offset = 
                    std::get<BAI::REFERENCE_BIN_INFO_INDEX::LINEAR_INDEX>(
                        index_info_[current_ref]);
                current_left_pos = 
                    (current_ref == left_ref) ? left_pos : 0;
                current_right_pos = 
                    (current_ref == right_ref) ? 
                        right_pos : LONGEST_REFERENCE_LENGTH;
                min_index = 
                    (current_left_pos >> BAM_LINEAR_SHIFT) > 
                        linear_offset.size() - 1 ?
                            linear_offset.size() - 1 : 
                            current_left_pos >> BAM_LINEAR_SHIFT;
                min_offset = linear_offset[min_index];
                overlap_bins = 
                    reg2bins(current_left_pos, current_right_pos);
                for (std::size_t i = 0;i < overlap_bins.size();++i)
                {
                    bins_iter = bins.find(overlap_bins[i]);
                    if (bins_iter == bins.end()) continue;
                    for (std::size_t j = 0;j < (bins_iter->second).size();++j)
                    {
                        if (std::get<BAI::CHUNK_INDEX::END>(
                                (bins_iter->second)[j]) > min_offset)
                        {
                            overlap_chunks_
                                .emplace_back((bins_iter->second)[j]);
                        }
                    }
                }
            }
            auto chunk_cmp = [](BAI::ChunkType a, BAI::ChunkType b)
            {
                return (std::get<BAI::CHUNK_INDEX::BEG>(a) < 
                        std::get<BAI::CHUNK_INDEX::BEG>(b));
            };
            std::sort(overlap_chunks_.begin()
                    , overlap_chunks_.end()
                    , chunk_cmp);
        };

        /**
         *  @brief Get current interested reference region.
         *  @return Interested reference region
         */
        bam::RegionType get_region()
        {
            return region_;
        };

        /**
         *  @brief Get number of unplaced reads which record in BAI file.
         *  @return The number of unplaced reads
         */
        std::uint64_t get_unplaced_read_num()
        {
            return n_no_coor_;
        }

        /**
         *  @brief Jump to specific reference position.
         *  @param refID Reference ID which position located
         *  @param pos Offset from the beginning of reference
         *  @see set_region()
         * 
         *  Equals to call set_region() which left bound is parameters 
         *  and with opened right bound, 
         *  acting like user jump to a position.
         */
        void jump(std::int64_t refID, std::int64_t pos)
        {
            set_region(RegionType {refID, pos, -1, -1});
        };

        /**
         * @brief Test whether current region has start position.
         * @return Whether current region has start position
         */
        bool is_start_region_specified()
        {
            return !(std::get<REGION_INDEX::LEFT_REF_ID>(region_) == -1 &&
                        std::get<REGION_INDEX::LEFT_POS>(region_) == -1);
        }

        /**
         * @brief Test whether current region has end position.
         * @return Whether current region has end position
         */
        bool is_end_region_specified()
        {
            return !(std::get<REGION_INDEX::RIGHT_REF_ID>(region_) == -1 &&
                        std::get<REGION_INDEX::RIGHT_POS>(region_) == -1);
        }

      private:
        /**
         * @brief Find overlapping bins with specified range.
         * @param beg Starting position
         * @param end Ending position
         * @return std::vector stores overlapping bin numbers
         *
         * Calculate the list of bins that may overlap with region 
         * [beg,end) (zero-based).<BR>
         * Function declared as static is because there has no data 
         * member being used. Everyone can access it without an object.
         */
        static std::vector<std::uint32_t> reg2bins(std::int64_t beg, std::int64_t end)
        {
            int i = 0, k;
            --end;
            std::vector<std::uint32_t> list;
            list.emplace_back(0);
            for (k =    1 + (beg >> 26); k <=    1 + (end >> 26); ++k)
                list.emplace_back(k);
            for (k =    9 + (beg >> 23); k <=    9 + (end >> 23); ++k)
                list.emplace_back(k);
            for (k =   73 + (beg >> 20); k <=   73 + (end >> 20); ++k)
                list.emplace_back(k);
            for (k =  585 + (beg >> 17); k <=  585 + (end >> 17); ++k)
                list.emplace_back(k);
            for (k = 4681 + (beg >> 14); k <= 4681 + (end >> 14); ++k)
                list.emplace_back(k);
            return list;
        }

        /// Indices of ChunkType to access chunk's begin and end
        enum CHUNK_INDEX {BEG, END};
        /// Indices of ReferenceBinInfoType to access 
        /// bins information ans linear index
        enum REFERENCE_BIN_INFO_INDEX {BINS, LINEAR_INDEX};
        /// Data structure to store chunk region
        using ChunkType = std::tuple<std::uint64_t, std::uint64_t>;
        /// Data structure to store all bin's chunk list in a reference
        using BinType = std::map<std::uint32_t, std::vector<ChunkType>>;
        /// Data structure to store linear index in a reference
        using LinearIndexType = std::vector<std::uint64_t>;
        /// Data structure to store summary in a reference 
        /// (include bins and linear index)
        using ReferenceBinInfoType = std::tuple<BinType
                                              , LinearIndexType>;

        /// Bits need to shift to get linear index's index
        const static std::uint8_t BAM_LINEAR_SHIFT = 14;
        /// Record longest reference length that BAM could accept
        const static std::int64_t LONGEST_REFERENCE_LENGTH = 0x1fffffff;
        /// Record BAI file's magic string
        static constexpr char BAI_MAGIC[4] = {'B', 'A', 'I', '\1'};

        /// Records the number of unmapped reads
        std::uint64_t n_no_coor_;
        /// Records all reference binning and linear index information
        std::vector<ReferenceBinInfoType> index_info_;
        /// Records current interested reference range
        RegionType region_;
        /// Records index of overlap_chunks_ now reached
        std::int64_t current_index_;
        /// Stores chunks that might overlapped with region
        std::vector<ChunkType> overlap_chunks_;
    };

    /**
     * @brief Class stores A BAM alignment.
     * @warning A BAM object has a reference to bam::Header.
     * the header MUST NOT be distruct before this object destruct.
     *
     * The main class to maintain a BAM alignment.<BR>
     * It provides functions to get/set alignment information 
     * and get alignment from BAM file format.
     */
    class BAM
    {
    friend class SAM;
    friend class bam::Header;
      public:
        using Header = bam::Header;

        /**
         * @brief Header constructor to reference to right header.
         * @param header Corresponding header with same BAM file
         * @see get_obj()
         *
         * Construct a BAM object which its header is referenced to 
         * the header which is loaded from desired BAM file.<BR>
         * This time BAM object still has no information about alignments
         * until called get_obj(), but header must be intialized first.
         */
        BAM(Header& header)
        : header_       ( header )
        , has_data_     ( false  )
        , data_members_ ( MemberType {
                            std::string()
                          , 0
                          , -1
                          , -1
                          , 255
                          , std::vector<CigarType>()
                          , -1
                          , -1
                          , 0
                          , std::string()
                          , std::string()
                          , std::vector<OptionalFieldType>()
                          , -1
                        })
        {
            std::get<MEMBER_INDEX::QNAME>(data_members_).reserve(23);
            std::get<MEMBER_INDEX::CIGAR>(data_members_).reserve(6);
            std::get<MEMBER_INDEX::SEQ>(data_members_).reserve(35);
            std::get<MEMBER_INDEX::QUAL>(data_members_).reserve(35);
            std::get<MEMBER_INDEX::OPTIONAL_FIELDS>(data_members_)
                .reserve(15);
        };

        /**
         * @brief Convert constructor from SAM object.
         * @param header BAM header constructed by rhs' header
         * @param rhs Source SAM object which want to convert.
         * @warning Parameter header should be constructed by rhs' header
         *          or might mismatch between alignment and header.
         *
         * A convert constructor from SAM to BAM.<BR>
         * Because SAM is a subset of BAM, BAM needs additional
         * information to construct it. So needs a bam::Header 
         * to construct it. Of course, header should available
         * until this object destruct.
         */
        BAM(Header& header, const SAM& rhs)
        : header_   ( header )
        , has_data_ (  true  )
        {
            std::get<MEMBER_INDEX::QNAME>(data_members_) = 
                std::get<sam::MEMBER_INDEX::QNAME>(rhs.data_members_);
            std::get<MEMBER_INDEX::FLAG>(data_members_) = 
                std::get<sam::MEMBER_INDEX::FLAG>(rhs.data_members_);
            // RNAME
            const std::string& ref_name = 
                std::get<sam::MEMBER_INDEX::RNAME>(rhs.data_members_);
            std::vector<ReferenceType>& refs = 
                std::get<sam::HEADER_INDEX::REFERENCE>(header_.header_);
            std::vector<ReferenceType>::iterator refs_it = refs.begin();
            std::int32_t& ref_id = 
                std::get<MEMBER_INDEX::REF_ID>(data_members_);
            if (ref_name == "*")
                ref_id = -1;
            else
            {
                while (refs_it != refs.end())
                {
                    if (std::get<REFERENCE_INDEX::REFERENCE_NAME>(*refs_it)
                         == ref_name)
                        break;
                    std::advance(refs_it, 1);
                }
                ref_id = (refs_it == refs.end()) ? 
                            -1 : std::distance(refs.begin(), refs_it);
            }
            std::int32_t& pos = 
                std::get<MEMBER_INDEX::POS>(data_members_);
            pos = std::get<sam::MEMBER_INDEX::POS>(rhs.data_members_) - 1;
            std::get<MEMBER_INDEX::MAPQ>(data_members_) = 
                std::get<sam::MEMBER_INDEX::MAPQ>(rhs.data_members_);
            // CIGAR
            std::vector<CigarType>& cigar_vec = 
                std::get<MEMBER_INDEX::CIGAR>(data_members_);
            std::string cigar_str = 
                std::get<sam::MEMBER_INDEX::CIGAR>(rhs.data_members_);
            std::vector<OptionalFieldType>& bam_of = 
                std::get<MEMBER_INDEX::OPTIONAL_FIELDS>(data_members_);
            std::uint32_t consume_ref_len = 0;
            cigar_vec = str2cigar_vec(cigar_str, consume_ref_len);
            if (cigar_vec.size() > MAX_CIGAR_OP_NUM)
            {
                LightString value;
                value.append('I');
                std::size_t cigar_vec_size = cigar_vec.size();
                Header::parse_char<std::int32_t>(&value.c_str[1], cigar_vec_size);
                value.size += sizeof(std::int32_t);
                for (std::size_t i = 0;i < cigar_vec_size;++i)
                {
                    if (value.size + sizeof(std::uint32_t) > value.capacity)
                        value.extend();
                    Header::parse_char<std::uint32_t>(&value.c_str[value.size]
                        , std::get<bam::CIGAR_INDEX::LENGTH>(cigar_vec[i]) << 4
                            | std::get<bam::CIGAR_INDEX::OP>(cigar_vec[i]));
                    value.size += sizeof(std::uint32_t);
                }
                bam_of.emplace_back(std::array<char, 2>{'C', 'G'}, 'B'
                                  , std::string(value.c_str, value.size));
                cigar_vec.clear();
                cigar_vec.emplace_back('S'
                                     , std::get<sam::MEMBER_INDEX::SEQ>(
                                        rhs.data_members_).size());
                cigar_vec.emplace_back('N', (ref_id == -1) ? 
                    0 : 
                    std::get<REFERENCE_INDEX::REFERENCE_LENGTH>(refs[ref_id]));
            }
            // RNEXT
            const std::string& nref_name = 
                std::get<sam::MEMBER_INDEX::RNEXT>(rhs.data_members_);
            if (nref_name == "*")
                std::get<MEMBER_INDEX::RNEXT_ID>(data_members_) = -1;
            else if (nref_name == "=")
                std::get<MEMBER_INDEX::RNEXT_ID>(data_members_) = ref_id;
            else
            {
                refs_it = refs.begin();
                while (refs_it != refs.end())
                {
                    if (std::get<REFERENCE_INDEX::REFERENCE_NAME>(*refs_it) 
                            == nref_name)
                        break;
                    std::advance(refs_it, 1);
                }
                std::get<MEMBER_INDEX::RNEXT_ID>(data_members_) = 
                    (refs_it == refs.end()) ? 
                        -1 : std::distance(refs.begin(), refs_it);
            }
            std::get<MEMBER_INDEX::PNEXT>(data_members_) = 
                std::get<sam::MEMBER_INDEX::PNEXT>(rhs.data_members_) - 1;
            std::get<MEMBER_INDEX::TLEN>(data_members_) = 
                std::get<sam::MEMBER_INDEX::TLEN>(rhs.data_members_);
            std::string& seq = 
                std::get<MEMBER_INDEX::SEQ>(data_members_);
            seq = std::get<sam::MEMBER_INDEX::SEQ>(rhs.data_members_);
            // QUAL
            std::string& qual = 
                std::get<MEMBER_INDEX::QUAL>(data_members_);
            const std::string& qual_33 = 
                std::get<sam::MEMBER_INDEX::QUAL>(rhs.data_members_);
            if (qual_33 == "*")
                qual = std::string(seq.size(), '\xFF');
            else
            {
                qual.reserve(qual_33.size());
                for (std::size_t i = 0;i < qual_33.size();++i)
                    qual.push_back(qual_33[i] - 33);
            }
            // Optional fields
            const std::vector<OptionalFieldType>& sam_of = 
                std::get<sam::MEMBER_INDEX::OPTIONAL_FIELDS>(rhs.data_members_);
            char value_type;
            std::int64_t number;
            for (std::size_t i = 0;i < sam_of.size();++i)
            {
                value_type = 
                    std::get<OPTIONAL_FIELD_INDEX::VALUE_TYPE>(sam_of[i]);
                const std::string& value = 
                    std::get<OPTIONAL_FIELD_INDEX::VALUE>(sam_of[i]);
                if (value_type == 'i')
                {
                    boost::spirit::qi::parse(value.begin()
                                           , value.end()
                                           , boost::spirit::long_
                                           , number);
                    if (number <= UCHAR_MAX && number > -1)
                        value_type = 'C';
                    else if(number <= SCHAR_MAX && number >= SCHAR_MIN)
                        value_type = 'c';
                    else if (number <= USHRT_MAX && number > -1)
                        value_type = 'S';
                    else if (number <= SHRT_MAX && number >= SHRT_MIN)
                        value_type = 's';
                    else if (number <= UINT_MAX && number > -1)
                        value_type = 'I';
                    else if (number <= INT_MAX && number >= INT_MIN)
                        value_type = 'i';
                    else
                        value_type = 'Z';
                }
                bam_of.emplace_back(
                    std::get<OPTIONAL_FIELD_INDEX::TAG>(sam_of[i])
                  , value_type
                  , pack_optional_field(value_type, value, number));
            }
            // bin
            std::get<MEMBER_INDEX::BIN>(data_members_) = 
                reg2bin(pos, pos + consume_ref_len);
        }

        BAM(const BAM& rhs) = default;
        BAM(BAM&& rhs) = default;
        ~BAM() = default;

        /**
         * @brief Get a specific alignment field.
         * @tparam n The field index from BAM_MEMBER_INDEX
         * @return Specific alignment field
         * @see BAM_MEMBER_INDEX
         */
        template <std::size_t n>
        const auto& get_member() const
        {
            return std::get<n>(data_members_);
        };

        /**
         * @brief Set a specific alignment field.
         * @tparam n The field index from BAM_MEMBER_INDEX
         * @see BAM_MEMBER_INDEX
         */
        template <std::size_t n>
        void set_member
            (const std::tuple_element_t<n, MemberType>& rhs)
        {
            std::get<n>(data_members_) = rhs;
        };

        /**
         * @brief Get the header which this object reference to.
         * @return Header object which reference to
         */
        const Header& get_header() const
        {
            return header_;
        };

        /**
         * @brief Test this object contains valid alignment information.
         * @return Whether this object contains valid alignment information
         * @warning If this object contains invalid alignment information,
         *          then no assumptions can be made about information that 
         *          get from get_member().
         * @see get_member()
         */
        bool is_valid() const
        {
            return has_data_;
        }

        /**
         * @brief Set this alignment information to valid.
         * @param valid The value user want to set to
         * @warning This function should use carefully or
         * might get invalid data leads to undefined behavior.
         */
        void set_valid(bool valid)
        {
            has_data_ = valid;
        }

        /**
         * @brief Convert BAM object to SAM string format.
         * @return String that has been converted
         * @warning If object is invalid will return empty string.
         * @see is_valid()
         */
        std::string to_string()
        {
            if (!has_data_)
                return std::string();
            LightString tmp;
            to_string_impl<0>(tmp);
            return std::string(tmp.c_str, tmp.size);
        }

        /// Clear all data in this object.
        void reset()
        {
            has_data_ = false;
            std::get<MEMBER_INDEX::FLAG>(data_members_) = 0;
            std::get<MEMBER_INDEX::REF_ID>(data_members_) = -1;
            std::get<MEMBER_INDEX::POS>(data_members_) = -1;
            std::get<MEMBER_INDEX::MAPQ>(data_members_) = 255;
            std::get<MEMBER_INDEX::RNEXT_ID>(data_members_) = -1;
            std::get<MEMBER_INDEX::PNEXT>(data_members_) = -1;
            std::get<MEMBER_INDEX::TLEN>(data_members_) = -1;
            std::get<MEMBER_INDEX::BIN>(data_members_) = 0;
            std::get<MEMBER_INDEX::QNAME>(data_members_).clear();
            std::get<MEMBER_INDEX::CIGAR>(data_members_).clear();
            std::get<MEMBER_INDEX::SEQ>(data_members_).clear();
            std::get<MEMBER_INDEX::QUAL>(data_members_).clear();
            std::get<MEMBER_INDEX::OPTIONAL_FIELDS>(data_members_)
                .clear();
        }

        /**
         * @brief Read alignment information from BAM file.
         * @param in istream which contains the BAM file
         * @param obj The object which want data tpo put into
         * @return Parameter in for continuous data flow
         * @warning Parameter in should be the istream object 
         * which has been used to parse into the header object, 
         * which is referenced by parameter obj.
         *
         * Read an alignment information in BAM file.<BR>
         * After all field are loaded sucessfully, this object 
         * will become valid state, or will marked as invalid.
         */
        static std::istream& get_obj(std::istream& in, BAM& obj)
        {
            if (obj.has_data_)
                obj.brief_reset();            
            // read alignment from file
            char* data = new char[sizeof(std::int32_t)];
            if (!obj.header_.read_byte_data(
                    in, data, sizeof(std::int32_t)))
            {
                delete[] data;
                return in;
            }
            std::int32_t block_size = 
                Header::convert_char<std::int32_t>(data);
            delete[] data;
            data = new char[block_size];
            obj.header_.read_byte_data(in, data, block_size);
            std::int32_t data_counter = 32;  // to TLEN field

            // directly parse simple data into data_members_
            std::get<MEMBER_INDEX::REF_ID>(obj.data_members_) = 
                Header::convert_char<std::int32_t>(&data[0]);
            std::get<MEMBER_INDEX::POS>(obj.data_members_) = 
                Header::convert_char<std::int32_t>(&data[4]);
            std::get<MEMBER_INDEX::MAPQ>(obj.data_members_) = 
                Header::convert_char<std::int32_t>(&data[9]);
            std::get<MEMBER_INDEX::BIN>(obj.data_members_) = 
                Header::convert_char<std::int32_t>(&data[10]);
            std::get<MEMBER_INDEX::FLAG>(obj.data_members_) = 
                Header::convert_char<std::uint16_t>(&data[14]);
            std::get<MEMBER_INDEX::RNEXT_ID>(obj.data_members_) = 
                Header::convert_char<std::int32_t>(&data[20]);
            std::get<MEMBER_INDEX::PNEXT>(obj.data_members_) = 
                Header::convert_char<std::int32_t>(&data[24]);
            std::get<MEMBER_INDEX::TLEN>(obj.data_members_) = 
                Header::convert_char<std::int32_t>(&data[28]);
            // read_name 8
            std::uint8_t l_read_name = 
                Header::convert_char<std::uint8_t>(&data[8]);
            std::get<MEMBER_INDEX::QNAME>(obj.data_members_) = 
                &data[data_counter];
            data_counter += l_read_name;
            // CIGAR
            std::uint16_t n_cigar_op = 
                Header::convert_char<std::uint16_t>(&data[12]);
            int2cigar(obj, &data[data_counter], n_cigar_op);
            data_counter += n_cigar_op * sizeof(std::uint32_t);
            // SEQ 16
            std::uint8_t int_seq;
            std::int32_t l_seq = Header::convert_char<std::int32_t>(&data[16]);
            LightString seq;
            for (std::size_t i = 0;i < (l_seq + 1) / 2;++i)
            {
                int_seq = 
                    Header::convert_char<std::uint8_t>(&data[data_counter]);
                seq.append(SEQ_NUM_TO_CHAR[int_seq >> 4]);
                seq.append(SEQ_NUM_TO_CHAR[int_seq & 0xf]);
                data_counter += sizeof(std::uint8_t);
            }
            if (seq.c_str[seq.size - 1] == '=')
                --seq.size;
            std::get<MEMBER_INDEX::SEQ>(obj.data_members_) = 
                std::string(seq.c_str, seq.size);
            // QUAL
            char* qual_string = new char[l_seq];
            std::memcpy(qual_string, &data[data_counter], l_seq);
            std::get<MEMBER_INDEX::QUAL>(obj.data_members_) = 
                std::string(qual_string, l_seq);
            data_counter += l_seq;
            // optional fields
            std::vector<OptionalFieldType>& of = 
                std::get<MEMBER_INDEX::OPTIONAL_FIELDS>(obj.data_members_);
            std::array<char, 2> tag;
            char val_type;
            while (data_counter != block_size)
            {
                tag[0] = data[data_counter++];
                tag[1] = data[data_counter++];
                val_type = data[data_counter++];
                switch (val_type)
                {
                    case 'A':
                    case 'c':
                    case 'C':
                    case 's':
                    case 'S':
                    case 'i':
                    case 'I':
                    case 'f':
                    {
                        std::uint8_t size = type2Size(val_type);
                        of.emplace_back(tag
                                      , val_type
                                      , std::string(&data[data_counter]
                                                  , size));
                        data_counter += size;
                        break;
                    }
                    case 'Z':
                    case 'H':
                    {
                        std::int32_t beg = data_counter;
                        while (data[data_counter++] != 0) ;
                        of.emplace_back(tag
                                      , val_type
                                      , std::string(&data[beg]
                                                  , data_counter - beg));
                        break;
                    }
                    case 'B':
                    {
                        std::uint8_t size = type2Size(data[data_counter]);
                        std::int32_t count = 
                            Header::convert_char<std::int32_t>(
                                &data[data_counter + 1]);
                        of.emplace_back(tag
                                      , val_type
                                      , std::string(&data[data_counter]
                                                  , size * count + 5));
                        data_counter += size * count + 5;
                        if (tag[0] == 'C' && 
                            tag[1] == 'G' && 
                            size == sizeof(std::uint32_t))
                                int2cigar(obj
                                        , &data[data_counter - size * count]
                                        , count);
                        break;
                    }
                }
            }
            obj.has_data_ = true;
            delete[] data;
            delete[] qual_string;
            return in;
        };

        /**
         * @brief Read alignment information from BAM file with 
         *        random access.
         * @param in istream which contains the BAM file
         * @param obj The object which want data tpo put into
         * @param bai The index information is related to the BAM file
         * @return Parameter in for continuous data flow
         * @warning Parameter in should be the istream object which 
         *          has been used to parse into the header object, 
         *          which is referenced by parameter obj.<BR>
         *          Parameter bai should contains the index 
         *          information related with corresponding BAM file, 
         *          or will lead to undefined behaviors.
         * @see get_obj(std::istream&, BAM&)
         *
         * Read an alignment information in the specific reference 
         * region which is in the BAM file.<BR>
         * First function will check index information is usable 
         * then start to find the alignment which is in the region. 
         * After found alignment and all field are loaded sucessfully, 
         * this object will become valid state, or will marked as invalid.
         */
        static std::istream& get_obj(std::istream& in, BAM& obj, BAI& bai)
        {
            if (bai.overlap_chunks_.empty())
                obj.brief_reset();
            else if (!bai.is_useable())
            {
                std::cerr << "ERROR: Index information is empty";
                std::cerr << std::endl;
            }
            else
            {
                // start to jump to first overlap chunk
                if (bai.current_index_ == -1)
                {
                    bai.current_index_ = 0;
                    obj.header_.seek(in
                                   , std::get<BAI::CHUNK_INDEX::BEG>(
                                        bai.overlap_chunks_[0]));
                }
                // file position is out of our traking, needs relocate
                else if (obj.header_.tell() > 
                            std::get<BAI::CHUNK_INDEX::END>(
                                bai.overlap_chunks_[bai.current_index_]))
                {
                    std::uint64_t current_pos = obj.header_.tell();
                    auto beg_it = 
                        bai.overlap_chunks_.begin() + ++bai.current_index_;
                    auto end_it = bai.overlap_chunks_.end();
                    auto step = std::distance(beg_it, end_it) / 2;
                    do
                    {
                        if (current_pos > 
                                std::get<BAI::CHUNK_INDEX::BEG>(
                                    *(beg_it + step)))
                        {
                            bai.current_index_ += step;
                            std::advance(beg_it, step);
                        }
                        else if (current_pos < 
                                    std::get<BAI::CHUNK_INDEX::BEG>(
                                        *(beg_it + step)))
                            std::advance(end_it, -step);
                        else
                            break;
                        step = std::distance(beg_it, end_it) / 2;
                    } while (step != 0);
                    if (current_pos >= 
                            std::get<BAI::CHUNK_INDEX::END>(*beg_it))
                    {
                        ++bai.current_index_;
                        obj.header_.seek(in
                                       , std::get<BAI::CHUNK_INDEX::BEG>(
                                            *end_it));
                    }
                }

                ALIGNMENT_STATUS current_state;
                do
                {
                    if (obj.header_.tell() >= 
                            std::get<BAI::CHUNK_INDEX::END>(
                                bai.overlap_chunks_[bai.current_index_])
                        && 
                        bai.current_index_ != 
                            bai.overlap_chunks_.size() - 1)
                    {
                        obj.header_.seek(in
                                       , std::get<BAI::CHUNK_INDEX::BEG>(
                                            bai.overlap_chunks_
                                            [++bai.current_index_]));
                    }
                    get_obj(in, obj);
                    current_state = check_alignment_status(obj, bai);
                } while (current_state == ALIGNMENT_STATUS::NO_OVERLAP);
                if (current_state == ALIGNMENT_STATUS::OUT_RANGE)
                    obj.brief_reset();
            }
            return in;
        };

        /**
         * @brief Convert this alignment information to 
         *        BAM binary format.
         * @return Converted binary data stream stored in std::string
         * @warning If object is invalid will return empty string.
         * @see is_valid()
         */
        std::string get_bamdata() const
        {
            char data[Header::CHUNK_SIZE + 4];
            std::uint32_t size = 0;
            if (!has_data_)
                return std::string();
            
            size = 4;
            Header::parse_char<std::int32_t>(&data[size]
                                      , std::get<MEMBER_INDEX::REF_ID>(
                                            data_members_)
                                      , &size);
            Header::parse_char<std::int32_t>(&data[size]
                                      , std::get<MEMBER_INDEX::POS>(
                                            data_members_)
                                      , &size);
            const std::string& read_name = 
                std::get<MEMBER_INDEX::QNAME>(data_members_);
            Header::parse_char<std::uint8_t>(&data[size]
                                      , read_name.size() + 1
                                      , &size);
            Header::parse_char<std::uint8_t>(&data[size]
                                      , std::get<MEMBER_INDEX::MAPQ>(
                                            data_members_)
                                      , &size);
            Header::parse_char<std::uint16_t>(&data[size]
                                       , std::get<MEMBER_INDEX::BIN>(
                                            data_members_)
                                       , &size);
            const std::vector<CigarType>& cigar = 
                std::get<MEMBER_INDEX::CIGAR>(data_members_);
            Header::parse_char<std::uint16_t>(&data[size]
                                       , cigar.size()
                                       , &size);
            Header::parse_char<std::uint16_t>(&data[size]
                                       , std::get<MEMBER_INDEX::FLAG>(
                                            data_members_)
                                       , &size);
            const std::string& seq = 
                std::get<MEMBER_INDEX::SEQ>(data_members_);
            Header::parse_char<std::int32_t>(&data[size]
                                      , seq.size()
                                      , &size);
            Header::parse_char<std::int32_t>(&data[size]
                                      , std::get<MEMBER_INDEX::RNEXT_ID>(
                                            data_members_)
                                      , &size);
            Header::parse_char<std::int32_t>(&data[size]
                                      , std::get<MEMBER_INDEX::PNEXT>(
                                            data_members_)
                                      , &size);
            Header::parse_char<std::int32_t>(&data[size]
                                      , std::get<MEMBER_INDEX::TLEN>(
                                            data_members_)
                                      , &size);
            // read_name
            std::memcpy(&data[size], read_name.c_str(), read_name.size() + 1);
            size += read_name.size() + 1;
            // cigar
            std::uint8_t cigar_op;
            std::uint32_t cigar_len;
            for (std::size_t i = 0;i < cigar.size();++i)
            {
                cigar_op = 
                    get_index(CIGAR_NUM_TO_CHAR, sizeof(CIGAR_NUM_TO_CHAR)
                            , std::get<bam::CIGAR_INDEX::OP>(cigar[i]));
                cigar_len = std::get<bam::CIGAR_INDEX::LENGTH>(cigar[i]);
                Header::parse_char<std::uint32_t>(&data[size]
                                           , cigar_len << 4 | cigar_op
                                           , &size);
            }
            // seq
            std::uint8_t seq1, seq2;
            for (std::size_t i = 0;i < seq.size();i += 2)
            {
                seq1 = get_index(SEQ_NUM_TO_CHAR, 
                                    sizeof(SEQ_NUM_TO_CHAR), seq[i]);
                seq2 = (i + 1 == seq.size()) ? 
                        0 : get_index(SEQ_NUM_TO_CHAR, 
                                        sizeof(SEQ_NUM_TO_CHAR), seq[i + 1]);
                Header::parse_char<std::uint8_t>(&data[size]
                                          , seq1 << 4 | seq2
                                          , &size);
            }
            // qual
            const std::string& qual = 
                std::get<MEMBER_INDEX::QUAL>(data_members_);
            std::memcpy(&data[size], qual.c_str(), qual.size());
            size += qual.size();
            // optional fields
            const std::vector<OptionalFieldType>& of = 
                std::get<MEMBER_INDEX::OPTIONAL_FIELDS>(data_members_);
            for (std::size_t i = 0;i < of.size();++i)
            {
                data[size++] = 
                    std::get<OPTIONAL_FIELD_INDEX::TAG>(of[i])[0];
                data[size++] = 
                    std::get<OPTIONAL_FIELD_INDEX::TAG>(of[i])[1];
                data[size++] = 
                    std::get<OPTIONAL_FIELD_INDEX::VALUE_TYPE>(of[i]);
                const std::string& str = 
                    std::get<OPTIONAL_FIELD_INDEX::VALUE>(of[i]);
                std::memcpy(&data[size], str.c_str(), str.size());
                size += str.size();
            }
            Header::parse_char<std::int32_t>(data, size - sizeof(std::int32_t));
            return std::string(data, size);
        }

        /**
         * @brief Convert list of BAM object to 
         *        BGZF compressed format and output.
         * @param out ostream which data will output to
         * @param obj List of BAM object to convert
         * @see get_bamdata()
         *
         * Convert list of BAM object to BGZF compressed format 
         * and output to the istream.<BR>
         * First convert each object to binary format then 
         * compress these data to BGZF file format.
         */
        static void dump(std::ostream& out, std::vector<BAM>& obj)
        {
            if (!obj.empty())
                out << obj[0].header_;

            char data[Header::CHUNK_SIZE];
            std::uint32_t size = 0;
            std::string buf;
            for (std::size_t i = 0;i < obj.size();++i)
            {
                buf = obj[i].get_bamdata();
                if (size + buf.size() > Header::CHUNK_SIZE)
                {
                    Header::deflate_block(out, data, size);
                    size = 0;
                }
                std::memcpy(&data[size], buf.c_str(), buf.size());
                size += buf.size();
            }
            if (size != 0)
                Header::deflate_block(out, data, size);
            out.write(EOF_marker, sizeof(EOF_marker));
        };

        /// default copy assignment
        BAM& operator=(const BAM& rhs)
        {
            has_data_ = rhs.has_data_;
            header_ = rhs.header_;
            data_members_ = rhs.data_members_;
            return *this;
        }

        /// default move assignment
        BAM& operator=(BAM&& rhs)
        {
            has_data_ = rhs.has_data_;
            header_ = rhs.header_;
            data_members_ = std::move(rhs.data_members_);
            return *this;
        }

        /**
         * @brief overload bool() operator.
         * @return Whether alignment information is valid
         */
        operator bool() const
        {
            return has_data_;
        }

        /**
         * @brief overload >> operator.
         * @param in istream which data comes from
         * @param rhs BAM object which will contain alignment information
         * @return Parameter in for continuous data flow
         * @see get_obj(std::istream&, BAM&)
         * 
         * Overload >> operator for user to use more intuitively.
         * Same as call get_obj() with sequencial access.
         */
        friend std::istream& operator>> (std::istream& in, BAM& rhs)
        {
            return get_obj(in, rhs);
        };

        /**
         * @brief overload << operator.
         * @param out ostream which data will output to
         * @param rhs BAM object which want to output
         * @return Parameter out for continuous data flow
         * @see get_bamdata()
         *
         * Output an alignment information in BGZF compressed format.<BR>
         * Same as call get_bamdata() and compress it to a block.
         */
        friend std::ostream& operator<< 
            (std::ostream& out, const BAM& rhs)
        {
            std::string buf = rhs.get_bamdata();
            Header::deflate_block(out
                                , const_cast<char*>(buf.c_str())
                                , buf.size());
            return out;
        };

        /// BGZF file format's EOF marker
        constexpr static char EOF_marker[] = 
            "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00";
      
      private:

        /// Only clear the field that is necessary
        void brief_reset()
        {
            has_data_ = false;
            std::get<MEMBER_INDEX::OPTIONAL_FIELDS>(data_members_)
                .resize(0);
        }

        /**
         * @brief Calculate value type to their corresponding size.
         * @param c The value type
         * @return Size of value type
         */
        constexpr static std::uint8_t type2Size(char c)
        {
            switch (c)
            {
                case 'A':
                case 'c':
                case 'C':
                    return 1; // sizeof(char, int8_t, std::uint8_t)
                case 's':
                case 'S':
                    return 2; // sizeof(int16_t, uint16_t);
                case 'i':
                case 'I':
                case 'f':
                    return 4; // sizeof(int32_t, uint32_t, float_t);
                default:
                    return 0;
            }
        }

        /**
         * @brief Get target's index in character array.
         * @param arr The character array.
         * @param target Target element want to search.
         * @return If found, return index. Else return -1
         */
        std::uint8_t get_index(const char* arr, std::size_t size, char target) const
        {
            return std::distance(arr
                               , std::find(arr
                                         , arr + size
                                         , target));
        }

        /**
         * @brief Convert Cigar in BAM compressed format to CigarType.
         * @param obj BAM object which CigarType will put into
         * @param data Data stream pointer points to 
         *             the start of Cigar data
         * @param count Number of Cigar operations
         */
        static void int2cigar(BAM& obj, char* data, std::int32_t count)
        {
            std::vector<CigarType>& cigar = 
                std::get<MEMBER_INDEX::CIGAR>(obj.data_members_);
            cigar.clear();
            std::uint32_t int_cigar;
            const std::uint32_t last_element = count * sizeof(std::uint32_t);
            for (std::size_t data_counter = 0; 
                 data_counter < last_element; 
                 data_counter += sizeof(std::uint32_t))
            {
                int_cigar = 
                    Header::convert_char<std::uint32_t>(&data[data_counter]);
                cigar.emplace_back(CIGAR_NUM_TO_CHAR[int_cigar & 0xf]
                                 , int_cigar >> 4);
            }
        }

        /// Enumerating different overlapping state 
        /// between alignment and region.
        enum class ALIGNMENT_STATUS
        {
            NO_OVERLAP
          , OVERLAPPED
          , OUT_RANGE
        };

        /**
         * @brief Determine whether BAM alignment is 
         *        overlapping with region.
         * @param obj BAM object which want to test
         * @param bai BAI object contains the interested region
         * @return Overlapping status between alignment and region
         */
        static ALIGNMENT_STATUS check_alignment_status(BAM& obj, BAI& bai)
        {
            if (!obj.has_data_)
                return ALIGNMENT_STATUS::OUT_RANGE;
            std::int32_t bam_ref = 
                std::get<MEMBER_INDEX::REF_ID>(obj.data_members_);
            std::int32_t bam_start = 
                std::get<MEMBER_INDEX::POS>(obj.data_members_);
            std::int32_t bam_end = 
                bam_start + std::get<MEMBER_INDEX::QUAL>(
                    obj.data_members_).size();
            RegionType region = bai.get_region();

            if (bam_ref == -1)
                return ALIGNMENT_STATUS::OUT_RANGE;
            else if (!bai.is_start_region_specified())
                return ALIGNMENT_STATUS::OVERLAPPED;
            else if (bam_ref < 
                        std::get<REGION_INDEX::LEFT_REF_ID>(region))
                return ALIGNMENT_STATUS::NO_OVERLAP;
            else if (bam_ref == 
                        std::get<REGION_INDEX::LEFT_REF_ID>(region))
            {
                if (bam_end < 
                        std::get<REGION_INDEX::LEFT_POS>(region))
                    return ALIGNMENT_STATUS::NO_OVERLAP;
                else if (bai.is_end_region_specified()
                         &&
                         bam_start > 
                            std::get<REGION_INDEX::RIGHT_POS>(region)
                         &&
                         bam_ref == 
                            std::get<REGION_INDEX::RIGHT_REF_ID>(region))
                    return ALIGNMENT_STATUS::OUT_RANGE;
                else
                    return ALIGNMENT_STATUS::OVERLAPPED;
            }
            else
            {
                if (bai.is_end_region_specified())
                {
                    if (bam_ref > 
                            std::get<REGION_INDEX::RIGHT_REF_ID>(region))
                        return ALIGNMENT_STATUS::OUT_RANGE;
                    else if (bam_ref < 
                                std::get<REGION_INDEX::RIGHT_REF_ID>(region))
                        return ALIGNMENT_STATUS::OVERLAPPED;
                    else
                    {
                        if (bam_start > 
                                std::get<REGION_INDEX::RIGHT_POS>(region))
                            return ALIGNMENT_STATUS::OUT_RANGE;
                        else
                            return ALIGNMENT_STATUS::OVERLAPPED;
                    }
                }
                else
                    return ALIGNMENT_STATUS::OVERLAPPED;
            }
        }

        /// Light string for internal use.
        struct LightString
        {
            LightString()
            : c_str     ( new char[1024] )
            , size      (       0        )
            , capacity  (      1024      )
            {
            }

            LightString (const LightString& rhs)
            {
                size = rhs.size;
                capacity = rhs.capacity;
                c_str = new char[capacity];
                std::memcpy(c_str, rhs.c_str, size);
            }

            LightString (LightString&& rhs) = default;

            ~LightString()
            {
                delete[] c_str;
            }

            void extend()
            {
                char* tmp = new char[2 * capacity];
                std::memcpy(tmp, c_str, size);
                delete[] c_str;
                c_str = tmp;
                capacity = 2 * capacity;
            }

            void append(const char* c, std::size_t s)
            {
                while (size + s > capacity)
                    extend();
                std::memcpy(&c_str[size], c, s);
                size += s;
            }

            void append(const std::string& str)
            {
                append(str.c_str(), str.size());
            }

            void append(const char c)
            {
                if (size + 1 > capacity)
                    extend();
                c_str[size++] = c;
            }

            char* c_str;
            std::size_t size;
            std::size_t capacity;
        };

        /**
         * @brief According to value type to read value from data stream.
         * @param result String that result need to append to
         * @param c The value type
         * @param pos Data stream pointer points to the start of value
         * @param offset When value is an array, 
         *               needs offset to help to locate right position
         */
        static void value_type_to_string
            (LightString& result, char c, const char* pos, std::size_t offset = 0)
        {
            std::size_t size;
            char tmp[45];
            switch (c)
            {
                case 'c':
                    size = std::sprintf(tmp, "%d", 
                                    Header::convert_char<std::int8_t>(
                                        pos + offset * sizeof(std::int8_t)));
                    break;
                case 'C':
                    size = std::sprintf(tmp, "%u", 
                                    Header::convert_char<std::uint8_t>(
                                        pos + offset * sizeof(std::uint8_t)));
                    break;
                case 's':
                    size = std::sprintf(tmp, "%d", 
                                    Header::convert_char<std::int16_t>(
                                        pos + offset * sizeof(std::int16_t)));
                    break;
                case 'S':
                    size = std::sprintf(tmp, "%u", 
                                    Header::convert_char<std::uint16_t>(
                                        pos + offset * sizeof(std::uint16_t)));
                    break;
                case 'i':
                    size = std::sprintf(tmp, "%d", 
                                    Header::convert_char<std::int32_t>(
                                        pos + offset * sizeof(std::int32_t)));
                    break;
                case 'I':
                    size = std::sprintf(tmp, "%u", 
                                    Header::convert_char<std::uint32_t>(
                                        pos + offset * sizeof(std::uint32_t)));
                    break;
                case 'f':
                    size = std::sprintf(tmp, "%f", 
                                    Header::convert_char<float>(
                                        pos + offset * sizeof(float)));
                    break;
                default:
                    assert(false && "Invalid value type in optional field, file might broken.");
            }
            result.append(tmp, size);
        }

        /**
         * @brief Convert list of CigarType to 
         *        Cigar string in SAM text format.
         * @param result String that result need to append to
         * @param vec List of CigarType to convert
         */
        static void cigar_vec2str
            (LightString& result, const std::vector<CigarType>& vec)
        {
            if (vec.size() == 0)
                result.append('*');
            std::size_t size;
            char tmp[32];
            for (std::size_t i = 0;i < vec.size();++i)
            {
                size = std::sprintf(tmp, "%u", std::get<bam::CIGAR_INDEX::LENGTH>(vec[i]));
                result.append(tmp, size);
                result.append(std::get<bam::CIGAR_INDEX::OP>(vec[i]));
            }
        }

        /**
         * @brief Convert Cigar string to list of CigarType.
         * @param cigar_str Target string to convert
         * @param consume_ref_len Reference counter for bin calculation
         * @return Result in vector of CigarType
         */
        static std::vector<CigarType> str2cigar_vec
            (std::string& cigar_str, std::uint32_t& consume_ref_len)
        {
            char cigar_op;
            std::uint32_t cigar_len, base;
            std::int64_t counter = cigar_str.size() - 1;
            std::vector<CigarType> result;
            result.reserve(5);
            if (cigar_str != "*")
            {
                while (counter != -1)
                {
                    cigar_op = cigar_str[counter--];
                    cigar_len = 0;
                    base = 1;
                    while (counter != -1 && 
                            cigar_str[counter] >= '0' && cigar_str[counter] <= '9')
                    {
                        cigar_len += (cigar_str[counter--] - '0') * base;
                        base *= 10;
                    }
                    if (std::find(CONSUME_REF_CIGAR
                                , CONSUME_REF_CIGAR + sizeof(CONSUME_REF_CIGAR)
                                , cigar_op) !=
                            CONSUME_REF_CIGAR + sizeof(CONSUME_REF_CIGAR))
                        consume_ref_len += cigar_len;
                    result.emplace(result.begin(), cigar_op, cigar_len);
                }
            }
            return result;
        }

        /**
         * @brief Convert quality in BAM format to 
         *        SAM visable text format.
         * @param result String that result need to append to
         * @param qual quality in BAM format
         */
        static void qual2ascii(LightString& result, const std::string& qual)
        {
            if (qual[0] == (char)255)
            {
                result.append('*');
                return;
            }

            std::size_t i = result.size;
            result.append(qual);
            for ( ;i < result.size;++i)
                result.c_str[i] += 33;
        }

        /**
         * @brief Convert value in BAM binary format to SAM text format.
         * @param tmp String that result need to append to
         * @param value_type The value type of value
         * @param target Target binary format to convert
         * @see value_type_to_string()
         */
        static void unpack_optional_field_string
            (LightString& tmp, char value_type, const std::string& target)
        {
            switch (value_type)
            {
                case 'A':
                    tmp.append(target.c_str(), target.size());
                    break;
                case 'Z':
                    tmp.append(target.c_str(), target.size() - 1);
                    break;
                case 'H':
                {
                    const char* data_str = target.c_str();
                    const char hex_char[] = "0123456789ABCDEF";
                    for (std::size_t j = 1;data_str[j] != '\0';++j)
                    {
                        tmp.append(hex_char[data_str[j] >> 4]);
                        tmp.append(hex_char[data_str[j] & 0xf]);
                    }
                    break;
                }
                case 'B':
                {
                    const char* data_str = target.c_str();
                    tmp.append(data_str[0]);
                    tmp.append(',');
                    std::int32_t count = 
                        bam::Header::convert_char<std::int32_t>(&data_str[1]);
                    for (std::size_t j = 0;j < count;++j)
                    {
                        value_type_to_string(tmp, data_str[0]
                                            , &data_str[5], j);
                        tmp.append(",", (j == count - 1 ? 0 : 1));
                    }
                    break;
                }
                default:
                    value_type_to_string(tmp, value_type, target.c_str());
            }
        }

        /**
         * @brief Convert value in SAM text format to BAM binary format.
         * @param c The value type
         * @param pos Data stream pointer points to 
         *            the next byte to write
         * @param target Target value needs to convert
         * @param offset (Optional) Data stream index remainder, 
         *               will pass to parse_char(). Default is nullptr
         * @see Header::parse_char()
         */
        void string_to_value_type(
            char c
          , char* pos
          , std::int64_t target
          , std::uint32_t* offset = nullptr
        )
        {
            switch (c)
            {
                case 'c':
                    Header::parse_char<std::int8_t>(pos, target, offset);
                    return;
                case 'C':
                    Header::parse_char<std::uint8_t>(pos, target, offset);
                    return;
                case 's':
                    Header::parse_char<std::int16_t>(pos, target, offset);
                    return;
                case 'S':
                    Header::parse_char<std::uint16_t>(pos, target, offset);
                    return;
                case 'i':
                    Header::parse_char<std::int32_t>(pos, target, offset);
                    return;
                case 'I':
                    Header::parse_char<std::uint32_t>(pos, target, offset);
                    return;
            }
        }

        /**
         * @brief Convert value in SAM text format to BAM binary format.
         * @param value_type The value type of value
         * @param target Target SAM text format to convert
         * @param converted_int Preparsed integer value
         * @return String of value in BAM binary format
         * @see string_to_value_type()
         */
        std::string pack_optional_field(
            char value_type
          , const std::string& target
          , std::int64_t converted_int
        )
        {
            char buffer[Header::CHUNK_SIZE];
            std::uint32_t buffer_size = 0;
            switch (value_type)
            {
                case 'A':
                    buffer[buffer_size++] = target[0];
                    break;
                case 'f':
                {
                    float f;
                    boost::spirit::qi::parse(target.begin()
                                           , target.end()
                                           , boost::spirit::float_
                                           , f);
                    Header::parse_char<float>(buffer
                                            , f
                                            , &buffer_size);
                    break;
                }
                case 'Z':
                case 'H':
                    std::memcpy(&buffer, target.c_str(), target.size() + 1);
                    buffer_size += target.size() + 1;
                    break;
                case 'B':
                {
                    std::vector<std::string> split_vec;
                    boost:split(split_vec
                              , target
                              , boost::is_any_of(", ")
                              , boost::token_compress_on);
                    buffer[buffer_size++] = split_vec[0][0];
                    std::int32_t vec_size = split_vec.size();
                    Header::parse_char<std::int32_t>(&buffer[buffer_size]
                                              , vec_size - 1
                                              , &buffer_size);
                    for (std::size_t i = 1;i < vec_size;++i)
                    {
                        if (buffer[0] != 'f')
                            boost::spirit::qi::parse(split_vec[i].begin()
                                                   , split_vec[i].end()
                                                   , boost::spirit::long_
                                                   , converted_int);
                        switch (buffer[0])
                        {
                            case 'f':
                            {
                                float f;
                                boost::spirit::qi::parse(split_vec[i].begin()
                                                       , split_vec[i].end()
                                                       , boost::spirit::float_
                                                       , f);
                                Header::parse_char<float>(&buffer[buffer_size]
                                                        , f
                                                        , &buffer_size);
                                break;
                            }
                            default:
                                string_to_value_type(buffer[0]
                                                   , &buffer[buffer_size]
                                                   , converted_int
                                                   , &buffer_size);
                        }
                    }
                }
                default:
                    string_to_value_type(value_type
                                       , buffer
                                       , converted_int
                                       , &buffer_size);
            }
            return std::string(buffer, buffer_size);
        }

        /**
         * @brief Check whether value type needs to convert to 'i' type.
         * @param value_type Original value type
         * @return Converted value type
         *
         * This function might be used when BAM convert to SAM, 
         * because some value type in BAM is belongs to 'i' type in SAM.
         */
        static char check_to_i_type(char value_type)
        {
            if (std::find(VALUE_TYPE_TO_I
                        , VALUE_TYPE_TO_I + sizeof(VALUE_TYPE_TO_I)
                        , value_type) !=
                    VALUE_TYPE_TO_I + sizeof(VALUE_TYPE_TO_I))
                return 'i';
            else
                return value_type;
        }

        /**
         * @brief Convert optional fields to SAM text format.
         * @param tmp String that result need to append to
         * @param target List of OptionalFieldType needs to convert
         */
        void optional_fields_to_string
            (LightString& tmp, std::vector<OptionalFieldType>& target)
        {
            if (target.empty())
                return;

            char value_type;
            for (std::size_t i = 0;i < target.size();++i)
            {
                tmp.append( 
                    std::get<OPTIONAL_FIELD_INDEX::TAG>(target[i])[0]);
                tmp.append(
                    std::get<OPTIONAL_FIELD_INDEX::TAG>(target[i])[1]);
                tmp.append(':');
                value_type = 
                    std::get<OPTIONAL_FIELD_INDEX::VALUE_TYPE>(target[i]);
                tmp.append(check_to_i_type(value_type));
                tmp.append(':');
                unpack_optional_field_string(tmp, value_type
                                            , std::get<
                                                OPTIONAL_FIELD_INDEX::VALUE
                                                >(target[i]));
                tmp.append('\t');
            }
            tmp.size -= 1;
        }

        /**
         * @brief Implememtation details of to_string.
         * @tparam Current field index needs to convert
         * @param result String that result need to append to
         */
        template <std::size_t N>
        void to_string_impl(LightString& result)
        {
            std::tuple_element_t<N, MemberType>& target = 
                std::get<N>(data_members_);
            if constexpr (N == MEMBER_INDEX::OPTIONAL_FIELDS)
            {
                optional_fields_to_string(result, target);
                if (result.c_str[result.size - 1] == '\t')
                    --result.size;
            }
            else
            {
                if constexpr (N == MEMBER_INDEX::REF_ID)
                {
                    if (target == -1)
                        result.append('*');
                    else
                        result.append(
                            std::get<sam::REFERENCE_INDEX::REFERENCE_NAME>(
                                header_.get_member<
                                    sam::HEADER_INDEX::REFERENCE>()[target]));
                }
                else if constexpr (N == MEMBER_INDEX::RNEXT_ID)
                {
                    if (target == -1)
                        result.append('*');
                    else if (target == std::get<
                                MEMBER_INDEX::REF_ID>(data_members_))
                        result.append('=');
                    else
                        result.append(std::get<
                            sam::REFERENCE_INDEX::REFERENCE_NAME>(
                                header_.get_member<
                                    sam::HEADER_INDEX::REFERENCE>()[target]));
                }
                else if constexpr (N == MEMBER_INDEX::CIGAR)
                    cigar_vec2str(result, target);
                else if constexpr (N == MEMBER_INDEX::QUAL)
                    qual2ascii(result, target);
                else if constexpr (std::is_same<
                                    decltype(target), std::string&>())
                    result.append(target);
                else
                {
                    std::size_t size;
                    char tmp[35];
                    if constexpr (N == MEMBER_INDEX::POS || 
                                  N == MEMBER_INDEX::PNEXT)
                        size = std::sprintf(tmp, "%d", target + 1);
                    else
                        size = std::sprintf(tmp, "%d", target);
                    result.append(tmp, size);
                }
                result.append('\t');
                to_string_impl<N + 1>(result);
            }
        }

        /**
         * @brief Calculate a region is belongs to which bin.
         * @param beg Starting position
         * @param end Ending position
         * @return Result bin number
         *
         * Calculate the list of bins that may overlap with region 
         * [beg,end) (zero-based).
         */
        static std::uint16_t reg2bin(std::int32_t beg, std::int32_t end)
        {
            --end;
            if (beg >> 14 == end >> 14) 
                return ((1 << 15) - 1) / 7 + (beg >> 14);
            if (beg >> 17 == end >> 17) 
                return ((1 << 12) - 1) / 7 + (beg >> 17);
            if (beg >> 20 == end >> 20) 
                return ((1 <<  9) - 1) / 7 + (beg >> 20);
            if (beg >> 23 == end >> 23) 
                return ((1 <<  6) - 1) / 7 + (beg >> 23);
            if (beg >> 26 == end >> 26) 
                return ((1 <<  3) - 1) / 7 + (beg >> 26);
            return 0;
        }

        /// Record maximum number of Cigar operations 
        /// that BAM could support.
        const static std::uint16_t MAX_CIGAR_OP_NUM = 65535;
        /// Record which value type in BAM should convert to 'i' in SAM.
        constexpr static char VALUE_TYPE_TO_I[] = 
        {
            'c', 'C', 's', 'S', 'i', 'I'
        };
        /// Record which Cigar operation will consumes reference position.
        constexpr static char CONSUME_REF_CIGAR[] = 
        {
            'M', 'D', 'N', '=', 'X'
        };
        /// A lookup table from Cigar operation code to character.
        constexpr static char CIGAR_NUM_TO_CHAR[] = 
        {
            'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'
        };
        /// A lookup table from sequence code to character.
        constexpr static char SEQ_NUM_TO_CHAR[] = 
        {
            '=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 
            'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
        };       

        /// Record whether current alignment information is valid.
        bool has_data_;
        /// Reference to corresponding header.
        Header& header_;
        /// Store current alignment information.
        MemberType data_members_;
    };

    SAM::SAM(const BAM& rhs)
    : header_   ( rhs.header_ )
    {
        BAM::LightString ls;
        std::get<sam::MEMBER_INDEX::QNAME>(data_members_) = 
            std::get<MEMBER_INDEX::QNAME>(rhs.data_members_);
        std::get<sam::MEMBER_INDEX::FLAG>(data_members_) = 
            std::get<MEMBER_INDEX::FLAG>(rhs.data_members_);
        // RNAME
        std::int32_t ref_id = 
            std::get<MEMBER_INDEX::REF_ID>(rhs.data_members_);
        std::vector<sam::ReferenceType>& refs = 
            std::get<sam::HEADER_INDEX::REFERENCE>(header_.header_);
        std::get<sam::MEMBER_INDEX::RNAME>(data_members_) = 
            (ref_id == -1) ? "*" : 
                std::get<sam::REFERENCE_INDEX::REFERENCE_NAME>(refs[ref_id]);
        std::get<sam::MEMBER_INDEX::POS>(data_members_) = 
            std::get<MEMBER_INDEX::POS>(rhs.data_members_) + 1;
        std::get<sam::MEMBER_INDEX::MAPQ>(data_members_) = 
            std::get<MEMBER_INDEX::MAPQ>(rhs.data_members_);
        // CIGAR
        BAM::cigar_vec2str(ls, 
            std::get<MEMBER_INDEX::CIGAR>(rhs.data_members_));
        std::get<sam::MEMBER_INDEX::CIGAR>(data_members_) = 
            std::string(ls.c_str, ls.size);
        // RNEXT
        ref_id = std::get<MEMBER_INDEX::RNEXT_ID>(rhs.data_members_);
        std::get<sam::MEMBER_INDEX::RNEXT>(data_members_) = 
            (ref_id == -1) ? "*" : 
                (ref_id == std::get<MEMBER_INDEX::REF_ID>(
                    rhs.data_members_)) ? "=" : 
                        std::get<sam::REFERENCE_INDEX::REFERENCE_NAME>(
                            refs[ref_id]);
        std::get<sam::MEMBER_INDEX::PNEXT>(data_members_) = 
            std::get<MEMBER_INDEX::PNEXT>(rhs.data_members_) + 1;
        std::get<sam::MEMBER_INDEX::TLEN>(data_members_) = 
            std::get<MEMBER_INDEX::TLEN>(rhs.data_members_);
        std::get<sam::MEMBER_INDEX::SEQ>(data_members_) = 
            std::get<MEMBER_INDEX::SEQ>(rhs.data_members_);
        // QUAL
        ls.size = 0;
        BAM::qual2ascii(ls, 
            std::get<MEMBER_INDEX::QUAL>(rhs.data_members_));
        std::get<sam::MEMBER_INDEX::QUAL>(data_members_) = 
            std::string(ls.c_str, ls.size);
        // Optional fields
        std::vector<OptionalFieldType>& sam_of = 
            std::get<sam::MEMBER_INDEX::OPTIONAL_FIELDS>(data_members_);
        const std::vector<OptionalFieldType>& bam_of = 
            std::get<MEMBER_INDEX::OPTIONAL_FIELDS>(rhs.data_members_);
        char value_type;
        for (std::size_t i = 0;i < bam_of.size();++i)
        {
            value_type = 
                std::get<OPTIONAL_FIELD_INDEX::VALUE_TYPE>(bam_of[i]);
            ls.size = 0;
            BAM::unpack_optional_field_string(ls, value_type
                                            , std::get<
                                                OPTIONAL_FIELD_INDEX::VALUE
                                                >(bam_of[i]));
            sam_of.emplace_back(
                std::get<OPTIONAL_FIELD_INDEX::TAG>(bam_of[i])
              , BAM::check_to_i_type(value_type)
              , std::string(ls.c_str, ls.size));
        }
    }
};
