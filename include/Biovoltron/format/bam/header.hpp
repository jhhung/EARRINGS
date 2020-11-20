/**
 *  @file sam.hpp
 *  @brief The program parsing SAM's header
 *  @author JHHlab corp
 */
#include <Biovoltron/format/sam/header.hpp>

namespace biovoltron::format {
    class BAM;
    std::ostream& operator<<(std::ostream&, const biovoltron::format::BAM&);
};

namespace biovoltron::format::bam {
    /// A shortcut to SORT_TYPE in namespace sam.
    using SORT_TYPE = ::biovoltron::format::sam::SORT_TYPE;
    /// A shortcut to GROUP_TYPE in namespace sam.
    using GROUP_TYPE = ::biovoltron::format::sam::GROUP_TYPE;
    /// A shortcut to PLATFORM in namespace sam.
    using PLATFORM = ::biovoltron::format::sam::PLATFORM;
    /// A shortcut to REFERENCE_INDEX in namespace sam.
    using REFERENCE_INDEX = ::biovoltron::format::sam::REFERENCE_INDEX;
    /// A shortcut to READ_GROUP_INDEX in namespace sam.
    using READ_GROUP_INDEX = ::biovoltron::format::sam::READ_GROUP_INDEX;
    /// A shortcut to PROGRAM_INDEX in namespace sam.
    using PROGRAM_INDEX = ::biovoltron::format::sam::PROGRAM_INDEX;
    /// A shortcut to HEADER_INDEX in namespace sam.
    using HEADER_INDEX = ::biovoltron::format::sam::HEADER_INDEX;
    /// Indices of RegionType which used in BAI object.
    enum REGION_INDEX
    {
        LEFT_REF_ID, LEFT_POS, RIGHT_REF_ID, RIGHT_POS
    };
    using OPTIONAL_FIELD_INDEX = ::biovoltron::format::sam::OPTIONAL_FIELD_INDEX;
    /// Indices of CigarType stored in BAM
    enum CIGAR_INDEX
    {
        OP, LENGTH
    };
    /// Indices of BAMMemberType to access information 
    /// in a BAM alignment.
    enum MEMBER_INDEX
    {
        QNAME
      , FLAG
      , REF_ID
      , POS
      , MAPQ
      , CIGAR
      , RNEXT_ID
      , PNEXT
      , TLEN
      , SEQ
      , QUAL
      , OPTIONAL_FIELDS
      , BIN
    };

    /// A shortcut to ReferenceType in namespace sam.
    using ReferenceType = sam::ReferenceType;
    /// A shortcut to ReadGroupType in namespace sam.
    using ReadGroupType = sam::ReadGroupType;
    /// A shortcut to ProgramType in namespace sam.
    using ProgramType = sam::ProgramType;
    /// A shortcut to HeaderType in namespace sam.
    using HeaderType = sam::HeaderType;
    /// Data structure to store desired reference region in BAI object.
    using RegionType = std::tuple<
        int64_t
      , int64_t
      , int64_t
      , int64_t
    >;
    /// A shortcut to OptionalFieldType in namespace sam.
    using OptionalFieldType = sam::OptionalFieldType;
    /// Data structure to store Cigar operations in BAM alignment.
    using CigarType = std::tuple<char, uint32_t>;
    /// Data structure to store information in BAM alignment.
    using MemberType = std::tuple<
        std::string                     ///< QNAME
      , uint16_t                        ///< FLAG
      , int32_t                         ///< RNAME replaced with ID
      , int32_t                         ///< POS
      , uint8_t                         ///< MAPQ
      , std::vector<CigarType>          ///< CIGAR
      , int32_t                         ///< RNEXT replaced with ID
      , int32_t                         ///< PNEXT
      , int32_t                         ///< TLEN
      , std::string                     ///< SEQ
      , std::string                     ///< QUAL
      , std::vector<OptionalFieldType>  ///< Optional fields
      , uint16_t                        ///< BIN number in 
                                        ///< BAM binning concepts
      >;

    /**
     * @brief Class to represents BAM file format's header.
     * @warning When there has any BAM object constructed from 
     *          the header object, or there has any SAM object 
     *          contruct from the BAM object which 
     *          reference to the header object, 
     *          the header MUST NOT destruct before they destruct.
     *
     * This class stores header in BAM file format.<BR>
     * Every BAM alignments from the same BAM file MUST reference 
     * to the same Header object, or might be undefined behavior.
     */
    class Header : public sam::Header
    {
    friend class ::biovoltron::format::BAM;
      public:
        /// Default constructor to initialize every data members.
        Header()
        : block_offset_     (   0   )
        , block_length_     (   0   )
        , block_address_    (   0   )
        {
            init_zargs();
        };

        /**
         * @brief istream preparse constructor to 
         *        parse header text in istream.
         * @param in istream object which 
         *        contains header field in BAM file
         * @see Header(), preparse()
         *  
         * This is a constructor can help you to parse 
         * header information in BAM file simultaneously 
         * when you construct Header object.<BR>
         * Same as call default constructor then call preparse function.
         */
        Header(std::istream& in)
        : block_offset_     (   0   )
        , block_length_     (   0   )
        , block_address_    (   0   )
        {
            init_zargs();
            preparse(in);
        };

        /**
         * @brief Convert constructor from sam::Header to bam::Header.
         * @param rhs Source sam::Header object.
         * @see sam::Header, bam::BAM(Header&, const sam::SAM&)
         * 
         * A constructor convert sam::Header to bam::Header.<BR>
         * This constructor might be called when user want to 
         * convert SAM object to BAM object because of the 
         * header in BAM object can not directly<BR>
         * reference to the header in SAM object. So needs to create 
         * a bam::Header to help new BAM object has 
         * bam::Header to reference to.
         */
        Header(const sam::Header& rhs)
        : sam::Header       (  rhs  )
        , block_offset_     (   0   )
        , block_length_     (   0   )
        , block_address_    (   0   )
        {
            init_zargs();
        }

        Header(const Header& rhs) = default;
        Header(Header&& rhs) = default;
        Header& operator=(const Header& rhs) = default;
        Header& operator=(Header&& rhs) = default;
        ~Header() = default;

        /**
         *  @brief Clear all data in the object.
         *  @warning MUST NOT call this function when any 
         *           SAM/BAM object is referencing to this header.
         * 
         *  Clear all data in this header object.<BR>
         *  This function might be called when user want to open a 
         *  new BAM file and want to reuse this object rather than 
         *  destruct it.
         */
        void reset()
        {
            block_offset_ = 0;
            block_length_ = 0;
            block_address_ = 0;
            std::get<sam::HEADER_INDEX::VERSION>(header_).clear();
            std::get<sam::HEADER_INDEX::ALIGNMENT_SORT_ORDER>(header_) = 
                sam::SORT_TYPE::UNKNOWN;
            std::get<sam::HEADER_INDEX::ALIGNMENT_GROUPING>(header_) = 
                sam::GROUP_TYPE::NONE;
            std::get<sam::HEADER_INDEX::REFERENCE>(header_).clear();
            std::get<sam::HEADER_INDEX::READ_GROUP>(header_).clear();
            std::get<sam::HEADER_INDEX::PROGRAM>(header_).clear();
            std::get<sam::HEADER_INDEX::COMMENT>(header_).clear();
            std::get<sam::HEADER_INDEX::PLAIN_TEXT>(header_).clear();
        }

        /**
         *  @brief Parse header field in BAM file.
         *  @param in istream which contains header in BAM file.
         * 
         *  Parse header information in istream. After parsing, 
         *  header become usable.
         */
        void preparse(std::istream& in)
        {
            reset();
            // magic & l_text
            std::size_t data_capacity = 8;
            char* data = new char[8];
            read_byte_data(in, data, 8);
            if (std::memcmp(data, BAM_MAGIC, 4) != 0)
            {
                std::cerr << "ERROR: BAM Magic string not match";
                std::cerr << std::endl;
                return;
            }
            // text
            std::int32_t l_text = convert_char<std::int32_t>(&data[4]);
            if (l_text + 1 > data_capacity)
            {
                delete []data;
                data = new char[l_text + 1];
                data_capacity = l_text + 1;
            }
            read_byte_data(in, data, l_text);
            data[l_text] = '\0';
            std::string text(data);
            std::vector<std::string> lines;
            boost::split(lines, text, boost::is_any_of("\n"), 
                boost::token_compress_on);
            lines.pop_back();
            preparse_impl(lines);
            // n_ref
            read_byte_data(in, data, sizeof(std::int32_t));
            std::int32_t n_ref = convert_char<std::int32_t>(data);
            // references
            std::int32_t l_name, l_ref;
            std::string name;
            std::vector<sam::ReferenceType>& refs = 
                std::get<sam::HEADER_INDEX::REFERENCE>(header_);
            std::vector<sam::ReferenceType>::iterator it;
            for (std::size_t i = 0;i < n_ref;++i)
            {
                // name
                read_byte_data(in, data, sizeof(std::int32_t));
                l_name = convert_char<std::int32_t>(data);
                if (l_name + 1 > data_capacity)
                {
                    delete[] data;
                    data = new char[l_name + 1];
                    data_capacity = l_name + 1;
                }
                read_byte_data(in, data, l_name);
                data[l_name] = '\0';
                name = data;
                // l_ref
                read_byte_data(in, data, sizeof(std::int32_t));
                l_ref = convert_char<std::int32_t>(data);
                it = refs.begin();
                while (it != refs.end())
                {
                    if (std::get<sam::REFERENCE_INDEX::REFERENCE_NAME>(*it) 
                            == name)
                        break;
                    std::advance(it, 1);
                }
                if (it == refs.end())
                    refs.emplace_back(name, l_ref, "", "", "", "");
            }
            delete[] data;
        }

        /**
         *  @brief Overload >> operator.
         *  @param in istream which data comes from
         *  @param rhs Header object which will contains data
         *  @return Parameter in for continuous data flow
         *  @see preparse()
         * 
         *  Overload >> operator for user to use more intuitively.
         *  Same as call rhs.preparse().
         */
        friend std::istream& operator>> 
            (std::istream& in, Header& rhs)
        {
            rhs.preparse(in);
            return in;
        };

        /**
         *  @brief Overload << operator.
         *  @param out ostream which data will direct to
         *  @param rhs Header object which want to output
         *  @return Parameter out for continuous data flow
         *  @see getMember()
         * 
         *  Output header with BAM file format (binary).<BR>
         *  If want to read header in plaintext, use getMember() to
         *  get header plaintext in header data structure.
         */
        friend std::ostream& operator<< 
            (std::ostream& out, const Header& rhs)
        {
            char data[CHUNK_SIZE];
            std::uint32_t offset = 0;
            // magic string
            std::memcpy(data, BAM_MAGIC, sizeof(BAM_MAGIC));
            offset += sizeof(BAM_MAGIC);
            // l_text & text
            const std::string& text = 
                std::get<sam::HEADER_INDEX::PLAIN_TEXT>(rhs.header_);
            std::int32_t buf = text.size();
            if (buf > INT_MAX)
            {
                std::cerr << "ERROR: Header text length is too long.";
                std::cerr << std::endl;
                return out;
            }
            parse_char<std::int32_t>(&data[offset], buf, &offset);
            if (offset + buf > CHUNK_SIZE)
            {
                const char* c_str = text.c_str();
                std::memcpy(&data[offset], c_str, CHUNK_SIZE - offset);
                rhs.deflate_block(out, data, CHUNK_SIZE);
                std::memcpy(data, &c_str[CHUNK_SIZE - offset], 
                                    buf - CHUNK_SIZE + offset);
                offset = buf - CHUNK_SIZE + offset;
            }
            else
            {
                std::memcpy(&data[offset], text.c_str(), buf);
                offset += buf;
            }
            // n_ref & refs
            const std::vector<sam::ReferenceType>& refs = 
                std::get<sam::HEADER_INDEX::REFERENCE>(rhs.header_);
            buf = refs.size();
            parse_char<std::int32_t>(&data[offset], buf, &offset);
            for (std::size_t i = 0;i < refs.size();++i)
            {
                const std::string& ref_name = 
                    std::get<sam::REFERENCE_INDEX::REFERENCE_NAME>(refs[i]);
                buf = ref_name.size();
                if (offset + 2 * sizeof(std::int32_t) + buf > CHUNK_SIZE)
                {
                    rhs.deflate_block(out, data, offset);
                    offset = 0;
                }
                parse_char<std::int32_t>(&data[offset], buf, &offset);
                std::memcpy(&data[offset], ref_name.c_str(), buf);
                offset += buf;
                parse_char<std::int32_t>
                    (&data[offset], 
                     std::get<sam::REFERENCE_INDEX::REFERENCE_LENGTH>(refs[i]), 
                     &offset);
            }
            rhs.deflate_block(out, data, offset);
            return out;
        };

      private:
        // Need deflate_block()
        friend std::ostream& ::biovoltron::format::operator<< 
            (std::ostream& out, const BAM& rhs);
        
        /// Initialize z_stream for future inflate_block() use.
        void init_zargs()
        {
            z_args.zalloc = Z_NULL;
            z_args.zfree = Z_NULL;
            z_args.opaque = Z_NULL;
            z_args.avail_in = 0;
            z_args.next_in = Z_NULL;
            inflateInit2(&z_args, GZIP_WINDOW_BITS);
        }
    
        /**
         *  @brief Convert data stream to desired data type.
         *  @tparam T The type which you want data to becomes
         *  @param data Data pointer points to the starting address 
         *         of the data which you want to convert
         *  @return Converted data in desired data type
         *  @see parse_char(), determine_big_endianess()
         * 
         *  A private function for internal use when read data stream 
         *  from BAM file, then convert byte into variable in 
         *  specific type.
         */
        template <typename T>
        static T convert_char(const char* data)
        {
            if (is_big_endian_)
            {
                std::size_t size = sizeof(T);
                std::unique_ptr<char[]> temp(new char[size]);
                for (std::size_t i = 0;i < size;++i)
                    temp[i] = data[size - 1 - i];
                return (*reinterpret_cast<T*>(temp.get()));
            }
            else
                return (*reinterpret_cast<const T*>(data));
        }

        /**
         *  @brief Convert variable into byte (little endian).
         *  @tparam T The type which you want data output to data stream
         *  @tparam U Original type which data is (Might different from T)
         *  @param out Data stream which you want to output
         *  @param data Target variable which you want to output
         *  @param offset (Optional) Data stream index remainder, 
         *                will increase by sizeof(T). Default is nullptr
         *  @see convert_char(), determine_big_endianess()
         * 
         *  A private function for internal use when write variable 
         *  to data stream in little endian form.<BR>
         *  Template parameters T and U might be different, but typically 
         *  is the same. T is the type BAM file format provides, 
         *  so U must converted to T to match the documentation.<BR>
         *  Offset is optional parameter that could provide to increase
         *  file pointer (out) outside the function for 
         *  convenience reasons.
         */
        template <typename T, typename U>
        static void parse_char
            (char* out, U data, std::uint32_t* offset = nullptr)
        {
            std::size_t size = sizeof(T);
            if (is_big_endian_)
            {
                for (std::size_t i = 0;i < size;++i)
                    out[i] = (char)(data >> 8 * i);
            }
            else
                std::memcpy(out, &data, size);
            
            if (offset)
                *offset += size;
        }

        /**
         *  @brief Check whether system is big endian
         *  @return Whether system is big endian
         * 
         *  Because BAM file format is liitle endian, 
         *  if system is big endian then program might work improperly.
         *  So needs this function to checking endianess 
         *  in working system.
         */
        static bool determine_big_endianess()
        {
            const std::uint16_t x = 0x0001;
            return (*(reinterpret_cast<const char*>(&x))) == 0;
        }

        /**
         *  @brief Check whether BGZF header is valid.
         *  @param header Byte contains BGZF header information
         *  @return Whether BGZF header is valid
         * 
         *  Because in BGZF file format, there has a specific 
         *  gzip header format. So needs this function to determine 
         *  the file is BGZF file or not.
         */
        bool check_header_is_valid(char header[18])
        {
            if (header[0]  != GZIP_ID1
            || header[1]  != GZIP_ID2
            || header[2]  != GZIP_CM
            || header[3]  != GZIP_FLAG
            || header[12] != GZIP_SI1
            || header[13] != GZIP_SI2
            || header[14] != GZIP_SLEN) return false;
                
            return true;
        }

        /**
         *  @brief Decompress BGZF block to buffer.
         *  @param in istream contains data in compressed format
         *  @return Byte count that decompresssed
         *  @see check_header_is_valid(), deflate_block()
         * 
         *  Use zlib inflate() to decompress a block 
         *  from BGZF file format.<BR>
         *  Decomressed data will put into block_buffer_ and update 
         *  block_offset_, block_length, and block_address_.
         */
        std::uint32_t inflate_block(std::istream& in)
        {
            char header[18];
            block_address_ = in.tellg();
            in.read(header, sizeof(header));
            if (in.eof())
                return -1;
            else if (!check_header_is_valid(header))
            {
                std::cerr << "ERROR: Header format not match";
                std::cerr << std::endl;
                return -1;
            }
                
            std::int8_t ret;
            char in_char[CHUNK_SIZE];
            std::uint16_t block_size = convert_char<std::uint16_t>(&header[16]);
            z_args.avail_in = block_size - GZIP_XLEN - 19;
            in.read(in_char, z_args.avail_in);
            z_args.next_in = reinterpret_cast<Bytef*>(in_char);
            z_args.avail_out = CHUNK_SIZE;
            z_args.next_out = reinterpret_cast<Bytef*>(block_buffer_);
            ret = inflate(&z_args, Z_FINISH);
            if (ret != Z_STREAM_END)
            {
                inflateEnd(&z_args);
                std::cerr << "ERROR: inflate() failed" << std::endl;
                return -1;
            }
            char buf[4];
            in.read(buf, sizeof(std::uint32_t)); // ignore CRC
            in.read(buf, sizeof(std::uint32_t));
            std::uint32_t i_size = convert_char<std::uint32_t>(buf);
            if (i_size != z_args.total_out)
            {
                std::cerr << "ERROR: ";
                std::cerr << "ISIZE doesn't match total_out in zlib";
                std::cerr << std::endl;
                return -1;
            }
            ret = inflateReset(&z_args);
            if (ret != Z_OK)
            {
                std::cerr << "ERROR: inflateReset() failed";
                std::cerr << std::endl;
                return -1;
            }

            block_offset_ = 0;
            block_length_ = i_size;
            return z_args.total_out;
        }

        /**
         *  @brief Compress data into BGZF block.
         *  @param outfile ostream which compressed data will forward to
         *  @param data Data stream which will be compressed
         *  @param len Data length in data stream
         *  @return Byte count after compresssion
         *  @see check_header_is_valid(), inflate_block()
         * 
         *  Use zlib deflate() to compress data into 
         *  a block in BGZF file format.<BR>
         *  This is a static function is because of there has no 
         *  "buffer" in data members, just pass data stream which 
         *  want to compress then every one can use this function.
         */
        static std::uint32_t deflate_block
            (std::ostream& outfile, char* data, std::uint32_t len)
        {
            std::int8_t ret, flush;
            std::uint32_t out_length;
            z_stream dz_args;
            char out[CHUNK_SIZE + GZIP_WRAP_DATA_SIZE];
            out[0] = GZIP_ID1;
            out[1] = GZIP_ID2;
            out[2] = GZIP_CM;
            out[3] = GZIP_FLAG;
            out[10] = GZIP_XLEN;
            out[11] = '\0';
            out[12] = GZIP_SI1;
            out[13] = GZIP_SI2;
            out[14] = GZIP_SLEN;
            out[15] = '\0';

            dz_args.zalloc = Z_NULL;
            dz_args.zfree = Z_NULL;
            dz_args.opaque = Z_NULL;
            dz_args.next_in = reinterpret_cast<Bytef*>(data);
            dz_args.avail_in = len;
            dz_args.next_out = reinterpret_cast<Bytef*>(&out[18]);
            dz_args.avail_out = CHUNK_SIZE;
            ret = deflateInit2(&dz_args
                             , Z_DEFAULT_COMPRESSION
                             , Z_DEFLATED
                             , GZIP_WINDOW_BITS
                             , GZIP_MEM_LEVEL
                             , Z_DEFAULT_STRATEGY);
            if (ret != Z_OK)
            {
                std::cerr << "ERROR: deflateinit2() failed";
                std::cerr << std::endl;
                return -1;
            }
            ret = deflate(&dz_args, Z_FINISH);
            if (ret != Z_STREAM_END)
            {
                deflateEnd(&dz_args);
                std::cerr << "ERROR: deflate() failed";
                std::cerr << std::endl;
                return -1;
            }
            ret = deflateEnd(&dz_args);
            if (ret != Z_OK)
            {
                std::cerr << "ERROR: deflateEnd() failed";
                std::cerr << std::endl;
                return -1;
            }

            out_length = dz_args.total_out + GZIP_WRAP_DATA_SIZE;
            if (out_length > CHUNK_SIZE)
            {
                std::cerr << "ERROR: Compressed block isze > ";
                std::cerr << "BGZF maximum block size" << std::endl;
                return -1;
            }
            parse_char<std::uint16_t>(&out[16], (std::uint16_t)out_length - 1);
            std::uint32_t gzip_crc = crc32(0, Z_NULL, 0);
            gzip_crc = crc32(gzip_crc
                             , reinterpret_cast<Bytef*>(data), len);
            parse_char<std::uint32_t>(&out[out_length] - 8, gzip_crc);
            parse_char<std::uint32_t>(&out[out_length] - 4, len);
            outfile.write(out, out_length);
            return out_length;
        }

        /**
         *  @brief Read specific number of byte from block_buffer_.
         *  @param in istream which contains file pointer to BAM file
         *  @param data A buffer to store result and return to caller
         *  @param size Byte count that caller want to read
         *  @return Whether read success
         *  @see inflate_block()
         * 
         *  Read binary data from decompressed data in block_buffer_.<BR>
         *  If buffer has not enough data, function will call 
         *  inflate_block() to get more data from BGZF file.
         */
        bool read_byte_data
            (std::istream& in, char* data, std::uint64_t size)
        {
            std::uint32_t min, available = block_length_ - block_offset_;
            do
            {
                if (available == 0)
                {
                    if (inflate_block(in) == -1)
                        return false;
                    available = block_length_;
                }
                min = (available > size) ? size : available;
                std::memcpy(data, block_buffer_ + block_offset_, min);
                block_offset_ += min;
                available -= min;
                size -= min;
                data += min;
            } while (size != 0);
            return true;
        }

        /**
         *  @brief Seek file pointer to given virtual file offset.
         *  @param in istream which contains file pointer to BGZF file
         *  @param offset virtual file offset in BGZF documentation
         *  @see tell()
         * 
         *  Seek istream's file pointer to virtual file offset.<BR>
         *  Might be called when random access is activated, 
         *  then istream's file pointer will move to starting address 
         *  of a block and inflate_block() will be invoked to 
         *  decompress block. Finallly block_offset_ will seek to 
         *  virtual file offset.
         */
        void seek(std::istream& in, std::uint64_t offset)
        {
            in.seekg(offset >> BAM_OFFSET_SHIFT, std::ios::beg);
            if (!in)
            {
                std::cerr << "ERROR: ";
                std::cerr << "coffset seek to wrong file position";
                std::cerr << std::endl;
                return;
            }
            else if (inflate_block(in) == -1)
                return;
            else if ((offset & 0xffff) > block_length_)
            {
                std::cerr << "ERROR: ";
                std::cerr << "uoffset seek to wrong file position";
                std::cerr << std::endl;
                return;
            }
            block_offset_ = offset & 0xffff;
        }

        /**
         *  @brief Tell caller the current virtual file offset.
         *  @return Current virtual file offset
         *  @see seek()
         */
        std::uint64_t tell()
        {
            return (block_address_ << BAM_OFFSET_SHIFT | block_offset_);
        }

        /// BGZF file format's ID1 field
        static const char GZIP_ID1 = 31;
        /// BGZF file format's ID2 field
        static const char GZIP_ID2 = static_cast<char>(139);
        /// BGZF file format's CM field
        static const char GZIP_CM = 8;
        /// BGZF file format's FLAG field
        static const char GZIP_FLAG = 4;
        /// BGZF file format's XLEN field
        static const char GZIP_XLEN = 6;
        /// BGZF file format's SI1 field
        static const char GZIP_SI1 = 66;
        /// BGZF file format's SI2 field
        static const char GZIP_SI2 = 67;
        /// BGZF file format's SLEN field
        static const char GZIP_SLEN = 2;
        /// Represents raw inflate/deflate for zlib
        static const std::int8_t GZIP_WINDOW_BITS = -15;
        /// Default memory level in zlib
        static const std::int8_t GZIP_MEM_LEVEL = 8;
        /// BGZF file format's block size without compresssed data size
        static const std::uint8_t GZIP_WRAP_DATA_SIZE = 26;
        /// BGZF file format's maximal block size
        static const std::uint32_t CHUNK_SIZE = 65536;
        /// The bits need to shift for block starting address in 
        /// BGZF's virtual file offset
        static const std::uint8_t BAM_OFFSET_SHIFT = 16;
        /// BAM file format's magic string
        static constexpr char BAM_MAGIC[4] = {'B', 'A', 'M', '\1'};

        /// Records offset in this block
        std::uint32_t block_offset_;
        /// Records this block's length
        std::uint32_t block_length_;
        /// Records this block's starting address
        std::uint64_t block_address_;
        /// Buffers decompressed data in this block
        char block_buffer_[CHUNK_SIZE];
        /// zlib arguments
        z_stream z_args;
        /// Records system's endianess
        static bool is_big_endian_;
    };
    bool Header::is_big_endian_ = Header::determine_big_endianess();

    // Template specification for float type
    // because of float endianess convertion is different from other type.
    template <>
    void Header::parse_char<float, float>
        (char* out, float data, std::uint32_t* offset)
    {
        if (is_big_endian_)
        {
            char* data_ptr = reinterpret_cast<char*>(&data);
            out[0] = data_ptr[3];
            out[1] = data_ptr[2];
            out[2] = data_ptr[1];
            out[3] = data_ptr[0];
        }
        else
            std::memcpy(out, &data, sizeof(float));

        if (offset)
            *offset += sizeof(float);
    }
};
