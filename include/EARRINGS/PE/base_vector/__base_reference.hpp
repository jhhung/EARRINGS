#ifndef BIOVOLTRON___BASE_REFERENCE
#define BIOVOLTRON___BASE_REFERENCE

#include <cstdint>
#include <type_traits>
#include <iterator>

// #define _RNA_SEQUENCE

/**
* @file __base_reference.hpp
* @brief defines the base implementation of base_vector.hpp file
*
* @author JHH corp
*/


namespace EARRINGS
{
/// the vector accept char as the unit of input
typedef char          char_type;
/// the vector accept unsigned char as the unit of input
typedef unsigned char uint_type;
/// this defines the base conversion between char and int
constexpr uint_type __i[] = {0, 0, 1, 0, 3, 0, 2};
#ifndef _RNA_SEQUENCE
/// the characters of RNA sequence if _RNA_SEQUENCE is defined
constexpr char_type __c[] = {'A', 'C', 'G', 'T'};
#else
/// the characters of DNA sequence
constexpr char_type __c[] = {'A', 'C', 'G', 'U'};
#endif

/**
* @brief convert char_type character to uint_type number
*
* @param a char_type character
* @return the unit_type number corresponding to char_type character
*
* In order to store charaters as numbers in a storage-efficiency way,
* we do some operations so that the four characters can be 
* stored using only 2 bits. To achieve this goal, an "&" operation of
* 6 will be done to the ascii number of the char_type character. 
* This operation gives us the index of the uint_type __i[] array.
* The value (let's say, i) after looking up the array is then stored
* in the memory. To do the conversion back to the origin character, 
* we check the ith position in char_type __c[] array. Then we can 
* simply get back our orgin character.
* 
* For example, the ascii number of 'A' is 65. After doing "&" operation
* with 6, we get 0. Giving a glimpse of the array uint_type __i, 
* we know that __i[0] = 0. The number, 0, in this case, is then stored
* into the memory. To get back the orgin character, we look up char_type
* __c[0], which gives us 'A'. And It is exactly the same as our input.
* Similarly, the ascii of the lowercase 'g' is 103. The operation 
* gives us 6. __i[6] = '2' and __c[2] = 'G'. We again get back to the
* orgin character.
* Notice that both uppercase and lowercase is converted back to the 
* same uppercase character. So the users are free to input either 
* uppercase or lowercase character.
*/
static constexpr uint_type to_uint_type(char_type c) noexcept {return __i[c & 6u];}

/**
* @brief convert uint_type number to char_type character
*
* @param an uint_type number
* @return the char_type character correspoding to the unit_type number
*
* @sa to_uint_type
*/
static constexpr char_type to_char_type(uint_type c) noexcept {return __c[c];}


template <class Cp, bool IsConst, class CharT = typename Cp::value_type, typename Cp::__storage_type = 0> class __base_iterator;
template <class Cp, class CharT = typename Cp::value_type> class __base_const_reference;

template <class Cp, class CharT = typename Cp::value_type>
class __base_reference
{
    typedef typename Cp::__storage_type    __storage_type;
    typedef typename Cp::__storage_pointer __storage_pointer;

    __storage_pointer    seg_;
    unsigned             shift_;
    const __storage_type mask_ = 0b11;

    friend typename Cp::__self;

    friend class __base_const_reference<Cp, CharT>;
    friend class __base_iterator<Cp, false, CharT>;
public:
    inline operator CharT() const noexcept
    {
        if constexpr (std::is_same_v<CharT, uint_type>)
            return              static_cast<CharT>(*seg_ >> shift_ & mask_);
        else
            return to_char_type(static_cast<CharT>(*seg_ >> shift_ & mask_));
    }

    inline CharT operator ~() const noexcept
    {
        if constexpr (std::is_same_v<CharT, uint_type>)
            return              static_cast<CharT>((*seg_ >> shift_ ^ mask_) & mask_);
        else
            return to_char_type(static_cast<CharT>((*seg_ >> shift_ ^ mask_) & mask_));
    }

    inline
    __base_reference& operator=(CharT x) noexcept
    {
        *seg_ &= ~(mask_ << shift_);
        if constexpr (std::is_same_v<CharT, uint_type>)
            *seg_ |= __storage_type(x) << shift_;
        else
            *seg_ |= __storage_type(to_uint_type(x)) << shift_;
        return *this;
    }

    inline
    __base_reference& operator=(const __base_reference& x) noexcept
    {return operator=(static_cast<CharT>(x));}

    inline void flip() noexcept {*seg_ ^= (mask_ << shift_);}
    inline __base_iterator<Cp, false, CharT> operator&() const noexcept
    {return __base_iterator<Cp, false, CharT>(seg_, shift_ / 2);}
    
    inline __storage_pointer get_seg() const noexcept {return seg_;}
private:
    inline
    __base_reference(__storage_pointer seg, unsigned pos) noexcept
            : seg_(seg), shift_(pos * 2) {}
};

template <class Cp, class CharT>
class __base_const_reference
{
    typedef typename Cp::__storage_type          __storage_type;
    typedef typename Cp::__const_storage_pointer __storage_pointer;

    __storage_pointer    seg_;
    unsigned             shift_;
    const __storage_type mask_ = 0b11;

    friend typename Cp::__self;
    friend class __base_iterator<Cp, true, CharT>;
public:
    inline
    __base_const_reference(const __base_reference<Cp, CharT>& x) noexcept
            : seg_(x.seg_), shift_(x.shift_) {}

    inline operator CharT() const noexcept
    {
        if constexpr (std::is_same_v<CharT, uint_type>)
            return              static_cast<CharT>(*seg_ >> shift_ & mask_);
        else
            return to_char_type(static_cast<CharT>(*seg_ >> shift_ & mask_));
    }

    inline __base_iterator<Cp, true, CharT> operator&() const noexcept
    {return __base_iterator<Cp, true, CharT>(seg_, shift_ / 2);}

    inline __storage_pointer get_seg() const noexcept {return seg_;}
private:
    inline
    constexpr
    __base_const_reference(__storage_pointer seg, unsigned pos) noexcept
            : seg_(seg), shift_(pos * 2) {}

    __base_const_reference& operator=(const __base_const_reference& x);
};

/**
* @class __base_iterator
* @brief 
*
* @tparam Cp: a vector
* @tparam IsConst: to check whether the vector is const or not
* @tparam CharT: the input type of the vector 
* @tparam Cp::__storage_type
* 
* Member variable:<br>
* size_t chunk_size_<br>
* SeqVec seq_<br>
* std::string name_<br>
*
*/

template <class Cp, bool IsConst, class CharT, typename Cp::__storage_type>
class __base_iterator
{
public:
    typedef typename Cp::difference_type                                               difference_type;
    typedef CharT                                                                      value_type;
    typedef __base_iterator                                                            pointer;
    typedef std::conditional_t
            <IsConst, __base_const_reference<Cp, CharT>, __base_reference<Cp, CharT> > reference;
    typedef std::random_access_iterator_tag                                            iterator_category;

private:
    typedef typename Cp::__storage_type                                 __storage_type;
    typedef std::conditional_t<IsConst, typename Cp::__const_storage_pointer,
                                        typename Cp::__storage_pointer> __storage_pointer;
    static const unsigned bases_per_word = Cp::bases_per_word;

    __storage_pointer seg_;
    unsigned          pos_;

public:
    inline __base_iterator() noexcept
            : seg_(nullptr), pos_(0)
    {}

    inline
    __base_iterator(const __base_iterator<Cp, false, CharT>& it) noexcept
            : seg_(it.seg_), pos_(it.pos_) {}

    inline reference operator*() const noexcept
    { return reference(seg_, pos_); }

    inline __base_iterator& operator++()
    {
        if (pos_ != bases_per_word - 1)
            ++pos_;
        else
        {
            pos_ = 0;
            ++seg_;
        }
        return *this;
    }

    inline __base_iterator operator++(int)
    {
        __base_iterator tmp = *this;
        ++(*this);
        return tmp;
    }

    inline __base_iterator& operator--()
    {
        if (pos_ != 0)
            --pos_;
        else
        {
            pos_ = bases_per_word - 1;
            --seg_;
        }
        return *this;
    }

    inline __base_iterator operator--(int)
    {
        __base_iterator tmp = *this;
        --(*this);
        return tmp;
    }

    inline __base_iterator& operator+=(difference_type n)
    {
        if (n >= 0)
            seg_ += (n + pos_) / bases_per_word;
        else
            seg_ += static_cast<difference_type>(n - bases_per_word + pos_ + 1)
                    / static_cast<difference_type>(bases_per_word);
        n &= (bases_per_word - 1);
        pos_ = static_cast<unsigned>((n + pos_) % bases_per_word);
        return *this;
    }

    inline __base_iterator& operator-=(difference_type n)
    {
        return *this += -n;
    }

    inline __base_iterator operator+(difference_type n) const
    {
        __base_iterator t(*this);
        t += n;
        return t;
    }

    inline __base_iterator operator-(difference_type n) const
    {
        __base_iterator t(*this);
        t -= n;
        return t;
    }

    inline
    friend __base_iterator operator+(difference_type n, const __base_iterator& it) { return it + n; }

    
    inline
    friend difference_type operator-(const __base_iterator& x, const __base_iterator& y)
    {
        return (x.seg_ - y.seg_) * bases_per_word + x.pos_ - y.pos_;
    }

    /// Get the nth element of __base_iterator
    inline reference operator[](difference_type n) const { return *(*this + n); }

    /// Check if x == y
    inline friend bool operator==(const __base_iterator& x, const __base_iterator& y)
    {
        return x.seg_ == y.seg_ && x.pos_ == y.pos_;
    }

    /// Check if x != y
    inline friend bool operator!=(const __base_iterator& x, const __base_iterator& y)
    {
        return !(x == y);
    }

    /// Check if x < y
    inline friend bool operator<(const __base_iterator& x, const __base_iterator& y)
    {
        return x.seg_ < y.seg_ || (x.seg_ == y.seg_ && x.pos_ < y.pos_);
    }

    /// Check if x > y
    inline friend bool operator>(const __base_iterator& x, const __base_iterator& y)
    {
        return y < x;
    }

    /// Check if x <= y
    inline friend bool operator<=(const __base_iterator& x, const __base_iterator& y)
    {
        return !(y < x);
    }

    /// Check if x >= y
    inline friend bool operator>=(const __base_iterator& x, const __base_iterator& y)
    {
        return !(x < y);
    }

private:
    inline
    __base_iterator(__storage_pointer seg, unsigned pos) noexcept
            : seg_(seg), pos_(pos) {}

    friend typename Cp::__self;

    friend class __base_reference<Cp, CharT>;
    friend class __base_const_reference<Cp, CharT>;
    friend class __base_iterator<Cp, true, CharT>;
};

/// swap the content of  two __base_reference
template <class Cp, class CharT>
inline
void
swap(__base_reference<Cp, CharT> x, __base_reference<Cp, CharT> y) noexcept
{
    CharT t = x;
    x = y;
    y = t;
}

template <class Cp, class Dp, class CharT>
inline
void
swap(__base_reference<Cp, CharT> x, __base_reference<Dp, CharT> y) noexcept
{
    CharT t = x;
    x = y;
    y = t;
}

template <class Cp, class CharT>
inline
void
swap(__base_reference<Cp, CharT> x, CharT& y) noexcept
{
    CharT t = x;
    x = y;
    y = t;
}

}

#endif //BIOVOLTRON___BASE_REFERENCE
