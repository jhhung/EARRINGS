#pragma once

/**
 * @file base_vector.hpp
 * @brief a special vector that optimize the storage of gene sequences
 *
 * @author JHH corp
 */

/*

namespace biovoltron
{

template <class CharT, class Allocator = std::allocator<CharT> >
class vector
{
public:
    typedef CharT                                    value_type;
    typedef Allocator                                allocator_type;
    typedef implementation-defined                   iterator;
    typedef implementation-defined                   const_iterator;
    typedef typename allocator_type::size_type       size_type;
    typedef typename allocator_type::difference_type difference_type;
    typedef iterator                                 pointer;
    typedef const_iterator                           const_pointer;
    typedef std::reverse_iterator<iterator>          reverse_iterator;
    typedef std::reverse_iterator<const_iterator>    const_reverse_iterator;

    class reference
    {
    public:
        reference(const reference&) noexcept;
        operator value_type() const noexcept;
        reference& operator=(const value_type x) noexcept;
        reference& operator=(const reference& x) noexcept;
        iterator operator&() const noexcept;
        void flip() noexcept;
    };

    class const_reference
    {
    public:
        const_reference(const reference&) noexcept;
        operator value_type() const noexcept;
        const_iterator operator&() const noexcept;
    };

    vector()
        noexcept(is_nothrow_default_constructible<allocator_type>::value);
    explicit vector(const allocator_type&);
    explicit vector(size_type n, const allocator_type& a = allocator_type()); // C++14
    vector(size_type n, const value_type& value, const allocator_type& = allocator_type());
    template <class InputIterator>
        vector(InputIterator first, InputIterator last, const allocator_type& = allocator_type());
    vector(const vector& x);
    vector(vector&& x)
        noexcept(is_nothrow_move_constructible<allocator_type>::value);
    vector(initializer_list<value_type> il);
    vector(initializer_list<value_type> il, const allocator_type& a);
    ~vector();
    vector& operator=(const vector& x);
    vector& operator=(vector&& x)
        noexcept(
             allocator_type::propagate_on_container_move_assignment::value ||
             allocator_type::is_always_equal::value); // C++17
    vector& operator=(initializer_list<value_type> il);
    template <class InputIterator>
        void assign(InputIterator first, InputIterator last);
    void assign(size_type n, const value_type& u);
    void assign(initializer_list<value_type> il);

    allocator_type get_allocator() const noexcept;

    iterator               begin() noexcept;
    const_iterator         begin()   const noexcept;
    iterator               end() noexcept;
    const_iterator         end()     const noexcept;

    reverse_iterator       rbegin() noexcept;
    const_reverse_iterator rbegin()  const noexcept;
    reverse_iterator       rend() noexcept;
    const_reverse_iterator rend()    const noexcept;

    const_iterator         cbegin()  const noexcept;
    const_iterator         cend()    const noexcept;
    const_reverse_iterator crbegin() const noexcept;
    const_reverse_iterator crend()   const noexcept;

    size_type size() const noexcept;
    size_type max_size() const noexcept;
    size_type capacity() const noexcept;
    bool empty() const noexcept;
    void reserve(size_type n);
    void shrink_to_fit() noexcept;

    reference       operator[](size_type n);
    const_reference operator[](size_type n) const;
    reference       at(size_type n);
    const_reference at(size_type n) const;

    reference       front();
    const_reference front() const;
    reference       back();
    const_reference back() const;

    storage_pointer       data() noexcept;
    const storage_pointer data() const noexcept;

    void push_back(const value_type& x);
    template <class... Args> reference emplace_back(Args&&... args);  // C++14; reference in C++17
    void pop_back();

    template <class... Args> iterator emplace(const_iterator position, Args&&... args);  // C++14
    iterator insert(const_iterator position, const value_type& x);
    iterator insert(const_iterator position, size_type n, const value_type& x);
    template <class InputIterator>
        iterator insert(const_iterator position, InputIterator first, InputIterator last);
    iterator insert(const_iterator position, initializer_list<value_type> il);

    iterator erase(const_iterator position);
    iterator erase(const_iterator first, const_iterator last);

    void clear() noexcept;

    void resize(size_type sz);
    void resize(size_type sz, value_type x);
    void swap(vector&)
        noexcept(allocator_traits<allocator_type>::propagate_on_container_swap::value ||
                 allocator_traits<allocator_type>::is_always_equal::value);  // C++17
    void flip() noexcept;

    bool __invariants() const;

    // for value_type = char_type:

    vector(const value_type* s, const allocator_type& a = allocator_type());
    vector(const std::string& s, const allocator_type& a = allocator_type());

    vector& operator=(const value_type* s);
    vector& operator=(const std::string& s);
    vector& operator=(value_type c);

    vector& operator+=(const vector& v);
    vector& operator+=(const std::string& s);
    vector& operator+=(const value_type* s);
    vector& operator+=(value_type c);
    vector& operator+=(std::initializer_list<value_type>);

};

template <class CharT, class Allocator> bool operator==(const vector<CharT,Allocator>& x, const vector<CharT,Allocator>& y);
template <class CharT, class Allocator> bool operator< (const vector<CharT,Allocator>& x, const vector<CharT,Allocator>& y);
template <class CharT, class Allocator> bool operator!=(const vector<CharT,Allocator>& x, const vector<CharT,Allocator>& y);
template <class CharT, class Allocator> bool operator> (const vector<CharT,Allocator>& x, const vector<CharT,Allocator>& y);
template <class CharT, class Allocator> bool operator>=(const vector<CharT,Allocator>& x, const vector<CharT,Allocator>& y);
template <class CharT, class Allocator> bool operator<=(const vector<CharT,Allocator>& x, const vector<CharT,Allocator>& y);

template <class CharT, class Allocator>
void swap(vector<CharT,Allocator>& x, vector<CharT,Allocator>& y)
    noexcept(noexcept(x.swap(y)));

// for value_type = char_type:

template<class CharT, class Allocator> bool operator==(const CharT* lhs, const vector<CharT, Allocator>& rhs) noexcept;
template<class CharT, class Allocator> bool operator==(const vector<CharT, Allocator>& lhs, const CharT* rhs) noexcept;
template<class CharT, class Allocator> bool operator==(const std::string& lhs, const vector<CharT, Allocator>& rhs) noexcept;
template<class CharT, class Allocator> bool operator==(const vector<CharT, Allocator>& lhs, const std::string& rhs) noexcept;

template<class CharT, class Allocator> bool operator!=(const CharT* lhs, const vector<CharT, Allocator>& rhs) noexcept;
template<class CharT, class Allocator> bool operator!=(const vector<CharT, Allocator>& lhs, const CharT* rhs) noexcept;
template<class CharT, class Allocator> bool operator!=(const std::string& lhs, const vector<CharT, Allocator>& rhs) noexcept;
template<class CharT, class Allocator> bool operator!=(const vector<CharT, Allocator>& lhs, const std::string& rhs) noexcept;

template<class CharT, class Allocator> bool operator< (const vector<CharT, Allocator>& lhs, const CharT* rhs) noexcept;
template<class CharT, class Allocator> bool operator< (const CharT* lhs, const vector<CharT, Allocator>& rhs) noexcept;
template<class CharT, class Allocator> bool operator< (const vector<CharT, Allocator>& lhs, const std::string& rhs) noexcept;
template<class CharT, class Allocator> bool operator< (const std::string& lhs, const vector<CharT, Allocator>& rhs) noexcept;

template<class CharT, class Allocator> bool operator> (const vector<CharT, Allocator>& lhs, const CharT* rhs) noexcept;
template<class CharT, class Allocator> bool operator> (const CharT* lhs, const vector<CharT, Allocator>& rhs) noexcept;
template<class CharT, class Allocator> bool operator> (const vector<CharT, Allocator>& lhs, const std::string& rhs) noexcept;
template<class CharT, class Allocator> bool operator> (const std::string& lhs, const vector<CharT, Allocator>& rhs) noexcept;

template<class CharT, class Allocator> bool operator<=(const vector<CharT, Allocator>& lhs, const CharT* rhs) noexcept;
template<class CharT, class Allocator> bool operator<=(const CharT* lhs, const vector<CharT, Allocator>& rhs) noexcept;
template<class CharT, class Allocator> bool operator<=(const vector<CharT, Allocator>& lhs, const std::string& rhs) noexcept;
template<class CharT, class Allocator> bool operator<=(const std::string& lhs, const vector<CharT, Allocator>& rhs) noexcept;

template<class CharT, class Allocator> bool operator>=(const vector<CharT, Allocator>& lhs, const CharT* rhs) noexcept;
template<class CharT, class Allocator> bool operator>=(const CharT* lhs, const vector<CharT, Allocator>& rhs) noexcept;
template<class CharT, class Allocator> bool operator>=(const vector<CharT, Allocator>& lhs, const std::string& rhs) noexcept;
template<class CharT, class Allocator> bool operator>=(const std::string& lhs, const vector<CharT, Allocator>& rhs) noexcept;

template<class CharT, class Allocator>
std::basic_ostream<CharT>&
operator<<(std::basic_ostream<CharT>& os, const vector<CharT, Allocator>& v);

}  // biovoltron

*/

#include <cassert>
#include "__base_reference.hpp"
#include <limits>
#include <climits>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <iostream>
#include <algorithm>
#include <string_view>
#include <vector>
#include <utility>
#include <cmath>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

namespace biovoltron
{

template <typename InputIterator>
inline constexpr bool __is_input_iterator_v =
        std::is_convertible_v<typename std::iterator_traits<InputIterator>::iterator_category,
                std::input_iterator_tag>;

template <typename InputIterator>
inline constexpr bool __is_forward_iterator_v =
        std::is_convertible_v<typename std::iterator_traits<InputIterator>::iterator_category,
                std::forward_iterator_tag>;

template <typename Alloc, typename Traits = std::allocator_traits<Alloc>>
inline constexpr bool __is_noexcept_move_assign_container_v =
        Traits::propagate_on_container_move_assignment::value || Traits::is_always_equal::value;

template <bool>
class __vector_base_common
{
protected:
    __vector_base_common() = default;
    void __throw_length_error() const;
    void __throw_out_of_range() const;
};

template <bool b>
void
__vector_base_common<b>::__throw_length_error() const
{
    throw std::length_error("vector");
}

template <bool b>
void
__vector_base_common<b>::__throw_out_of_range() const
{
    throw std::out_of_range("vector");
}

template <class CharT, class Allocator = std::allocator<CharT> >
class vector : private __vector_base_common<true>
{
    static_assert(std::is_same_v<CharT, char_type> || std::is_same_v<CharT, uint_type>, "CharT must be char_type or uint_type");

public:
    /// define vector as __self
    typedef vector                                   __self;
    /// CharT is either char_type or uint_type
    typedef CharT                                    value_type;
    typedef std::char_traits<CharT>                  traits_type;
    typedef Allocator                                allocator_type;
    typedef std::allocator_traits<allocator_type>    __alloc_traits;
    typedef typename __alloc_traits::size_type       size_type;
    typedef typename __alloc_traits::difference_type difference_type;
    typedef size_type                                __storage_type;
    typedef __base_iterator<vector, false>           pointer;
    typedef __base_iterator<vector, true>            const_pointer;
    typedef pointer                                  iterator;
    typedef const_pointer                            const_iterator;
    typedef std::reverse_iterator<iterator>          reverse_iterator;
    typedef std::reverse_iterator<const_iterator>    const_reverse_iterator;

    static const unsigned bases_per_word = static_cast<unsigned>(sizeof(__storage_type) * CHAR_BIT / 2);
private:
    typedef typename __alloc_traits::template rebind_alloc<__storage_type> __storage_allocator;
    typedef std::allocator_traits<__storage_allocator> __storage_traits;
    typedef typename __storage_traits::pointer         __storage_pointer;
    typedef typename __storage_traits::const_pointer   __const_storage_pointer;

    __storage_pointer                                      begin_;
    size_type                                              size_;
    std::pair<size_type, __storage_allocator> cap_alloc_;
public:
    typedef __base_reference      <vector> reference;
    typedef __base_const_reference<vector> const_reference;
private:
    size_type& __cap() noexcept
    {return cap_alloc_.first;}

    const size_type& __cap() const noexcept
    {return cap_alloc_.first;}

    __storage_allocator& __alloc() noexcept
    {return cap_alloc_.second;}

    const __storage_allocator& __alloc() const noexcept
    {return cap_alloc_.second;}

    static size_type __internal_cap_to_external(size_type n) noexcept
    {return n * bases_per_word;}

    static size_type __external_cap_to_internal(size_type n) noexcept
    {return (n - 1) / bases_per_word + 1;}

public:
    vector() noexcept(std::is_nothrow_default_constructible_v<allocator_type>);

    explicit vector(const allocator_type& a) noexcept;
    ~vector();
    explicit vector(size_type n);
    explicit vector(size_type n, const allocator_type& a);
    vector(size_type n, const value_type& x);
    vector(size_type n, const value_type& x, const allocator_type& a);
    template <class InputIterator>
    vector(InputIterator first, InputIterator last,
           std::enable_if_t<__is_input_iterator_v  <InputIterator> &&
                           !__is_forward_iterator_v<InputIterator>>* = 0);
    template <class InputIterator>
    vector(InputIterator first, InputIterator last, const allocator_type& a,
           std::enable_if_t<__is_input_iterator_v  <InputIterator> &&
                           !__is_forward_iterator_v<InputIterator>>* = 0);
    template <class ForwardIterator>
    vector(ForwardIterator first, ForwardIterator last,
           std::enable_if_t<__is_forward_iterator_v<ForwardIterator>>* = 0);
    template <class ForwardIterator>
    vector(ForwardIterator first, ForwardIterator last, const allocator_type& a,
           std::enable_if_t<__is_forward_iterator_v<ForwardIterator>>* = 0);

    vector(const vector& v);
    vector(const vector& v, const allocator_type& a);
    vector& operator=(const vector& v);

    vector(std::initializer_list<value_type> il);
    vector(std::initializer_list<value_type> il, const allocator_type& a);

    vector(vector&& v) noexcept;
    vector(vector&& v, const allocator_type& a);

    vector& operator=(vector&& v)
    noexcept(__is_noexcept_move_assign_container_v<Allocator, __alloc_traits>);

    vector& operator=(std::initializer_list<value_type> il)
    {assign(il.begin(), il.end()); return *this;}

    template <class InputIterator>
    std::enable_if_t
    <
        __is_input_iterator_v  <InputIterator> &&
       !__is_forward_iterator_v<InputIterator>
    >
    assign(InputIterator first, InputIterator last);
    template <class ForwardIterator>
    std::enable_if_t
    <
        __is_forward_iterator_v<ForwardIterator>
    >
    assign(ForwardIterator first, ForwardIterator last);

    void assign(size_type n, const value_type& x);

    void assign(std::initializer_list<value_type> il)
    {assign(il.begin(), il.end());}

    allocator_type get_allocator() const noexcept
    {return allocator_type(this->__alloc());}

    size_type max_size() const noexcept;

    size_type capacity() const noexcept
    {return __internal_cap_to_external(__cap());}

    size_type size() const noexcept
    {return size_;}

    bool empty() const noexcept
    {return size_ == 0;}
    void reserve(size_type n);
    void shrink_to_fit() noexcept;

    iterator begin() noexcept
    {return __make_iter(0);}

    const_iterator begin() const noexcept
    {return __make_iter(0);}

    iterator end() noexcept
    {return __make_iter(size_);}

    const_iterator end() const noexcept
    {return __make_iter(size_);}

    reverse_iterator rbegin() noexcept
    {return       reverse_iterator(end());}

    const_reverse_iterator rbegin() const noexcept
    {return const_reverse_iterator(end());}

    reverse_iterator rend() noexcept
    {return       reverse_iterator(begin());}

    const_reverse_iterator rend()   const noexcept
    {return const_reverse_iterator(begin());}

    const_iterator         cbegin()  const noexcept
    {return __make_iter(0);}

    const_iterator         cend()    const noexcept
    {return __make_iter(size_);}

    const_reverse_iterator crbegin() const noexcept
    {return rbegin();}

    const_reverse_iterator crend()   const noexcept
    {return rend();}

    reference       operator[](size_type n)       {return __make_ref(n);}
    const_reference operator[](size_type n) const {return __make_ref(n);}
    reference       at(size_type n);
    const_reference at(size_type n) const;

    reference       front()       {return __make_ref(0);}
    const_reference front() const {return __make_ref(0);}
    reference       back()        {return __make_ref(size_ - 1);}
    const_reference back()  const {return __make_ref(size_ - 1);}

    __storage_pointer       data() noexcept
    {return begin_;}
    const __storage_pointer data() const noexcept
    {return begin_;}

    void push_back(const value_type& x);
    template <class... Args>
    reference emplace_back(Args&&... args)
    {
        push_back(value_type(std::forward<Args>(args)...));
        return this->back();
    }

    void pop_back() {--size_;}

    template <class... Args>
    iterator emplace(const_iterator position, Args&&... args)
    {return insert(position, value_type(std::forward<Args>(args)...));}

    iterator insert(const_iterator position, const value_type& x);
    iterator insert(const_iterator position, size_type n, const value_type& x);
    template <class InputIterator>
    std::enable_if_t
    <
        __is_input_iterator_v  <InputIterator> &&
       !__is_forward_iterator_v<InputIterator>,
        iterator
    >
    insert(const_iterator position, InputIterator first, InputIterator last);
    template <class ForwardIterator>
    std::enable_if_t
    <
        __is_forward_iterator_v<ForwardIterator>,
        iterator
    >
    insert(const_iterator position, ForwardIterator first, ForwardIterator last);

    iterator insert(const_iterator position, std::initializer_list<value_type> il)
    {return insert(position, il.begin(), il.end());}

    iterator erase(const_iterator position);
    iterator erase(const_iterator first, const_iterator last);

    void clear() noexcept {size_ = 0;}

    void swap(vector&) noexcept;

    static void swap(reference x, reference y) noexcept {biovoltron::swap(x, y);}

    void resize(size_type sz, value_type x = 0);
    void flip() noexcept;

    bool __invariants() const;

    // for value_type = char_type:
    template <typename charT = CharT>
    vector(const value_type* s,
           std::enable_if_t<std::is_same_v<charT, char_type>>* = 0);
    template <typename charT = CharT>
    vector(const value_type* s, const allocator_type& a,
           std::enable_if_t<std::is_same_v<charT, char_type>>* = 0);
    template <typename charT = CharT>
    vector(const std::string& s,
           std::enable_if_t<std::is_same_v<charT, char_type>>* = 0);
    template <typename charT = CharT>
    vector(const std::string& s, const allocator_type& a,
           std::enable_if_t<std::is_same_v<charT, char_type>>* = 0);

    template <typename charT = CharT>
    std::enable_if_t<std::is_same_v<charT, char_type>, vector&>
    operator=(const value_type* s);
    template <typename charT = CharT>
    std::enable_if_t<std::is_same_v<charT, char_type>, vector&>
    operator=(const std::string& s);
    template <typename charT = CharT>
    std::enable_if_t<std::is_same_v<charT, char_type>, vector&>
    operator=(value_type c);

    template <typename charT = CharT>
    std::enable_if_t<std::is_same_v<charT, char_type>, vector&>
    operator+=(const vector& v);
    template <typename charT = CharT>
    std::enable_if_t<std::is_same_v<charT, char_type>, vector&>
    operator+=(const std::string& s);
    template <typename charT = CharT>
    std::enable_if_t<std::is_same_v<charT, char_type>, vector&>
    operator+=(const value_type* s);
    template <typename charT = CharT>
    std::enable_if_t<std::is_same_v<charT, char_type>, vector&>
    operator+=(value_type c);
    template <typename charT = CharT>
    std::enable_if_t<std::is_same_v<charT, char_type>, vector&>
    operator+=(std::initializer_list<value_type>);

private:
    void __invalidate_all_iterators();
    void __vallocate(size_type n);
    void __vdeallocate() noexcept;

    static size_type __align_it(size_type new_size) noexcept
    {return new_size + (bases_per_word-1) & ~((size_type)bases_per_word-1);}
    size_type __recommend(size_type new_size) const;
    void __construct_at_end(size_type n, value_type x);
    template <class ForwardIterator>
    std::enable_if_t
    <
        __is_forward_iterator_v<ForwardIterator>
    >
    __construct_at_end(ForwardIterator first, ForwardIterator last);

    reference __make_ref(size_type pos) noexcept
    {return reference(begin_ + pos / bases_per_word, static_cast<unsigned>(pos % bases_per_word));}

    const_reference __make_ref(size_type pos) const noexcept
    {return const_reference(begin_ + pos / bases_per_word, static_cast<unsigned>(pos % bases_per_word));}

    iterator __make_iter(size_type pos) noexcept
    {return iterator(begin_ + pos / bases_per_word, static_cast<unsigned>(pos % bases_per_word));}

    const_iterator __make_iter(size_type pos) const noexcept
    {return const_iterator(begin_ + pos / bases_per_word, static_cast<unsigned>(pos % bases_per_word));}

    iterator __const_iterator_cast(const_iterator p) noexcept
    {return begin() + (p - cbegin());}

    void __copy_assign_alloc(const vector& v)
    {__copy_assign_alloc(v, std::bool_constant<
                __storage_traits::propagate_on_container_copy_assignment::value>());}

    void __copy_assign_alloc(const vector& c, std::true_type)
    {
        if (__alloc() != c.__alloc())
            __vdeallocate();
        __alloc() = c.__alloc();
    }

    void __copy_assign_alloc(const vector&, std::false_type) {}

    void __move_assign(vector& c, std::false_type);
    void __move_assign(vector& c, std::true_type)
    noexcept(std::is_nothrow_move_assignable_v<allocator_type>);

    void __move_assign_alloc(vector& c)
    noexcept(
    !__storage_traits::propagate_on_container_move_assignment::value ||
    std::is_nothrow_move_assignable_v<allocator_type>)
    {__move_assign_alloc(c, std::bool_constant<
                __storage_traits::propagate_on_container_move_assignment::value>());}

    void __move_assign_alloc(vector& c, std::true_type)
    noexcept(std::is_nothrow_move_assignable_v<allocator_type>)
    {
        __alloc() = std::move(c.__alloc());
    }

    void __move_assign_alloc(vector&, std::false_type) noexcept {}

	template<class Archive>
	void save(Archive& ar, const unsigned int version) const
	{
		std::vector<size_type> buf;
		size_type seg_size((size_ - 1) / bases_per_word + 1);

		buf.reserve(seg_size);
		for (size_type i(0); i < seg_size; i++)
			buf.emplace_back(*(begin_ + i));

		ar << size_ << buf << cap_alloc_.first;
	}
	
	template<class Archive>
	void load(Archive& ar, const unsigned int version)
	{
		size_type size, seg_size;
		std::vector<size_type> buf;
		
		ar >> size;
		ar >> buf >> cap_alloc_.first;

		reserve(__recommend(size));
		size_ = size;
		seg_size = (size - 1) / bases_per_word + 1;
		for (size_type i(0); i < seg_size; i++)
			*(begin_ + i) = buf[i];
	}

	BOOST_SERIALIZATION_SPLIT_MEMBER()

    friend class __base_reference<vector>;
    friend class __base_const_reference<vector>;
    friend class __base_iterator<vector, false>;
    friend class __base_iterator<vector, true>;
	friend class boost::serialization::access;
};

template <class CharT, class Allocator>
void
vector<CharT, Allocator>::__invalidate_all_iterators()
{
}

//  Allocate space for n objects
//  throws length_error if n > max_size()
//  throws (probably bad_alloc) if memory run out
//  Precondition:  begin_ == end_ == __cap() == 0
//  Precondition:  n > 0
//  Postcondition:  capacity() == n
//  Postcondition:  size() == 0

template <class CharT, class Allocator>
void
vector<CharT, Allocator>::__vallocate(size_type n)
{
    if (n > max_size())
        this->__throw_length_error();
    n = __external_cap_to_internal(n);
    this->begin_ = __storage_traits::allocate(this->__alloc(), n);
    this->size_ = 0;
    this->__cap() = n;
}

template <class CharT, class Allocator>
void
vector<CharT, Allocator>::__vdeallocate() noexcept
{
    if (this->begin_ != nullptr)
    {
        __storage_traits::deallocate(this->__alloc(), this->begin_, __cap());
        __invalidate_all_iterators();
        this->begin_ = nullptr;
        this->size_ = this->__cap() = 0;
    }
}

template <class CharT, class Allocator>
typename vector<CharT, Allocator>::size_type
vector<CharT, Allocator>::max_size() const noexcept
{
    size_type amax = __storage_traits::max_size(__alloc());
    size_type nmax = std::numeric_limits<size_type>::max() / 2;  // end() >= begin(), always
    if (nmax / bases_per_word <= amax)
        return nmax;
    return __internal_cap_to_external(amax);
}

//  Precondition:  new_size > capacity()
template <class CharT, class Allocator>
typename vector<CharT, Allocator>::size_type
vector<CharT, Allocator>::__recommend(size_type new_size) const
{
    const size_type ms = max_size();
    if (new_size > ms)
        this->__throw_length_error();
    const size_type cap = capacity();
    if (cap >= ms / 2)
        return ms;
    return std::max(2*cap, __align_it(new_size));
}

//  Default constructs n objects starting at end_
//  Precondition:  n > 0
//  Precondition:  size() + n <= capacity()
//  Postcondition:  size() == size() + n
template <class CharT, class Allocator>
void
vector<CharT, Allocator>::__construct_at_end(size_type n, value_type x)
{
    size_type old_size = this->size_;
    this->size_ += n;
    std::fill_n(__make_iter(old_size), n, x);
}

template <class CharT, class Allocator>
template <class ForwardIterator>
std::enable_if_t
<
    __is_forward_iterator_v<ForwardIterator>
>
vector<CharT, Allocator>::__construct_at_end(ForwardIterator first, ForwardIterator last)
{
    size_type old_size = this->size_;
    this->size_ += std::distance(first, last);
    std::copy(first, last, __make_iter(old_size));
}

/**
* @brief default constructor
* 
* @tparam CharT: char_type or uint_type
* @tparam Allocator
*
*/
template <class CharT, class Allocator>
vector<CharT, Allocator>::vector()
noexcept(std::is_nothrow_default_constructible_v<allocator_type>)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_()
{
}

/**
* @brief a constructor allocates with user input allocator
* 
* @tparam CharT: char_type or uint_type
* @tparam Allocator
*
* @param a - allocator to use for all memory allocations of this container
*
*/
template <class CharT, class Allocator>
vector<CharT, Allocator>::vector(const allocator_type& a) noexcept
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_(0, static_cast<__storage_allocator>(a))
{
}

/**
* @brief a constructor that allocates size_type*n spaces
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator
*
* @param n - the size of the container
*
*/
template <class CharT, class Allocator>
vector<CharT, Allocator>::vector(size_type n)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_()
{
    if (n > 0)
    {
        __vallocate(n);
        __construct_at_end(n, 0);
    }
}

/**
* @brief a constructor that allocates size_type*n spaces with user input allocator
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param n - the size of the container
* @param a - allocator to use for all memory allocations of this container
*
*/
template <class CharT, class Allocator>
vector<CharT, Allocator>::vector(size_type n, const allocator_type& a)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_(0, static_cast<__storage_allocator>(a))
{
    if (n > 0)
    {
        __vallocate(n);
        __construct_at_end(n, 0);
    }
}

/**
* @brief a constructor that allocates size_type*n spaces with value_type x
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param n - the size of the container
* @param x - the value to initialize elements of the container with
*
*/
template <class CharT, class Allocator>
vector<CharT, Allocator>::vector(size_type n, const value_type& x)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_()
{
    if (n > 0)
    {
        __vallocate(n);
        __construct_at_end(n, x);
    }
}

/**
* @brief a constructor that allocates size_type*n spaces with value_type x and user input allocator
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param n - the size of the container
* @param x - the value to initialize elements of the container with
* @param first, last - the range to copy the elements from
*
*/
template <class CharT, class Allocator>
vector<CharT, Allocator>::vector(size_type n, const value_type& x, const allocator_type& a)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_(0, static_cast<__storage_allocator>(a))
{
    if (n > 0)
    {
        __vallocate(n);
        __construct_at_end(n, x);
    }
}

/**
* @brief a constructor that initializes with InputIterator
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param first, last - the range to copy the elements from
*
*/
template <class CharT, class Allocator>
template <class InputIterator>
vector<CharT, Allocator>::vector(InputIterator first, InputIterator last,
                                 std::enable_if_t<__is_input_iterator_v  <InputIterator> &&
                                                 !__is_forward_iterator_v<InputIterator>>*)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_()
{
    try
    {
        for (; first != last; ++first)
            push_back(*first);
    }
    catch (...)
    {
        if (begin_ != nullptr)
            __storage_traits::deallocate(__alloc(), begin_, __cap());
        __invalidate_all_iterators();
        throw;
    }
}

/**
* @brief a constructor that initializes with InputIterator and allocates size with user input allocator
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param first, last - the range to copy the elements from
* @param a - allocator to use for all memory allocations of this container
*
*/
template <class CharT, class Allocator>
template <class InputIterator>
vector<CharT, Allocator>::vector(InputIterator first, InputIterator last, const allocator_type& a,
                                 std::enable_if_t<__is_input_iterator_v  <InputIterator> &&
                                                 !__is_forward_iterator_v<InputIterator>>*)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_(0, static_cast<__storage_allocator>(a))
{
    try
    {
        for (; first != last; ++first)
            push_back(*first);
    }
    catch (...)
    {
        if (begin_ != nullptr)
            __storage_traits::deallocate(__alloc(), begin_, __cap());
        __invalidate_all_iterators();
        throw;
    }
}

/**
* @brief a constructor that initializes with ForwardIterator
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param first, last - the range to copy the elements from
*
*/
template <class CharT, class Allocator>
template <class ForwardIterator>
vector<CharT, Allocator>::vector(ForwardIterator first, ForwardIterator last,
                                 std::enable_if_t<__is_forward_iterator_v<ForwardIterator>>*)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_()
{
    size_type n = static_cast<size_type>(std::distance(first, last));
    if (n > 0)
    {
        __vallocate(n);
        __construct_at_end(first, last);
    }
}

/**
* @brief a constructor that initializes with ForwardIterator and allocates size with user input allocator
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param first, last - the range to copy the elements from
* @param allocator to use for all memory allocations of this container
*
*/
template <class CharT, class Allocator>
template <class ForwardIterator>
vector<CharT, Allocator>::vector(ForwardIterator first, ForwardIterator last, const allocator_type& a,
                                 std::enable_if_t<__is_forward_iterator_v<ForwardIterator>>*)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_(0, static_cast<__storage_allocator>(a))
{
    size_type n = static_cast<size_type>(std::distance(first, last));
    if (n > 0)
    {
        __vallocate(n);
        __construct_at_end(first, last);
    }
}

/**
* @brief a constructor that initializes with initializer_list
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param il - initializer list to initialize the elements of the container with
*
*/
template <class CharT, class Allocator>
vector<CharT, Allocator>::vector(std::initializer_list<value_type> il)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_()
{
    size_type n = static_cast<size_type>(il.size());
    if (n > 0)
    {
        __vallocate(n);
        __construct_at_end(il.begin(), il.end());
    }
}

/**
* @brief a constructor that initializes with initializer_list and allocates with user input allocator
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param il - initializer list to initialize the elements of the container with
* @param a - allocator to use for all memory allocations of this container
*
*/
template <class CharT, class Allocator>
vector<CharT, Allocator>::vector(std::initializer_list<value_type> il, const allocator_type& a)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_(0, static_cast<__storage_allocator>(a))
{
    size_type n = static_cast<size_type>(il.size());
    if (n > 0)
    {
        __vallocate(n);
        __construct_at_end(il.begin(), il.end());
    }
}

/**
* @brief destructor
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
*/
template <class CharT, class Allocator>
vector<CharT, Allocator>::~vector()
{
    if (begin_ != nullptr)
        __storage_traits::deallocate(__alloc(), begin_, __cap());
    __invalidate_all_iterators();
}

template <class CharT, class Allocator>
vector<CharT, Allocator>::vector(const vector& v)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_(0, __storage_traits::select_on_container_copy_construction(v.__alloc()))
{
    if (v.size() > 0)
    {
        __vallocate(v.size());
        __construct_at_end(v.begin(), v.end());
    }
}

/**
* @brief a copy constructor that copies from another container and allocates with user input allocator
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param v - the container to be copied
* @param a - allocator to use for all memory allocations of this container
*
*/
template <class CharT, class Allocator>
vector<CharT, Allocator>::vector(const vector& v, const allocator_type& a)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_(0, a)
{
    if (v.size() > 0)
    {
        __vallocate(v.size());
        __construct_at_end(v.begin(), v.end());
    }
}

/**
* @brief a copy assignment copies from another const container
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param v - the container to be copied
*
*/
template <class CharT, class Allocator>
vector<CharT, Allocator>&
vector<CharT, Allocator>::operator=(const vector& v)
{
    if (this != &v)
    {
        __copy_assign_alloc(v);
        if (v.size_)
        {
            if (v.size_ > capacity())
            {
                __vdeallocate();
                __vallocate(v.size_);
            }
            std::copy(v.begin_, v.begin_ + __external_cap_to_internal(v.size_), begin_);
        }
        size_ = v.size_;
    }
    return *this;
}

/**
* @brief a move constructor that moves another container
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param v - the container to be moved
*
*/
template <class CharT, class Allocator>
vector<CharT, Allocator>::vector(vector&& v) noexcept
    : begin_    (v.begin_),
      size_     (v.size_),
      cap_alloc_(std::move(v.cap_alloc_))
{
    v.begin_ = nullptr;
    v.size_ = 0;
    v.__cap() = 0;
}

/**
* @brief a move constructor that moves container and allocates with user input allocator
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param v - the container to be moved
* @param a - allocator to use for all memory allocations of this container
*
*/
template <class CharT, class Allocator>
vector<CharT, Allocator>::vector(vector&& v, const allocator_type& a)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_(0, a)
{
    if (a == allocator_type(v.__alloc()))
    {
        this->begin_ = v.begin_;
        this->size_ = v.size_;
        this->__cap() = v.__cap();
        v.begin_ = nullptr;
        v.__cap() = v.size_ = 0;
    }
    else if (v.size() > 0)
    {
        __vallocate(v.size());
        __construct_at_end(v.begin(), v.end());
    }
}

/**
* @brief a move assignment that moves another container
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param - the container to be moved
*
*/
template <class CharT, class Allocator>
vector<CharT, Allocator>&
vector<CharT, Allocator>::operator=(vector&& v)
noexcept(__is_noexcept_move_assign_container_v<Allocator, __alloc_traits>)
{
    __move_assign(v, std::bool_constant<
            __storage_traits::propagate_on_container_move_assignment::value>());
    return *this;
}

template <class CharT, class Allocator>
void
vector<CharT, Allocator>::__move_assign(vector& c, std::false_type)
{
    if (__alloc() != c.__alloc())
        assign(c.begin(), c.end());
    else
        __move_assign(c, std::true_type());
}

template <class CharT, class Allocator>
void
vector<CharT, Allocator>::__move_assign(vector& c, std::true_type)
noexcept(std::is_nothrow_move_assignable_v<allocator_type>)
{
    __vdeallocate();
    __move_assign_alloc(c);
    this->begin_ = c.begin_;
    this->size_ = c.size_;
    this->__cap() = c.__cap();
    c.begin_ = nullptr;
    c.__cap() = c.size_ = 0;
}

/**
* @brief an assign constructor that initializes with value_type x
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param n - the new size of the container
* @param x - the value to initialize elements of the container with
*
*/
template <class CharT, class Allocator>
void
vector<CharT, Allocator>::assign(size_type n, const value_type& x)
{
    size_ = 0;
    if (n > 0)
    {
        size_type c = capacity();
        if (n <= c)
            size_ = n;
        else
        {
            vector v(__alloc());
            v.reserve(__recommend(n));
            v.size_ = n;
            swap(v);
        }
        std::fill_n(begin(), n, x);
    }
    __invalidate_all_iterators();
}

/**
* @brief a assign constructor that initializes with InputIterator
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param first, last - 	the range to copy the elements from
*
*/
template <class CharT, class Allocator>
template <class InputIterator>
std::enable_if_t
<
    __is_input_iterator_v  <InputIterator> &&
   !__is_forward_iterator_v<InputIterator>
>
vector<CharT, Allocator>::assign(InputIterator first, InputIterator last)
{
    clear();
    for (; first != last; ++first)
        push_back(*first);
}

/**
* @brief a assign constructor that initializes with ForwardIterator
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param first, last - 	the range to copy the elements from
*
*/
template <class CharT, class Allocator>
template <class ForwardIterator>
std::enable_if_t
<
    __is_forward_iterator_v<ForwardIterator>
>
vector<CharT, Allocator>::assign(ForwardIterator first, ForwardIterator last)
{
    clear();
    difference_type ns = std::distance(first, last);
    assert(ns >= 0 && "invalid range specified");
    const size_t n = static_cast<size_type>(ns);
    if (n)
    {
        if (n > capacity())
        {
            __vdeallocate();
            __vallocate(n);
        }
        __construct_at_end(first, last);
    }
}

/**
* @brief reserve size_type*n spaces
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param n - new capacity of the vector
*
*/
template <class CharT, class Allocator>
void
vector<CharT, Allocator>::reserve(size_type n)
{
    if (n > capacity())
    {
        vector v(this->__alloc());
        v.__vallocate(n);
        v.__construct_at_end(this->begin(), this->end());
        swap(v);
        __invalidate_all_iterators();
    }
}

/**
* @brief requests the removal of unused capacity
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
*/
template <class CharT, class Allocator>
void
vector<CharT, Allocator>::shrink_to_fit() noexcept
{
    if (__external_cap_to_internal(size()) > __cap())
    {
        try
        {
            vector(*this, allocator_type(__alloc())).swap(*this);
        }
        catch (...)
        {
        }
    }
}

/**
* @brief access specified element with bounds checking
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param n - position of the element to return
*
* @return reference to the requested element
*
*/
template <class CharT, class Allocator>
typename vector<CharT, Allocator>::reference
vector<CharT, Allocator>::at(size_type n)
{
    if (n >= size())
        this->__throw_out_of_range();
    return (*this)[n];
}

/**
* @brief access specified const element with bounds checking
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param n - position of the element to return
*
* @return reference to the const requested element
*
*/
template <class CharT, class Allocator>
typename vector<CharT, Allocator>::const_reference
vector<CharT, Allocator>::at(size_type n) const
{
    if (n >= size())
        this->__throw_out_of_range();
    return (*this)[n];
}


/**
* @brief adds an element to the end
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param x - the element to be added
*
*/
template <class CharT, class Allocator>
void
vector<CharT, Allocator>::push_back(const value_type& x)
{
    if (this->size_ == this->capacity())
        reserve(__recommend(this->size_ + 1));
    ++this->size_;
    back() = x;
}

/**
* @brief inserts n elements at specific position with value_type x
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param position - const iterator. For example, v.begin()+3
* @param x - the value to be inserted.
*
* @return iterator pointing to the inserted value
*/
template <class CharT, class Allocator>
typename vector<CharT, Allocator>::iterator
vector<CharT, Allocator>::insert(const_iterator position, const value_type& x)
{
    iterator r;
    if (size() < capacity())
    {
        const_iterator old_end = end();
        ++size_;
        std::copy_backward(position, old_end, end());
        r = __const_iterator_cast(position);
    }
    else
    {
        vector v(__alloc());
        v.reserve(__recommend(size_ + 1));
        v.size_ = size_ + 1;
        r = std::copy(cbegin(), position, v.begin());
        std::copy_backward(position, cend(), v.end());
        swap(v);
    }
    *r = x;
    return r;
}


/**
* @brief inserts n elements at specific position with value_type x
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param position: const iterator. For example, v.begin()+3
* @param n - the number of elements to be inserted
* @param x - the value to be inserted.
*
* @return iterator pointing to the inserted value
*/
template <class CharT, class Allocator>
typename vector<CharT, Allocator>::iterator
vector<CharT, Allocator>::insert(const_iterator position, size_type n, const value_type& x)
{
    iterator r;
    size_type c = capacity();
    if (n <= c && size() <= c - n)
    {
        const_iterator old_end = end();
        size_ += n;
        std::copy_backward(position, old_end, end());
        r = __const_iterator_cast(position);
    }
    else
    {
        vector v(__alloc());
        v.reserve(__recommend(size_ + n));
        v.size_ = size_ + n;
        r = std::copy(cbegin(), position, v.begin());
        std::copy_backward(position, cend(), v.end());
        swap(v);
    }
    std::fill_n(r, n, x);
    return r;
}

/**
* @brief inserts n elements at specific position with value_type x
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param position - const iterator. For example, v.begin()+3
* @param first, last - insert elements start from range [first, last) position
* @param x - the value to be inserted
*
* @return iterator pointing to the inserted value
*/
template <class CharT, class Allocator>
template <class InputIterator>
std::enable_if_t
<
    __is_input_iterator_v  <InputIterator> &&
   !__is_forward_iterator_v<InputIterator>,
    typename vector<CharT, Allocator>::iterator
>
vector<CharT, Allocator>::insert(const_iterator position, InputIterator first, InputIterator last)
{
    difference_type off = position - begin();
    iterator p = __const_iterator_cast(position);
    iterator old_end = end();
    for (; size() != capacity() && first != last; ++first)
    {
        ++this->size_;
        back() = *first;
    }
    vector v(__alloc());
    if (first != last)
    {
        try
        {
            v.assign(first, last);
            difference_type old_size = static_cast<difference_type>(old_end - begin());
            difference_type old_p = p - begin();
            reserve(__recommend(size() + v.size()));
            p = begin() + old_p;
            old_end = begin() + old_size;
        }
        catch (...)
        {
            erase(old_end, end());
            throw;
        }
    }
    p = std::rotate(p, old_end, end());
    insert(p, v.begin(), v.end());
    return begin() + off;
}

/**
* @brief inserts n elements at specific position with value_type x
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param position - const iterator. For example, v.begin()+3
* @param first, last - insert elements start from range [first, last) position
* @param x - the value to be inserted
*
* @return iterator pointing to the inserted value
*/
template <class CharT, class Allocator>
template <class ForwardIterator>
std::enable_if_t
<
    __is_forward_iterator_v<ForwardIterator>,
    typename vector<CharT, Allocator>::iterator
>
vector<CharT, Allocator>::insert(const_iterator position, ForwardIterator first, ForwardIterator last)
{
    const difference_type n_signed = std::distance(first, last);
    assert(n_signed >= 0 && "invalid range specified");
    const size_type n = static_cast<size_type>(n_signed);
    iterator r;
    size_type c = capacity();
    if (n <= c && size() <= c - n)
    {
        const_iterator old_end = end();
        size_ += n;
        std::copy_backward(position, old_end, end());
        r = __const_iterator_cast(position);
    }
    else
    {
        vector v(__alloc());
        v.reserve(__recommend(size_ + n));
        v.size_ = size_ + n;
        r = std::copy(cbegin(), position, v.begin());
        std::copy_backward(position, cend(), v.end());
        swap(v);
    }
    std::copy(first, last, r);
    return r;
}

/**
* @brief erases element at certain position
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param position - the position of value to be erased 
*
* @return iterator pointing to the last removed element
*/
template <class CharT, class Allocator>
typename vector<CharT, Allocator>::iterator
vector<CharT, Allocator>::erase(const_iterator position)
{
    iterator r = __const_iterator_cast(position);
    std::copy(position + 1, this->cend(), r);
    --size_;
    return r;
}

/**
* @brief erases elements between two iterator
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
* 
* @param first, end - erased value from range [first, last) position
*
* @return iterator pointing to the last removed element
*/
template <class CharT, class Allocator>
typename vector<CharT, Allocator>::iterator
vector<CharT, Allocator>::erase(const_iterator first, const_iterator last)
{
    iterator r = __const_iterator_cast(first);
    difference_type d = last - first;
    std::copy(last, this->cend(), r);
    size_ -= d;
    return r;
}

/**
* @brief swaps the content of two vector
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
* 
* @param x - the container to be swapped
*
*/
template <class CharT, class Allocator>
void
vector<CharT, Allocator>::swap(vector& x) noexcept
{
    std::swap(this->begin_, x.begin_);
    std::swap(this->size_, x.size_);
    std::swap(this->__cap(), x.__cap());
    if constexpr (std::bool_constant<__alloc_traits::propagate_on_container_swap::value>::value)
        std::swap(this->__alloc(), x.__alloc());
}

/**
* @brief changes the number of elements stored
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param sz - new size of the container
* @param x - the value to initialize the new elements with
*
*/
template <class CharT, class Allocator>
void
vector<CharT, Allocator>::resize(size_type sz, value_type x)
{
    size_type cs = size();
    if (cs < sz)
    {
        iterator r;
        size_type c = capacity();
        size_type n = sz - cs;
        if (n <= c && cs <= c - n)
        {
            r = end();
            size_ += n;
        }
        else
        {
            vector v(__alloc());
            v.reserve(__recommend(size_ + n));
            v.size_ = size_ + n;
            r = std::copy(cbegin(), cend(), v.begin());
            swap(v);
        }
        std::fill_n(r, n, x);
    }
    else
        size_ = sz;
}

/**
* @brief changes the bases to complementary bases
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
*/
template <class CharT, class Allocator>
void
vector<CharT, Allocator>::flip() noexcept
{
    // do middle whole words
    size_type n = size_;
    __storage_pointer p = begin_;
    for (; n >= bases_per_word; ++p, n -= bases_per_word)
        *p = ~*p;

    // do last partial word
    if (n > 0)
    {
        __storage_type m = ~__storage_type(0) >> (bases_per_word - n);
        __storage_type b = *p & m;
        *p &= ~m;
        *p |= ~b & m;
    }
}

/**
* @brief changes the bases to complementary bases
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
*/
template <class CharT, class Allocator>
bool
vector<CharT, Allocator>::__invariants() const
{
    if (this->begin_ == nullptr)
    {
        if (this->size_ != 0 || this->__cap() != 0)
            return false;
    }
    else
    {
        if (this->__cap() == 0)
            return false;
        if (this->size_ > this->capacity())
            return false;
    }
    return true;
}

/**
* @brief compares whether two vectors are identical or not
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param x, y - lhs and rhs container to be compared with
*
*/
template <class CharT, class Allocator>
bool
operator==(const vector<CharT, Allocator>& x, const vector<CharT, Allocator>& y)
{
    const typename vector<CharT, Allocator>::size_type sz = x.size();
    return sz == y.size() && std::equal(x.begin(), x.end(), y.begin());
}

template <class CharT, class Allocator>
bool
operator!=(const vector<CharT, Allocator>& x, const vector<CharT, Allocator>& y)
{
    return !(x == y);
}

/**
* @brief compares whether two vectors are identical or not
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param x, y - lhs and rhs container to be compared with
*
*/
template <class CharT, class Allocator>
bool
operator< (const vector<CharT, Allocator>& x, const vector<CharT, Allocator>& y)
{
    return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
}

/**
* @brief compares whether lhs > rhs
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param x, y - lhs and rhs container to be compared with
*
*/
template <class CharT, class Allocator>
bool
operator> (const vector<CharT, Allocator>& x, const vector<CharT, Allocator>& y)
{
    return y < x;
}

/**
* @brief compares whether lhs >= rhs
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param x, y - lhs and rhs container to be compared with
*
*/
template <class CharT, class Allocator>
bool
operator>=(const vector<CharT, Allocator>& x, const vector<CharT, Allocator>& y)
{
    return !(x < y);
}

/**
* @brief compares whether lhs <= rhs
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param x, y - lhs and rhs container to be compared with
*
*/
template <class CharT, class Allocator>
bool
operator<=(const vector<CharT, Allocator>& x, const vector<CharT, Allocator>& y)
{
    return !(y < x);
}

/**
* @brief switches the content of lhs and rhs
* 
* @tparam CharT - char_type or uint_type
* @tparam Allocator 
*
* @param x, y - lhs and rhs container to be compared with
*
*/
template <class CharT, class Allocator>
void
swap(vector<CharT, Allocator>& x, vector<CharT, Allocator>& y)
noexcept(noexcept(x.swap(y)))
{
    x.swap(y);
}

// for value_type = char_type:
/**
* @brief a constructor that initializes with value s
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
* @param s - const value_type pointer used as source to initialize
*
*/
template <class CharT, class Allocator>
template <class charT>
vector<CharT, Allocator>::vector(const value_type* s,
                                 std::enable_if_t<std::is_same_v<charT, char_type>>*)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_()
{
    assert(s != nullptr && "vector(const char*) detected nullptr");
    std::string_view sv(s);
    size_type n = sv.size();
    if (n > 0)
    {
        __vallocate(n);
        __construct_at_end(sv.begin(), sv.end());
    }
}

/**
* @brief a constructor that initializes with value_type s and allocates with user input allocator
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
* @param s - const value_type pointer used to initialize
* @param a - allocator to use for all memory allocations of this container
*
*/
template <class CharT, class Allocator>
template <class charT>
vector<CharT, Allocator>::vector(const value_type* s, const allocator_type& a,
                                 std::enable_if_t<std::is_same_v<charT, char_type>>*)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_(0, static_cast<__storage_allocator>(a))
{
    assert(s != nullptr && "vector(const char*, allocator) detected nullptr");
    std::string_view sv(s);
    size_type n = sv.size();
    if (n > 0)
    {
        __vallocate(n);
        __construct_at_end(sv.begin(), sv.end());
    }
}

/**
* @brief a constructor that initializes with string
* 
* @tparam CharT: char_type
* @tparam Allocator 
*
* @param s - string used to initialize
*
*/
template <class CharT, class Allocator>
template <class charT>
vector<CharT, Allocator>::vector(const std::string& s,
                                 std::enable_if_t<std::is_same_v<charT, char_type>>*)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_()
{
    size_type n = s.size();
    if (n > 0)
    {
        __vallocate(n);
        __construct_at_end(s.begin(), s.end());
    }
}

/**
* @brief a constructor that initializes with const value_type s and allocates with user input allocator
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
* @param s - string used to initialize
* @param a - allocator to use for all memory allocations of this container
*
*/
template <class CharT, class Allocator>
template <class charT>
vector<CharT, Allocator>::vector(const std::string& s, const allocator_type& a,
                                 std::enable_if_t<std::is_same_v<charT, char_type>>*)
    : begin_    (nullptr),
      size_     (0),
      cap_alloc_(0, static_cast<__storage_allocator>(a))
{
    size_type n = s.size();
    if (n > 0)
    {
        __vallocate(n);
        __construct_at_end(s.begin(), s.end());
    }
}

/**
* @brief an assignment that assign the container with value_type s
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
* @param s - const value_type pointer used to initialize
*
*/
template <class CharT, class Allocator>
template <class charT>
std::enable_if_t<std::is_same_v<charT, char_type>, vector<CharT, Allocator>&>
vector<CharT, Allocator>::operator=(const value_type* s)
{
    assert(s != nullptr && "vector::operator= received nullptr");
    std::string_view sv(s);
    assign(sv.begin(), sv.end());
    return *this;
}

/**
* @brief an assignment that assign the container with string s
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
* @param s - string used to initialize
*
*/
template <class CharT, class Allocator>
template <class charT>
std::enable_if_t<std::is_same_v<charT, char_type>, vector<CharT, Allocator>&>
vector<CharT, Allocator>::operator=(const std::string& s)
{
    assign(s.begin(), s.end());
    return *this;
}

/**
* @brief an assignment that assign the container with value_type c
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
* @param c - value_type used to initialize
*
*/
template <class CharT, class Allocator>
template <class charT>
std::enable_if_t<std::is_same_v<charT, char_type>, vector<CharT, Allocator>&>
vector<CharT, Allocator>::operator=(value_type c)
{
    clear();
    push_back(c);
    return *this;
}

/**
* @brief concatenates two containers
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
* @param v - container to be concatenated
*
*/
template <class CharT, class Allocator>
template <class charT>
std::enable_if_t<std::is_same_v<charT, char_type>, vector<CharT, Allocator>&>
vector<CharT, Allocator>::operator+=(const vector& v)
{
    for (value_type c : v)
        push_back(c);
    return *this;
}

/**
* @brief appends a string to the container
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
* @param s - string to be appended
*
*/
template <class CharT, class Allocator>
template <class charT>
std::enable_if_t<std::is_same_v<charT, char_type>, vector<CharT, Allocator>&>
vector<CharT, Allocator>::operator+=(const std::string& s)
{
    for (value_type c : s)
        push_back(c);
    return *this;
}

/**
* @brief appends a value_type pointer to the container
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
* @param s - const value_type pointer to be appended
*
*/
template <class CharT, class Allocator>
template <class charT>
std::enable_if_t<std::is_same_v<charT, char_type>, vector<CharT, Allocator>&>
vector<CharT, Allocator>::operator+=(const value_type* s)
{
    assert(s != nullptr && "vector::operator+= received nullptr");
    std::string_view sv(s);
    for (value_type c : sv)
        push_back(c);
    return *this;
}

/**
* @brief appends a value_type to the container
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
* @param c - value_type to be appended
*
*/
template <class CharT, class Allocator>
template <class charT>
std::enable_if_t<std::is_same_v<charT, char_type>, vector<CharT, Allocator>&>
vector<CharT, Allocator>::operator+=(value_type c)
{
    push_back(c);
    return *this;
}

/**
* @brief appends a initializer_list to the container
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
* @param il - initializer_list to be appended
*
*/
template <class CharT, class Allocator>
template <class charT>
std::enable_if_t<std::is_same_v<charT, char_type>, vector<CharT, Allocator>&>
vector<CharT, Allocator>::operator+=(std::initializer_list<value_type> il)
{
    for (value_type c : il)
        push_back(c);
    return *this;
}

/**
* @brief compares whether charT lhs == container rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator==(const CharT* lhs, const vector<CharT, Allocator>& rhs) noexcept
{
    std::string_view sv(lhs);
    const std::string_view::size_type sz = sv.size();
    return sz == rhs.size() && std::equal(sv.begin(), sv.end(), rhs.begin());
}

/**
* @brief compares whether container lhs == charT rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator==(const vector<CharT, Allocator>& lhs, const CharT* rhs) noexcept
{
    std::string_view sv(rhs);
    const std::string_view::size_type sz = sv.size();
    return sz == lhs.size() && std::equal(lhs.begin(), lhs.end(), sv.begin());
}

/**
* @brief compares whether string lhs == container rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator==(const std::string& lhs, const vector<CharT, Allocator>& rhs) noexcept
{
    const std::string::size_type sz = lhs.size();
    return sz == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

/**
* @brief compares whether container lhs == string rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator==(const vector<CharT, Allocator>& lhs, const std::string& rhs) noexcept
{
    const typename vector<CharT, Allocator>::size_type sz = lhs.size();
    return sz == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

/**
* @brief compares whether CharT lhs != container rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator!=(const CharT* lhs, const vector<CharT, Allocator>& rhs) noexcept
{
    return !(lhs == rhs);
}

/**
* @brief compares whether container lhs != CharT rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator!=(const vector<CharT, Allocator>& lhs, const CharT* rhs) noexcept
{
    return !(lhs == rhs);
}

/**
* @brief compares whether string lhs != container rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator!=(const std::string& lhs, const vector<CharT, Allocator>& rhs) noexcept
{
    return !(lhs == rhs);
}

/**
* @brief compares whether container lhs != string rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator!=(const vector<CharT, Allocator>& lhs, const std::string& rhs) noexcept
{
    return !(lhs == rhs);
}

/**
* @brief compares whether container lhs < CharT rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator< (const vector<CharT, Allocator>& lhs, const CharT* rhs) noexcept
{
    std::string_view sv(rhs);
    return std::lexicographical_compare(lhs.begin(), lhs.end(), sv.begin(), sv.end());
}

/**
* @brief compares whether charT lhs < container rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator< (const CharT* lhs, const vector<CharT, Allocator>& rhs) noexcept
{
    std::string_view sv(lhs);
    return std::lexicographical_compare(sv.begin(), sv.end(), rhs.begin(), rhs.end());
}

/**
* @brief compares whether container lhs < string rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator< (const vector<CharT, Allocator>& lhs, const std::string& rhs) noexcept
{
    return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

/**
* @brief compares whether string lhs < container rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator< (const std::string& lhs, const vector<CharT, Allocator>& rhs) noexcept
{
    return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

/**
* @brief compares whether container lhs > CharT rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator> (const vector<CharT, Allocator>& lhs, const CharT* rhs) noexcept
{
    return rhs < lhs;
}

/**
* @brief compares whether CharT lhs > container rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator> (const CharT* lhs, const vector<CharT, Allocator>& rhs) noexcept
{
    return rhs < lhs;
}

/**
* @brief compares whether container lhs > string rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator> (const vector<CharT, Allocator>& lhs, const std::string& rhs) noexcept
{
    return rhs < lhs;
}

/**
* @brief compares whether string lhs > container rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator> (const std::string& lhs, const vector<CharT, Allocator>& rhs) noexcept
{
    return rhs < lhs;
}

/**
* @brief compares whether container lhs <= CharT rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator<=(const vector<CharT, Allocator>& lhs, const CharT* rhs) noexcept
{
    return !(rhs < lhs);
}

/**
* @brief compares whether CharT lhs <= container rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator<=(const CharT* lhs, const vector<CharT, Allocator>& rhs) noexcept
{
    return !(rhs < lhs);
}

/**
* @brief compares whether container lhs <= string rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator<=(const vector<CharT, Allocator>& lhs, const std::string& rhs) noexcept
{
    return !(rhs < lhs);
}

/**
* @brief compares whether string lhs <= container rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator<=(const std::string& lhs, const vector<CharT, Allocator>& rhs) noexcept
{
    return !(rhs < lhs);
}

/**
* @brief compares whether container lhs >= CharT rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator>=(const vector<CharT, Allocator>& lhs, const CharT* rhs) noexcept
{
    return !(lhs < rhs);
}

/**
* @brief compares whether CharT lhs >= container rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator>=(const CharT* lhs, const vector<CharT, Allocator>& rhs) noexcept
{
    return !(lhs < rhs);
}

/**
* @brief compares whether container lhs >= string rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator>=(const vector<CharT, Allocator>& lhs, const std::string& rhs) noexcept
{
    return !(lhs < rhs);
}

/**
* @brief compares whether string lhs >= container rhs
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, bool>
operator>=(const std::string& lhs, const vector<CharT, Allocator>& rhs) noexcept
{
    return !(lhs < rhs);
}

/**
* @brief performs stream output on container
* 
* @tparam CharT - char_type
* @tparam Allocator 
*
*/
template<class CharT, class Allocator>
std::enable_if_t<std::is_same_v<CharT, char_type>, std::ostream&>
operator<<(std::ostream& os, const vector<CharT, Allocator>& v)
{
    for (CharT c : v)
        os << c;
    return os;
}

using Sequence = vector<char_type>;

inline
Sequence operator "" _seq(const char* str, std::size_t len)
{
    return Sequence(str);
}

}

