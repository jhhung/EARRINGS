#pragma once
#include <utility>

namespace nucleona{ 

template<class CHAR>
class ConstStringProto
{
public:
    using value_type = CHAR;
    using size_type = std::size_t;
    using iterator = CHAR*;
    using const_iterator = const CHAR*;
    constexpr ConstStringProto( const CHAR* str )
    : str_ ( str )
    {}
    constexpr size_type size() const
    {
        size_type i(0); 
        for( i = 0; str_[i] != '\0'; i ++ );
        return i;
    }
    constexpr decltype(auto) c_str() const
    {
        return str_;
    }
    constexpr auto& at( size_type i ) const 
    {
        return *(str_ + i);
    }
    constexpr auto& back() const 
    {
        return *(str_ + size() - 1);
    }
    constexpr const_iterator begin() const
    {
        return str_;
    }
    constexpr const_iterator end() const
    {
        return str_ + size();
    }
private:
    const CHAR* str_;
};
template<class CHAR>
constexpr auto make_const_string(const CHAR* cstr )
{
    return ConstStringProto<CHAR>(cstr);
}

template<class T, class CHAR>
decltype(auto) operator<<(T&& o, const ConstStringProto<CHAR>& str )
{
    return o << str.c_str();
}
template<class T, class CHAR>
decltype(auto) operator==(const ConstStringProto<CHAR>& str, T&& o )
{
    return o == str.c_str();
}
template<class T, class CHAR>
decltype(auto) operator==(T&& o, const ConstStringProto<CHAR>& str)
{
    return operator==( str, FWD(o) );
}
using ConstString = ConstStringProto<char>;

}
