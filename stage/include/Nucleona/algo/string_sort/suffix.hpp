#pragma once
#include <cinttypes>
namespace nucleona::algo::string_sort {
template<class String, class Idx = std::size_t>
struct Suffix {
    using value_type = char;
    Suffix() = default;
    Suffix(String& s, std::size_t i = 0)
    : str_(&s)
    , idx_(i)
    {}
    String& data() {
        return *str_;
    }
    auto& idx() {
        return idx_;
    }
    const auto& idx() const {
        return idx_;
    }
    decltype(auto) operator[](std::size_t i) const {
        return str_->operator[](idx_ + i);
    }
    decltype(auto) operator[](std::size_t i) {
        return str_->operator[](idx_ + i);
    }
    auto& output(std::ostream& os) const {
        for(std::size_t i = idx_; i < str_->size(); i ++ ) {
            os << (*str_)[i];
        }
        return os;
    }
private:
    String* str_;
    Idx idx_;
};
template<class String, class Idx>
auto make_suffix(String& str, const Idx& idx) {
    return Suffix<String>{str, idx};
}
}