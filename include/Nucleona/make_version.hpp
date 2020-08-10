#pragma once
#include <Nucleona/algo/split.hpp>
#include <Nucleona/language.hpp>
namespace nucleona {
    
struct Version {
private:
    static auto parse_micro_code_and_ext_tag( const std::string& str) {
        auto mic_ext_tag = nucleona::algo::split(str, "-");
        auto micro_code = std::stoi(mic_ext_tag.at(0));
        std::string ext_tag;
        if(mic_ext_tag.size() > 1) {
            ext_tag = mic_ext_tag.at(1);
        }
        return std::make_tuple(micro_code, ext_tag);
    }
public:
    Version() = default;
    Version(int _major_code, int _minor_code, const std::string& _micro_code_ext_tag )
    : major_code(_major_code)
    , minor_code(_minor_code)
    {
        std::tie(micro_code, ext_tag) = parse_micro_code_and_ext_tag(_micro_code_ext_tag);
    }
    bool operator==(const Version& v) const {
        return 
            major_code == v.major_code &&
            minor_code == v.minor_code &&
            micro_code == v.micro_code
        ;
    }
    bool operator<=( const Version& v) const {
        if( major_code == v.major_code ) {
            if(minor_code == v.minor_code) {
                if( micro_code == v.micro_code ) {
                    return true;
                } else return micro_code < v.micro_code;
            } else return minor_code < v.minor_code;
        } else return major_code < v.major_code;
    }
    bool operator>( const Version& v) const {
        return !(operator<=(v));
    }
    std::string to_string() const {
        auto res = 
            std::to_string(major_code) + "." 
            + std::to_string(minor_code) + "."
            + std::to_string(micro_code)
        ;
        if(!ext_tag.empty()) {
            res += ("-" + ext_tag);
        }
        return res;
    }
    static Version parse( const std::string& str) {
        Version vsn;
        auto version_str_vec = nucleona::algo::split(str, ".");

        vsn.major_code = std::stoi(version_str_vec.at(0));
        vsn.minor_code = std::stoi(version_str_vec.at(1));
        std::tie(
            vsn.micro_code, 
            vsn.ext_tag
        ) = parse_micro_code_and_ext_tag(version_str_vec.at(2));
        return vsn;
    }
    int major_code;
    int minor_code;
    int micro_code;
    std::string ext_tag;
};

}

#define MAKE_VERSION(name,major,minor,micro) \
nucleona::Version name() { \
    return nucleona::Version( \
        major, minor, UNWRAP_SYM_STR(micro) \
    ); \
}