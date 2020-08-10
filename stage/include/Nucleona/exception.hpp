#pragma once
#include <string>
#include <Nucleona/language/forward.hpp>
namespace nucleona {
    template<class T>
    struct DebugException : public T
    {
        template<class... ARGS>
        DebugException(
            const std::string& file,
            int line_n,
            const std::string& func,
            ARGS&&... args 
        )
        : T(FWD(args)...)
        {
            std::stringstream ss;
            ss 
                << "exception throw at " 
                << file << ":" << line_n << ":[" << func << "]\n";
            ss 
                << "what(): " << T::what() << std::endl;
            dbg_msg_ = ss.str();
        }
        virtual const char* what() const noexcept {
            return dbg_msg_.c_str();
        }
    private:

        std::string dbg_msg_;
    };
    template<class T>
    auto make_debug_exception(
        const std::string& file,
        int line,
        const std::string& func,
        T&& e
    ){
        return DebugException<std::decay_t<T>>(file, line, func, FWD(e));
    }
}