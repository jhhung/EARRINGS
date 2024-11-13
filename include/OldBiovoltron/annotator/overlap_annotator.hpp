#pragma once
#include <vector>
#include <tuple>
#include <fstream>
#include <type_traits>
#include <OldBiovoltron/annotator/annotation.hpp>

namespace biovoltron::annotator {
template<class ALGORITHM>
class Overlap_annotator {
public:
    template<bool...>
    struct bool_pack;

    template<bool... v>
    using all_true = std::is_same<bool_pack<true, v...>, bool_pack<v..., true>>;

    Overlap_annotator() {
        algorithm_ = new ALGORITHM();
    }

    template<class ...DBs>
    void set_db(DBs &&...dbs) {
        // check if all the input container is a vector type
        static_assert(
                all_true<std::is_same<
                        std::vector<typename std::decay<DBs>::type::value_type>,
                        typename std::decay<DBs>::type>::value...>::value,
                "all the input type must be std::vector<T>");
        
        algorithm_->set_db(std::forward<DBs>(dbs)...);
    }

    template<class ...QUERIES>
    void set_query(QUERIES &&...queries) {
        static_assert(
                all_true<std::is_same<
                        std::vector<typename std::decay<QUERIES>::type::value_type>,
                        typename std::decay<QUERIES>::type
                >::value...>::value,
                "all the input type must be std::vector<T>");
        
        algorithm_->set_query(std::forward<QUERIES>(queries)...);
    }

    void reset() {
        algorithm_->reset();
    }

    auto annotate()
    -> decltype(std::declval<ALGORITHM>().annotate()) {
        return algorithm_->annotate();
    }

    ~Overlap_annotator() {}

private:
    ALGORITHM *algorithm_;
};
}
