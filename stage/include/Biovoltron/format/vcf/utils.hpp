#pragma once
#include <vector>
#include <Biovoltron/utils/string_view.hpp>
namespace biovoltron::format::vcf{

/**
 * @brief Refer from boost::split, but for string_view and can
 * ignore anything inside the quote
 *
 * @param results A vector of string_view to store splitted 
 * results
 * @param str A string_view which want to split
 * @param delim A string_view used as delimeter
 *
 * Reference from boost::split(). But we don't need the 
 * function of boost::is_any_of because we always split inputs by 
 * only a char. And since segments surround by quote may contain 
 * char we used to split, we have to ignore these char when splitting 
 * to reach the expect results as we thought
 *
 * Time Complexity: O(n)<br>
 *	\e n the length of input string_view
 */
void split_with_quote_handle(
    std::vector<utils::StringView>& results,
    const utils::StringView& str,
    const utils::StringView& delim
);
}
