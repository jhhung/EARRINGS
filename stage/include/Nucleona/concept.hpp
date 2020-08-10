/**
 * @file    Nucleona/concept.hpp
 * @author  Chia-Hua Chang
 * @brief   The concept check library
 * @details This is a compiling time type check system implement with SFINAE
 * feature of C++ language.
 * The main idea is use decltype(...) to create a type inference path,
 * and if the path doesn't satisfied for compiler, then the function or type
 * will be removed from candidate symbol list in compiling time.
 */
#pragma once
#include <Nucleona/concept/core.hpp>
#include <Nucleona/concept/iterator.hpp>
#include <Nucleona/concept/same.hpp>
#include <Nucleona/concept/lvalue_reference.hpp>
#include <Nucleona/concept/rvalue_reference.hpp>
#include <Nucleona/concept/reference.hpp>

namespace nucleona {
}
