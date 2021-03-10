/*
 * This file is part of HQSpre.
 *
 * Copyright 2016/17 Ralf Wimmer, Sven Reimer, Paolo Marin, Bernd Becker
 * Albert-Ludwigs-Universitaet Freiburg, Freiburg im Breisgau, Germany
 *
 * HQSpre is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HQSpre is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with HQSpre. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HQSPRE_AUX_HPP_
#define HQSPRE_AUX_HPP_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <utility>

/**
 * \file auxil.hpp
 * \brief Auxiliary macros and functions (mainly optimized set operations)
 * \author Ralf Wimmer
 * \date 2015-16
 */

/**
 * \def val_assert_msg(x, m)
 * \brief Macro for asserting conditions
 *
 * If the condition \f$x\f$ is violated the program is aborted and a
 * message \f$m\f$ is printed.
 * Using valgrind, a stack trace can be obtained (which is not
 * supported by the standard assert macro).
 */

/**
 * \def val_assert(x)
 * \brief Macro for asserting conditions
 *
 * If the condition is violated the program is aborted.
 * Using valgrind, a stack trace can be obtained (which is not
 * supported by the standard assert macro).
 */

#if defined(val_assert)
#    warning Macro val_assert already defined!
#    undef val_assert
#endif

#if defined(val_assert_msg)
#    undef val_assert_msg
#endif

#if defined(HAVE_VALGRIND) && !defined(NDEBUG)
#    include <valgrind/valgrind.h>
#    define val_assert(x)                                                                      \
        do {                                                                                   \
            if (!(x)) {                                                                        \
                VALGRIND_PRINTF_BACKTRACE("Assertion failed at %s:%i!\n", __FILE__, __LINE__); \
                assert(x);                                                                     \
            }                                                                                  \
        } while (0)

#    define val_assert_msg(x, y)                                                                              \
        do {                                                                                                  \
            if (!(x)) {                                                                                       \
                VALGRIND_PRINTF_BACKTRACE("Assertion failed at %s:%i!\nReason: %s\n", __FILE__, __LINE__, y); \
                assert(x);                                                                                    \
            }                                                                                                 \
        } while (0)
#else
#    define val_assert(x) assert(x)

#    define val_assert_msg(x, y) assert((x) && (y))
#endif

namespace hqspre {

/**
 * Represents the truth value of a variable or a formula.
 */
enum class TruthValue : int
{
    FALSE   = 0,  ///< false or unsatisfiable
    TRUE    = 1,  ///< true or satisfiable
    UNKNOWN = 2   ///< unassigned or unknown truth value
};

/**
 * \brief Inserts the contents of set 'source' into the set 'target'
 * \pre Both ranges need to be sorted.
 */
template <typename T1, typename T2>
inline void
fast_set_union(const T1& source, T2& target)
{
    auto pos = target.begin();
    for (const auto& element : source) {
        pos = target.insert(pos, element);
    }
}

/**
 * Given two sorted ranges A = [a_iter, a_end)  and B = [b_iter, b_end) without
 * duplicates, this function checks if \f$A\setminus B\f$ is empty. The result
 * is true iff \f$A\setminus B\f$ is empty.
 *
 * \param a_iter Iterator pointing to the first element of the first range
 * \param a_end Iterator pointing behind the last element of the first range
 * \param b_iter Iterator pointing to the first element of the second range
 * \param b_end Iterator pointing behind the last element of the second range
 * \return \f$A\setminus B=\emptyset\f$
 */
template <typename InputIterator1, typename InputIterator2>
bool
set_difference_empty(InputIterator1 a_iter, InputIterator1 a_end, InputIterator2 b_iter, InputIterator2 b_end)
{
    while (a_iter != a_end) {
        if (b_iter == b_end || *a_iter < *b_iter) {
            return false;
        } else if (*a_iter == *b_iter) {
            ++a_iter;
            ++b_iter;
        } else {  // *a_iter > *b_iter
            ++b_iter;
        }
    }

    return true;
}

/**
 * Given two sorted ranges A = [a_iter, a_end)  and B = [b_iter, b_end) without
 * duplicates, this function checks if \f$A\setminus B\f$ and \f$B\setminus A\f$
 * are empty. The result is a pair (e_a,e_b) such that e_a = true iff
 * \f$A\setminus B\f$ is empty and e_b = true iff \f$B\setminus A\f$ is empty.
 *
 * \param a_iter Iterator pointing to the first element of the first range
 * \param a_end Iterator pointing behind the last element of the first range
 * \param b_iter Iterator pointing to the first element of the second range
 * \param b_end Iterator pointing behind the last element of the second range
 * \return \f$(A\setminus B=\emptyset, B\setminus A=\emptyset)\f$
 */
template <typename InputIterator1, typename InputIterator2>
std::pair<bool, bool>
two_sided_difference_empty(InputIterator1 a_iter, InputIterator1 a_end, InputIterator2 b_iter, InputIterator2 b_end)
{
    bool a_minus_b_empty = true;
    bool b_minus_a_empty = true;

    while (true) {
        if (a_iter == a_end && b_iter == b_end)
            break;
        else if (a_iter == a_end) {
            b_minus_a_empty = false;
            break;
        } else if (b_iter == b_end) {
            a_minus_b_empty = false;
            break;
        } else if (*a_iter == *b_iter) {
            ++a_iter;
            ++b_iter;
        } else if (*a_iter < *b_iter) {
            ++a_iter;
            a_minus_b_empty = false;
        } else if (*a_iter > *b_iter) {
            ++b_iter;
            b_minus_a_empty = false;
        }

        if (!a_minus_b_empty && !b_minus_a_empty) break;
    }

    return std::make_pair(a_minus_b_empty, b_minus_a_empty);
}

/**
 * Given two sorted ranges A = [a_iter, a_end) and B = [b_iter, b_end) without
 * duplicates, this function computes \f$A\setminus B\f$ and \f$B\setminus A\f$
 * and stores them in the sets a_minus_b and b_minus_a, respectively.
 *
 * \param a_iter Iterator pointing ot the first element of the first range
 * \param a_end Iterator pointing behind the last element of the first range
 * \param b_iter Iterator pointing ot the first element of the second range
 * \param b_end Iterator pointing behind the last element of the second range
 * \param a_minus_b Iterator to the container in which the elements of
 * \f$A\setminus B\f$ are stored \param b_minus_a Iterator to the container in
 * which the elements of \f$B\setminus A\f$ are stored \return a pair of
 * iterators pointing to the end of the constructed ranges.
 */
template <typename InputIterator1, typename InputIterator2, typename OutputIterator1, typename OutputIterator2>
std::pair<OutputIterator1, OutputIterator2>
two_sided_difference(InputIterator1 a_iter, InputIterator1 a_end, InputIterator2 b_iter, InputIterator2 b_end,
                     OutputIterator1 a_minus_b, OutputIterator2 b_minus_a)
{
    while (true) {
        if (a_iter == a_end && b_iter == b_end)
            break;
        else if (a_iter == a_end && b_iter != b_end) {
            std::copy(b_iter, b_end, b_minus_a);
            break;
        } else if (a_iter != a_end && b_iter == b_end) {
            std::copy(a_iter, a_end, a_minus_b);
            break;
        } else if (*a_iter < *b_iter) {
            *a_minus_b = *a_iter;
            ++a_iter;
            ++a_minus_b;
        } else if (*a_iter > *b_iter) {
            *b_minus_a = *b_iter;
            ++b_iter;
            ++b_minus_a;
        } else {
            ++a_iter;
            ++b_iter;
        }  // end if
    }      // end while
    return std::make_pair(a_minus_b, b_minus_a);
}

/**
 * Given two sorted ranges A = [a_iter, a_end) and B = [minus_iter, minus_end)
 * without duplicates, this function returns \f$| A \setminus B |\f$.
 */
template <typename InputIterator1, typename InputIterator2>
std::size_t
set_difference_count(InputIterator1 a_iter, InputIterator1 a_end, InputIterator2 b_iter, InputIterator2 b_end)
{
    std::size_t a_count = 0;
    while (a_iter != a_end) {
        if (b_iter == b_end) {
            ++a_iter;
            ++a_count;
        } else if (*a_iter == *b_iter) {
            ++a_iter;
            ++b_iter;
        } else if (*a_iter < *b_iter) {
            ++a_iter;
            ++a_count;
        } else {
            ++b_iter;
        }  // *a_iter > *minus_iter
    }

    return a_count;
}

}  // end namespace hqspre

#endif
