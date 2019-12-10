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

 * HQSpre is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with HQSpre. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HQSPRE_LITERAL_HPP_
#define HQSPRE_LITERAL_HPP_

/**
 * \file literal.hpp
 * \author Ralf Wimmer
 * \date 05/2016
 * \brief Manipulation of literals
*/

#include "aux.hpp"

namespace hqspre {


using Literal  = std::uint32_t;
using Variable = std::uint32_t;
using ClauseID = std::size_t;


/**
 * \brief Converts a literal to its underlying variable.
 */
inline constexpr Variable lit2var(const Literal lit) noexcept
{
    return (lit >> 1);
}


/**
 * \brief Converts a variable to a literal with the given polarity.
 */
inline constexpr Literal var2lit(const Variable var, const bool negated = false) noexcept
{
    return ((var << 1) | (negated ? 1 : 0));
}


/**
 * \brief Converts a literal to DIMACS format.
 *
 * Positive literals correspond to the variable index,
 * negative literals to the negative variable index.
 */
inline constexpr int lit2dimacs(const Literal lit) noexcept
{
    return ((lit & 1) == 1) ? (-1 * static_cast<int>(lit2var(lit))) : (static_cast<int>(lit2var(lit)));
}

/**
 * \brief Negates a given literal.
 */
inline constexpr Literal negate(const Literal lit) noexcept
{
    return (lit ^ 1);
}


/**
 * \brief Negates a given literal iff 'negate' is true.
 */
inline constexpr Literal negate_if(const Literal lit, const bool negate) noexcept
{
    return (negate ? (lit ^ 1) : lit);
}


/**
 * \brief Checks whether a given literal is positive.
 */
inline constexpr bool isPositive(const Literal lit) noexcept
{
    return ((lit & 1) == 0);
}


/**
 * \brief Checks whether a given literal is negative.
 */
inline constexpr bool isNegative(const Literal lit) noexcept
{
    return ((lit & 1) != 0);
}


/**
 * \brief Returns true iff the given value satisfied the literal.
 */
inline constexpr bool isSatisfied(const Literal lit, const TruthValue val) noexcept
{
    return (val == TruthValue::TRUE && isPositive(lit)) || (val == TruthValue::FALSE && isNegative(lit));
}

/**
 * \brief Returns true iff the given value violates the literal.
 */
inline constexpr bool isUnsatisfied(const Literal lit, const TruthValue val) noexcept
{
    return (val == TruthValue::TRUE && isNegative(lit)) || (val == TruthValue::FALSE && isPositive(lit));
}

/**
 * \brief Returns true iff the given literal is unassigned.
 */
inline constexpr bool isUnassigned(const Literal /* lit */, const TruthValue val) noexcept
{
    return (val == TruthValue::UNKNOWN);
}



inline constexpr TruthValue makeSatisfied(const Literal lit)
{
    if (isNegative(lit)) return TruthValue::FALSE;
    else return TruthValue::TRUE;
}


inline constexpr TruthValue makeUnsatisfied(const Literal lit)
{
    if (isNegative(lit)) return TruthValue::TRUE;
    else return TruthValue::FALSE;
}


inline constexpr TruthValue makeUnassigned(const Literal /* lit */)
{
    return TruthValue::UNKNOWN;
}


} // end namespace hqspre

#endif
