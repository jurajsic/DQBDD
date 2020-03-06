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

#include "clause.hpp"

#include <algorithm>
#include <iterator>
#include <ostream>
#include <set>
#include <vector>

#include <easylogging++.hpp>
#include "aux.hpp"
#include "literal.hpp"

/**
 * \file clause.cpp
 * \brief Implementation of operations on clauses
 * \author Ralf Wimmer
 * \date 05-06/2016
 */

namespace hqspre {

/**
 * \brief Creates a clause from a vector of literals
 *
 * \param literals the vector of literals
 * \param needs_check if true the literals get sorted and duplicate literals are
 * removed \param status status of the clause (mandatory, optional, deleted)
 * \pre If needs_check is false, the user has to guarantee that the clause is
 *      sorted and does not contain duplicate literals.
 */
Clause::Clause(ClauseData&& literals, const bool needs_check, const ClauseStatus status) noexcept :
    _literals(std::move(literals)),
    _status(status)
{
    if (needs_check) {
        // Sort the literals and remove duplicates
        std::sort(_literals.begin(), _literals.end());
        _literals.erase(std::unique(_literals.begin(), _literals.end()), _literals.end());
    }
    val_assert(std::is_sorted(_literals.cbegin(), _literals.cend()));
    val_assert(is_unique(_literals.cbegin(), _literals.cend()));

    computeSignature();
}

/**
 * \brief Checks if the clause is a tautology.
 *
 * \return true if the clause contains both a variable and its negation at the
 * same time.
 */
bool
Clause::isTautology() const noexcept
{
    if (_literals.empty()) return false;

    for (std::size_t pos = 0; pos < _literals.size() - 1; ++pos) {
        val_assert(_literals[pos] != _literals[pos + 1]);
        if (lit2var(_literals[pos]) == lit2var(_literals[pos + 1])) {
            return true;
        }
    }
    return false;
}

/**
 * \brief Checks if a clause contains a given literal.
 * \return true if the literal appears in the clause, false otherwise.
 */
bool
Clause::containsLiteral(const Literal lit) const
{
    return std::binary_search(_literals.cbegin(), _literals.cend(), lit);
}

/**
 * \brief Computes a signature of the clause.
 *
 * This signature can be used to make subsumption checking more
 * efficient. For clauses \f$C\f$ and \f$C'\f$, if
 * \f$\mathrm{sig}(C)\land\neg\mathrm{sig}(C') \neq 0\f$, then
 * \f$C\f$ cannot subsume \f$C'\f$, i.e., \f$C\not\subseteq C'\f$.
 *
 * \note If a clause is modified, the signature has to be recomputed.
 * \sa Formula::subsumption()
 * \sa Clause::subsetOf(const Clause&)
 */
void
Clause::computeSignature()
{
    _signature = 0u;

    for (const Literal lit : _literals) {
        addSignatureLit(_signature, lit);
        addSignatureLit(_var_signature, lit2var(lit));
    }
}

/**
 * \brief Checks if the current clause is a subset of `other`.
 *
 * Assumes that all involved clauses are sorted
 * \sa Clause::computeSignature()
 * \param other the clause that is checked if it contains the current clause
 * \param other_signature the signature of the clause `other`
 */
template<typename Container>
bool
Clause::subsetOf(const Container& other, const std::uint64_t other_signature) const
{
    if (this->size() > other.size()) return false;
    if ((this->getSignature() & ~(other_signature)) != 0) return false;

    return std::includes(this->cbegin(), this->cend(), other.cbegin(), other.cend());
}

// Excplicit template initialization
template bool Clause::subsetOf(const std::set<Literal>& other, const std::uint64_t other_signature) const;
template bool Clause::subsetOf(const std::vector<Literal>& other, const std::uint64_t other_signature) const;
template bool Clause::subsetOf(const Clause& other, const std::uint64_t other_signature) const;

/**
 * \brief Computes the resolvent of the two clauses w.r.t. the given variable.
 *
 * \param c_pos clause which contains 'var' positive
 * \param c_neg clause which contains 'var' negative
 * \param var pivot variable
 * \return the resolvent of the two clauses.
 * \note The two clauses must be sorted, non-tautological and may not
 *       contain duplicate literals. The resolvent then has the same properties.
 */
Clause
resolve(const Clause& c_pos, const Clause& c_neg, const Variable var)
{
    val_assert(c_pos.size() >= 1 && c_neg.size() >= 1);
    val_assert(c_pos.containsLiteral(var2lit(var, false)));
    val_assert(c_neg.containsLiteral(var2lit(var, true)));

    Clause::ClauseData result;
    result.reserve(c_pos.size() + c_neg.size() - 2);

    const Literal lit_pos = var2lit(var, false);
    const Literal lit_neg = var2lit(var, true);

    auto iter_pos = c_pos.cbegin();
    auto iter_neg = c_neg.cbegin();
    while (true) {
        // Have we reached the end of the positive clause?
        if (iter_pos == c_pos.cend()) {
            while (iter_neg != c_neg.cend()) {
                if (*iter_neg != lit_neg) result.push_back(*iter_neg);
                ++iter_neg;
            }
            break;
        }

        // Have we reached the end of the negative clause?
        if (iter_neg == c_neg.cend()) {
            while (iter_pos != c_pos.cend()) {
                if (*iter_pos != lit_pos) result.push_back(*iter_pos);
                ++iter_pos;
            }
            break;
        }

        // We haven't reached the end of any clause
        if (*iter_pos == *iter_neg) {
            result.push_back(*iter_pos);
            ++iter_pos;
            ++iter_neg;
        } else if (*iter_pos < *iter_neg) {
            if (*iter_pos != lit_pos) result.push_back(*iter_pos);
            ++iter_pos;
        } else if (*iter_neg < *iter_pos) {
            if (*iter_neg != lit_neg) result.push_back(*iter_neg);
            ++iter_neg;
        }
    }

    return Clause(std::move(result), false);
}

/**
 * \brief Prints a clause in DIMACS format.
 */
std::ostream&
operator<<(std::ostream& stream, const Clause& clause)
{
    for (const Literal lit : clause) {
        stream << lit2dimacs(lit) << ' ';
    }
    stream << '0';

    return stream;
}

/**
 * \brief Checks if the clause is in a consistent state.
 *
 * Consistent means that it is sorted in ascending order,
 * does not contain duplicate literals and is not a tautology.
 * \return true iff the clause is consistent.
 */
bool
Clause::checkConsistency() const
{
    for (const Literal lit : _literals) {
        if (lit < 2) {
            LOG(ERROR) << "Invalid literal in clause: " << lit;
            return false;
        }
    }

    if (!std::is_sorted(_literals.cbegin(), _literals.cend())) {
        LOG(ERROR) << "Clause is not sorted: " << *this;
        return false;
    }

    if (!is_unique(_literals.cbegin(), _literals.cend())) {
        LOG(ERROR) << "Clause contains duplicate literals: " << *this;
        return false;
    }

    if (isTautology()) {
        LOG(ERROR) << "Clause is a tautology: " << *this;
        return false;
    }

    return true;
}

}  // end namespace hqspre
