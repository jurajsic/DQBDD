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

#ifndef HQSPRE_FORMULA_BLOCKED_ET_AL_IPP_
#define HQSPRE_FORMULA_BLOCKED_ET_AL_IPP_

/**
 * \file formula_blocked_et_al.ipp
 * \brief Implementation of template functions for blocked clause elimination
 * \author Ralf Wimmer
 * \date 02/2016
 */

namespace hqspre {

/**
 * \brief Check if the given clause is blocked.
 * \param current_clause the clause that is to be checked.
 * \return the blocking lit if the clause is blocked, otherwise 0
 * \pre The clause needs to be sorted and may not contain duplicate literals.
 * \pre Assumes that "_seen" vector was initialized with current_clause
 * \note The clause does not need to be contained in the formula.
 */
template <typename Container>
inline Literal
Formula::clauseBlocked(const Container& current_clause) const
{
    val_assert(_unit_stack.empty());

    const auto blocking_lit = std::find_if(current_clause.cbegin(), current_clause.cend(), [this](const Literal lit) {
        return _prefix->isExistential(lit2var(lit)) && clauseBlockedByLit(lit);
    });

    if (blocking_lit != current_clause.cend()) {
        return *blocking_lit;
    }

    return 0;
}

/**
 * \brief Check if the given clause is blocked by given lit.
 * \param blocking_lit the literal for which we check which clauses it blocks.
 * \return true if any clause is blocked by `lit`
 * \pre Assumes that "_seen" vector was initialized with corresponding clause
 * \note The clause does not need to be contained in the formula.
 * \note The clause does not to be sorted
 */
inline bool
Formula::clauseBlockedByLit(const Literal blocking_lit) const
{
    const Variable blocking_var = lit2var(blocking_lit);
    const Literal  neg_lit      = negate(blocking_lit);

    for (const ClauseID other_c_nr : _occ_list[neg_lit]) {
        val_assert(!clauseDeleted(other_c_nr));

        const Clause& other_clause = _clauses[other_c_nr];

        _process_limit.decreaseLimitBy(2, static_cast<int>(other_clause.size()));

        bool tautology = false;
        // Now check for every other literal in other_clause whether we gain a
        // tautology
        for (const Literal lit : other_clause) {
            const Variable current_var = lit2var(lit);

            // "clause" contains a negated lit -> resolvent is potential tautology
            if (_seen[negate(lit)] && current_var != blocking_var) {
                // check for correct dependencies
                if (_prefix->dependenciesSubset(current_var, blocking_var)) {
                    tautology = true;
                    break;
                }
            }
        }

        // If at least one resolvent does not become a tautology, "clause" is not
        // blocked
        if (!tautology) {
            return false;
        }
    }  // end for (other_clause)
    return true;
}

}  // end namespace hqspre

#endif
