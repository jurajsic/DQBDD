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

#include <algorithm>
#include <set>
#include <vector>

#include <easylogging++.hpp>

#include "aux.hpp"
#include "clause.hpp"
#include "formula.hpp"
#include "literal.hpp"
#include "prefix.hpp"
#include "timer.hpp"

/**
 * \file formula_univ_red.cpp
 * \brief Implementation of universal reduction
 * \author Ralf Wimmer
 * \date 02/2016
 */

namespace hqspre {

/**
 * \brief Performs universal reduction on a single clause.
 *
 * Universal reduction removes a universal literal \f$u\f$ from
 * a clause \f$C\f$ if \f$C\f$ does not contain any existential literal
 * which depends on \f$u\f$.
 *
 * \note For the computation of Skolem functions, universal reduction can be
 *       ignored as it preserves Skolem functions.
 *
 * \param clause the clause that is to be reduced
 * \param c_nr the ID of the clause, if the clause is already in the database,
 * and -1 otherwise. \return true if the clause was modified \throws
 * UNSATException if universal reduction yields the empty clause (and `c_nr` is
 * not negative) \todo Use dependency schemes to make universal reduction more
 * powerful (see paper by Bayersdorff et al.)
 */
bool
Formula::universalReduction(Clause& clause, const int c_nr)
{
    val_assert(_prefix);
    ScopeTimer univ_red(getTimer(WhichTimer::UNIV_REDUCTION));

    _removed_lits.clear();

    if (_prefix->type() == PrefixType::QBF) {
        val_assert(_qbf_prefix && !_dqbf_prefix);

        // For QBF, determine the maximal existential level in the clause
        int max_level = -1;
        for (const Literal lit : clause) {
            if (isExistential(lit2var(lit))) {
                max_level = std::max(max_level, static_cast<int>(_qbf_prefix->getLevel(lit2var(lit))));
            }
        }

        // Perform the reduction: all universal literals with a higher level
        // than 'max_level' can be removed.
        std::size_t next_free = 0;
        for (std::size_t curr_index = 0; curr_index < clause.size(); ++curr_index) {
            const Variable var = lit2var(clause[curr_index]);
            if (isExistential(var) || static_cast<int>(_qbf_prefix->getLevel(var)) <= max_level) {
                // clause[curr_index] cannot be deleted
                if (next_free != curr_index) {
                    clause[next_free] = clause[curr_index];
                }
                ++next_free;
            } else {
                _removed_lits.push_back(clause[curr_index]);
            }
        }

        // If necessary, shrink the clause and update the occurrence/implication
        // lists
        if (!_removed_lits.empty()) {
            VLOG(3) << __FUNCTION__ << " removed " << _removed_lits.size() << " literal(s) from clause " << clause
                    << '.';

            if (c_nr >= 0) {
                _gates.touchClause(c_nr);
                if (clause.size() == 2) {
                    // remove implications
                    if (_removed_lits.size() == 2) {
                        removeImplication(negate(_removed_lits[0]), _removed_lits[1]);
                        removeImplication(negate(_removed_lits[1]), _removed_lits[0]);
                    } else {
                        removeImplication(negate(clause[0]), _removed_lits[0]);
                        removeImplication(negate(_removed_lits[0]), clause[0]);
                    }
                }
                // update occurrence lists
                for (const Literal lit : _removed_lits) {
                    removeFromOccList(lit, c_nr);
                }
                _clause_sizes[c_nr] = next_free;
            }

            clause.resize(next_free);
            clause.computeSignature();

            if (c_nr >= 0) {
                // Check for unit clauses and implications
                if (clause.empty()) {
                    throw UNSATException("Universal reduction created an empty clause.");
                } else if (clause.size() == 1) {
                    const Literal unit = clause[0];
                    removeClause(c_nr);
                    pushUnit(unit, PureStatus::UNIT);
                } else if (clause.size() == 2) {
                    addImplications(clause[0], clause[1], c_nr);
                }
            }
            stat(Statistics::UNIV_REDUCTION) += _removed_lits.size();
            return true;
        } else {
            return false;
        }

    } else {
        // formula is a DQBF
        val_assert(_prefix->type() == PrefixType::DQBF && _dqbf_prefix);

        const std::size_t num_univ = numUVars();

        // If all literals are existential, we can skip the clause
        if (std::all_of(clause.cbegin(), clause.cend(),
                        [this](const Literal lit) -> bool { return isExistential(lit2var(lit)); })) {
            return false;
        }

        if (std::any_of(clause.cbegin(), clause.cend(),
                        [this](const Literal lit) -> bool { return _dqbf_prefix->inRMB(lit2var(lit)); })) {
            return false;
        }

        // Compute the union of the dependency sets of all existential variables
        // in the clause.
        std::set<Variable> dependencies;
        for (const Literal lit : clause) {
            const Variable var = lit2var(lit);
            if (isExistential(var)) {
                fast_set_union(_dqbf_prefix->getDependencies(var), dependencies);
            }
            if (dependencies.size() == num_univ) {
                return false;
            }
        }

        // Perform the reduction
        std::size_t next_free = 0;
        for (std::size_t curr_index = 0; curr_index < clause.size(); ++curr_index) {
            const Variable var = lit2var(clause[curr_index]);
            if (isExistential(var) || dependencies.find(var) != dependencies.cend()) {
                // clause[curr_index] cannot be deleted
                if (next_free != curr_index) {
                    clause[next_free] = clause[curr_index];
                }
                ++next_free;
            } else {
                val_assert(isUniversal(lit2var(clause[curr_index])));
                _removed_lits.push_back(clause[curr_index]);
            }
        }

        if (next_free == 0) {
            throw UNSATException("Universal reduction created an empty clause.");
            return true;
        }

        // If necessary, shrink the clause
        if (!_removed_lits.empty()) {
            VLOG(3) << __FUNCTION__ << " removed " << _removed_lits.size() << " literals from clause " << clause << '.';

            if (c_nr >= 0) {
                _gates.touchClause(c_nr);
                if (clause.size() == 2) {
                    // remove implications
                    if (_removed_lits.size() == 2) {
                        removeImplication(negate(_removed_lits[0]), _removed_lits[1]);
                        removeImplication(negate(_removed_lits[1]), _removed_lits[0]);
                    } else {
                        removeImplication(negate(clause[0]), _removed_lits[0]);
                        removeImplication(negate(_removed_lits[0]), clause[0]);
                    }
                }
                // update occurrence lists
                for (const Literal lit : _removed_lits) {
                    removeFromOccList(lit, c_nr);
                }
                _clause_sizes[c_nr] = next_free;
            }
            clause.resize(next_free);
            clause.computeSignature();
            stat(Statistics::UNIV_REDUCTION) += _removed_lits.size();

            if (c_nr >= 0) {
                // Check for unit clauses and implications
                if (clause.empty()) {
                    throw UNSATException("Universal reduction created an empty clause.");
                    return true;
                } else if (clause.size() == 1) {
                    pushUnit(clause[0], PureStatus::UNIT);
                    removeClause(c_nr);
                } else if (clause.size() == 2) {
                    addImplications(clause[0], clause[1], c_nr);
                }
            }
            return true;
        } else {
            return false;
        }
    }
    return false;
}

/**
 * \brief Performs universal reduction on the whole formula.
 *
 * Universal reduction removes a universal literal \f$u\f$ from
 * a clause \f$C\f$ if \f$C\f$ does not contain any existential literal
 * which depends on \f$u\f$.
 *
 * \note For the computation of Skolem functions, universal reduction can be
 *       ignored as it preserves Skolem functions.
 *
 * \return true if the formula was modified
 * \throws UNSATException if universal reduction yields the empty clause
 */
bool
Formula::universalReduction()
{
    VLOG(1) << __FUNCTION__;

    bool modified = false;

    for (ClauseID c_nr = 0; c_nr <= maxClauseIndex(); ++c_nr) {
        if (clauseDeleted(c_nr)) {
            continue;
        }
        if (universalReduction(_clauses[c_nr], static_cast<int>(c_nr))) {
            modified = true;
        }
    }

    return modified;
}

}  // end namespace hqspre
