/*
 * This file is part of HQSpre.
 *
 * Copyright 2016/18 Ralf Wimmer, Sven Reimer, Paolo Marin, Bernd Becker
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
#include <vector>

#include <easylogging++.hpp>

#include "formula.hpp"
#include "literal.hpp"
#include "propagator.hpp"

/**
 * \file formula_vivification.cpp
 * \brief Implementation of clause vivification
 * \author Ralf Wimmer
 * \date 11/2017
 */

namespace hqspre {

/**
 * \brief Applies vivification to the clauses of the formula.
 *
 * For more details on vivification see:
 * C. Piette, Y. Hamadi, L. Sais: <i>Vivifying Propositional Clausal
 * Formulae.</i><br> Proc. of ECAI 2008, IOS Press.
 * DOI: 10.3233/978-1-58603-891-5-525
 *
 * \return true if the formula has been modified.
 */
bool
Formula::applyVivification()
{
    if (!_unit_stack.empty()) {
        unitPropagation();
    }

    ScopeTimer v_timer(getTimer(WhichTimer::VIVIFICATION));

    VLOG(1) << __FUNCTION__;

    _process_limit.setLimit(PreproMethod::VIVIFICATION);
    SatPropagator prop(_prefix, _clauses, _implications);

    const auto old_viv_lits    = stat(Statistics::VIVIFIED_LITERALS);
    const auto old_viv_clauses = stat(Statistics::VIVIFIED_CLAUSES);

    bool result = false;
    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (_process_limit.reachedLimit()) {
            VLOG(1) << "Vivification terminated due to process limit.";
            break;
        }
        prop.clear();
        if (_interrupt) {
            break;
        }
        if (clauseDeleted(c_nr)) {
            continue;
        }
        if (_clauses[c_nr].size() < _settings.vivify_min_size) {
            continue;
        }
        if (vivifyClause(c_nr, prop)) {
            result = true;
        }
    }

    VLOG(2) << "Vivification removed " << (stat(Statistics::VIVIFIED_CLAUSES) - old_viv_clauses) << " clauses and "
            << (stat(Statistics::VIVIFIED_LITERALS) - old_viv_lits) << " literals.";

    return result;
}

/**
 * \brief Applies vivification to an individual clause.
 *
 * Parameter _settings.vivify_delete determines whether vivification is used to
 * delete clauses or to shorten clauses.
 *
 * \param c_nr the ID of the clause to be vivified
 * \param prop a Propagator object for performing Boolean constraint propagation
 * \returns true if the clause was modified or deleted
 * \sa Formula::applyVivification()
 */
bool
Formula::vivifyClause(const ClauseID c_nr, SatPropagator& prop)
{
    val_assert(c_nr < _clauses.size());

    if (clauseDeleted(c_nr)) {
        return false;
    }
    if (!_unit_stack.empty()) {
        unitPropagation();
        prop.reset();
    }

    const std::vector<TruthValue>& assignment  = prop.getAssignment();
    const Clause&                  clause      = _clauses[c_nr];
    const std::size_t              clause_size = clause.size();

    // Consider the literals in a randomized order
    std::vector<Literal> buffer(clause.cbegin(), clause.cend());
    val_assert(buffer.size() == clause.size());

    std::random_shuffle(buffer.begin(), buffer.end());
    val_assert(buffer.size() == clause.size());

    std::vector<Literal> lits;
    lits.reserve(clause.size());

    for (int pos = 0; pos < static_cast<int>(clause_size) - 1; ++pos) {
        if (_interrupt) {
            return false;
        }

        const Literal current_lit = buffer[pos];
        lits.push_back(negate(current_lit));
        _process_limit.decreaseLimitBy(_occ_list[negate(current_lit)].size());

        try {
            prop.bcp(negate(current_lit), static_cast<int>(c_nr));
        } catch (const ConflictException&) {
            std::transform(lits.cbegin(), lits.cend(), lits.begin(), negate);

            if (!_settings.vivify_delete) {
                VLOG(3) << "Vivifying #" << c_nr << " led to a conflict. Removing " << (clause.size() - lits.size())
                        << " literals.";
                stat(Statistics::VIVIFIED_LITERALS) += (clause_size - lits.size());

                prop.removeClause(c_nr);
                removeClause(c_nr);

                const auto nr = addClause(std::move(lits), true, false);
                if (nr >= 0) {
                    prop.addClause(nr);
                } else {
                    unitPropagation();
                    prop.reset();
                }
            } else {
                VLOG(3) << "Vivifying #" << c_nr << " led to a conflict. Removing clause.";
                prop.removeClause(c_nr);
                removeClause(c_nr);
                stat(Statistics::VIVIFIED_CLAUSES) += 1;
            }

            return true;
        }

        // Check if another literal of the clause (or its negation) is implied
        for (int other_pos = pos + 1; other_pos < static_cast<int>(clause_size); ++other_pos) {
            if (pos == static_cast<int>(clause.size()) - 2 && other_pos == pos + 1) {
                break;
            }

            const Literal  other_lit = buffer[other_pos];
            const Variable other_var = lit2var(other_lit);
            if (isSatisfied(other_lit, assignment[other_var])) {
                // We can shorten or delete the clause
                if (_settings.vivify_delete) {
                    stat(Statistics::VIVIFIED_CLAUSES) += 1;
                    VLOG(3) << "Vivifying #" << c_nr << " led to an implied satisfied literal. Removing clause.";

                    prop.removeClause(c_nr);
                    removeClause(c_nr);
                } else {
                    std::transform(lits.cbegin(), lits.cend(), lits.begin(), negate);
                    lits.push_back(other_lit);

                    stat(Statistics::VIVIFIED_LITERALS) += (clause_size - lits.size());
                    VLOG(3) << "Vivifying #" << c_nr << " led to an implied satisfied literal. Removing "
                            << (clause.size() - lits.size()) << " literals.";

                    prop.removeClause(c_nr);
                    removeClause(c_nr);
                    const auto nr = addClause(std::move(lits), true, false);
                    if (nr >= 0) {
                        prop.addClause(nr);
                    } else {
                        unitPropagation();
                        prop.reset();
                    }
                }
                return true;
            } else if (isUnsatisfied(other_lit, assignment[other_var])) {
                // We can remove at least 'other_lit' from the clause

                // We abuse 'lits' to store those literals that can be deleted from the
                // clause
                lits.clear();
                lits.push_back(other_lit);
                for (std::size_t next_pos = other_pos + 1; next_pos < clause_size; ++next_pos) {
                    const Literal next_lit = buffer[next_pos];
                    if (isUnsatisfied(next_lit, assignment[lit2var(next_lit)])) {
                        lits.push_back(next_lit);
                    }
                }

                VLOG(3) << "Vivifying #" << c_nr << " led to implied unsatisfied literals. Removing " << lits.size()
                        << " literals.";
                prop.removeClause(c_nr);
                for (const auto l : lits) {
                    removeLiteral(c_nr, l);
                }
                prop.addClause(c_nr);

                stat(Statistics::VIVIFIED_LITERALS) += lits.size();
                return true;
            }
        }
    }
    return false;
}

}  // end namespace hqspre
