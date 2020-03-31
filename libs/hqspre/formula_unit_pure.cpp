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

#include <vector>

#include <easylogging++.hpp>

#include "auxil.hpp"
#include "clause.hpp"
#include "exceptions.hpp"
#include "formula.hpp"
#include "gate.hpp"
#include "literal.hpp"

/**
 * \file formula_unit_pure.cpp
 * \brief Implementation of unit propagation and pure literal detection
 * \author Ralf Wimmer
 * \date 02/2016
 */

namespace hqspre {

/**
 * \brief Puts a literal onto the unit stack.
 *
 * If the same literal has already been with assigned opposite sign, an
 * UNSATException is thrown.
 *
 * \param lit the unit literal
 * \param is_pure If false, a universal unit literal leads to an
 *        UNSATException. However, for universal pure literals this is not
 *        correct. In this case, the parameter must be set to true.
 * \throw UNSATException if a universal unit literal is added and
 * allow_universal is false. \throw UNSATException if the same unit literal just
 * with opposite sign has beed added before.
 */
void
Formula::pushUnit(const Literal lit, const PureStatus status)
{
    VLOG(3) << "Literal " << lit2dimacs(lit) << " is " << (status == PureStatus::UNIT ? "unit" : "pure");

    const Variable var = lit2var(lit);
    if (status == PureStatus::UNIT && isUniversal(var)) {
        throw UNSATException("Universal literal is unit.");
    }

    if (isUnsatisfied(lit, _assignment[var])) {
        throw UNSATException("Contradicting unit literals.");
    }

    _unit_stack.push_back(lit);
    _assignment[var] = makeSatisfied(lit);

    if (status == PureStatus::PURE)
        ++stat(Statistics::PURE);
    else
        ++stat(Statistics::UNIT);
}

/**
 * \brief Processes unit clauses until none is left.
 *
 * Satisfied clauses are deleted, the negated unit literal is
 * removed from all clauses in which it appears, potentially leading
 * to new unit clauses. This is iterated until the formula does not
 * change anymore.
 * \return true iff the formula was modified.
 * \throws UNSATException if unitPropagation() determins a conflict.
 */
bool
Formula::unitPropagation()
{
    if (_unit_stack.empty()) return false;

    std::size_t count{0};

    // Process all variables on the 'unit_stack'.
    while (!_unit_stack.empty()) {
        const Literal top_lit = _unit_stack.back();
        _unit_stack.pop_back();
        const Variable top_var = lit2var(top_lit);
        if (varDeleted(top_var)) continue;
        const Literal neg_top_lit = negate(top_lit);

#ifdef SKOLEM
        if (_settings.skolem && isExistential(top_var)) {
            _skolem_data.push_back(std::make_unique<SkolemUnitPure>(top_lit));
        }
#endif

        // Delete all clauses in which 'top_lit' occurs
        while (!_occ_list[top_lit].empty()) {
            const ClauseID c_nr = _occ_list[top_lit].front();
            val_assert(!clauseDeleted(c_nr));
            // check for new pures after deletion of the clause
            checkPure(_clauses[c_nr], top_lit);
            removeClause(c_nr);
        }

        // Delete ~top_lit from all clauses where it occurs
        // If the clause is binary, it becomes unit -> call pushUnit()
        // If the clause is ternary, it becomes binary and we must
        //   update the implications.

        // Because 'removeClause(c_nr)' modifies _occ_list,
        // we may not simply iterate over it to delete the
        // satisfied clauses.
        std::vector<ClauseID> to_be_removed;

        for (const ClauseID c_nr : _occ_list[neg_top_lit]) {
            val_assert(!clauseDeleted(c_nr));
            auto& c = _clauses[c_nr];
            // binary clauses become unit
            if (c.size() == 2) {
                if (c[0] == neg_top_lit)
                    pushUnit(c[1], PureStatus::UNIT);
                else
                    pushUnit(c[0], PureStatus::UNIT);

                to_be_removed.push_back(c_nr);
            } else {
                removeLiteral(c_nr, neg_top_lit);
            }
        }
        for (const ClauseID c_nr : to_be_removed) {
            removeClause(c_nr);
        }

        // all lists of the propagated variable must be empty.
        val_assert(_occ_list[top_lit].empty());
        val_assert(_occ_list[neg_top_lit].empty());
        val_assert(_implications[top_lit].empty());
        val_assert(_implications[neg_top_lit].empty());

        // remove variable
        removeVar(top_var);
        ++count;
    }

    VLOG_IF(count > 0, 2) << __FUNCTION__ << " removed " << count << " unit and pure literals.";
    return count > 0;
}

/**
 * \brief Searches for pure literals and eliminates them.
 *
 * A literal is pure if it appears in the whole formula only in one polarity.
 * Existential pure literals can be set to true, universal pure literals to
 * false. This is done by calling pushUnit (settings allow_universal to true).
 * \throw UNSATException if propagating the pure literals leads to a conflict.
 * \return true if the formula was modified.
 */
bool
Formula::findPure()
{
    val_assert_msg(_unit_stack.empty(), "You must call Formula::unitPropagation() first!");

    const std::size_t old_stat_pure = stat(Statistics::PURE);
    bool              again         = true;

    do {
        again = false;

        for (Literal lit = minLitIndex(); lit <= maxLitIndex(); ++lit) {
            if (varDeleted(lit2var(lit))) continue;

            if (_occ_list[lit].empty()) {
                // 'lit' does not occur, so '~lit' is pure

                if (isExistential(lit2var(lit))) {
                    pushUnit(negate(lit), PureStatus::PURE);
                } else {
                    pushUnit(lit, PureStatus::PURE);
                }
                again = true;
                unitPropagation();
            } else if (_settings.impl_chains) {
                checkImplicationChain(lit);
            }
        }
    } while (again);

    return (stat(Statistics::PURE) > old_stat_pure);
}

/**
 * \brief Checks whether a literal is pure if one additional clause where
 * lit appears positively would be deleted
 *
 * This method is used in formula_bce for a quick check whether a literal
 * is pure after a blocked clause is eliminated
 *
 * \return true if literal would be pure
 */
bool
Formula::checkPure(Literal lit) const
{
    return ((_occ_list[lit].size() == 1) && !_occ_list[negate(lit)].empty());
}

/**
 * \brief Checks whether there will be new pure literals after "clause" is
 * deleted
 *
 * Checks every literal in clause except the "except_lit"
 */
void
Formula::checkPure(const Clause& clause, Literal except_lit)
{
    for (const Literal lit : clause) {
        if (lit != except_lit && _assignment[lit2var(lit)] == TruthValue::UNKNOWN && checkPure(lit)) {
            if (isExistential(lit2var(lit))) {
                pushUnit(negate(lit), PureStatus::PURE);
            } else {
                pushUnit(lit, PureStatus::PURE);
            }
        }
    }
}

}  // end namespace hqspre
