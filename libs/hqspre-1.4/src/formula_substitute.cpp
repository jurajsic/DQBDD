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

#include <cstdint>
#include <iterator>
#ifdef SKOLEM
#include <memory>
#endif
#include <utility>
#include <vector>

#include <easylogging++.hpp>
#include "aux.hpp"
#include "clause.hpp"
#include "formula.hpp"
#include "gate.hpp"
#include "literal.hpp"
#include "prefix.hpp"
#include "timer.hpp"

/**
 * \file formula_substitute.cpp
 * \brief Implementation of gate substitution and gate rewriting
 * \author Ralf Wimmer
 * \author Paolo Marin
 * \date 06/2016
 */

namespace hqspre {

/**
 * \brief Replaces a Tseitin variable by its gate definition.
 *
 * This function can be used to eliminate existential variables
 * which result from Tseitin transformation. It typically creates
 * fewer clauses than variable elimination by resolution.<br>
 * See: Een, Biere: <i>Effective Preprocessing in SAT through
 *      Variable and Clause Elimination</i>, SAT 2005.
 * \param g the gate defining the substituted variable
 */
bool Formula::substituteGate(Gate& g)
{
    val_assert(_unit_stack.empty());
    val_assert(_gates.gateValid(g));

    VLOG(3) << __FUNCTION__ << " substituting gate " << g;


    const auto approx_costs = computeSubstitutionCosts(g);
    if (approx_costs > _settings.max_substitution_cost + 80) {
        VLOG(3) << " - aborted due to estimated cost " << approx_costs << '.';
        return false;
    }

    VLOG(3) << " - estimated cost " << approx_costs;

    const Literal o_lit = g._output_literal;
    const Literal no_lit = negate(o_lit);
    std::vector<Clause> to_add;

    signed int actual_costs = 0;

    for (const ClauseID c_nr: g._encoding_clauses) {
        const auto& clause = _clauses[c_nr];
        if (clause.containsLiteral(o_lit)) {
            for (const ClauseID other_c_nr: _occ_list[no_lit]) {
                if (isNegative(o_lit)) {
                    to_add.push_back(resolve(_clauses[other_c_nr], clause, lit2var(o_lit)));
                } else {
                    to_add.push_back(resolve(clause, _clauses[other_c_nr], lit2var(o_lit)));
                }
                auto& lastclause = to_add.back();
                if (lastclause.isTautology()) {
                    to_add.pop_back();
                } else {
                    actual_costs += static_cast<int>(lastclause.size());
                }
            }
        } else if (clause.containsLiteral(no_lit)) {
            for (const ClauseID other_c_nr: _occ_list[o_lit]) {
                if (isPositive(o_lit)) {
                    to_add.push_back(resolve(_clauses[other_c_nr], clause, lit2var(o_lit)));
                } else {
                    to_add.push_back(resolve(clause, _clauses[other_c_nr], lit2var(o_lit)));
                }
                auto& lastclause = to_add.back();

                if (lastclause.isTautology()) {
                    to_add.pop_back();
                } else {
                    actual_costs += static_cast<int>(lastclause.size());
                }
            }
        } else {
            LOG(ERROR) << "Clause " << clause << " does not contain output literal " << lit2dimacs(o_lit);
            LOG(ERROR) << "  The clause is: " << clause << ", the substituted variable: " << lit2var(o_lit);
            LOG(ERROR) << "  The gate is: " << g;
            LOG(ERROR) << "  Its encoding clauses are:";
            for (const auto c_nr: g._encoding_clauses) LOG(ERROR) << "   " << clause[c_nr];
            val_assert(false);
        }

    }

    for (const ClauseID c_nr: _occ_list[o_lit]) {
        actual_costs -= static_cast<int>(_clauses[c_nr].size());
    }
    for (const ClauseID c_nr: _occ_list[no_lit]) {
        actual_costs -= static_cast<int>(_clauses[c_nr].size());
    }

    if (actual_costs > _settings.max_substitution_cost) {
        VLOG(3) << " - aborted due to actual costs of " << actual_costs << '.';
        return false;
    }


    // Delete the clauses that are no longer needed
    while (!_occ_list[o_lit].empty()) {
        removeClause(_occ_list[o_lit].front());
    }
    while (!_occ_list[no_lit].empty()) {
        removeClause(_occ_list[no_lit].front());
    }


    for (const Clause& c: to_add) {
        setSeen(c);
        if (!isForwardSubsumed(c, -1)) {
            clearSeen(c);
            addClause(std::move(c));
        } else {
            clearSeen(c);
            actual_costs -= static_cast<int>(c.size());
        }
    }

#ifdef SKOLEM
    if (_settings.skolem) {
        _skolem_data.push_back(std::make_unique<SkolemGate>(g));
    }
#endif

    removeVar(lit2var(o_lit));
    ++stat(Statistics::SUBSTITUTION);

    VLOG(3) << " - succeeded with actual costs " << actual_costs << '.';

    return true;
}

/**
 * \brief Replaces a Tseitin variable by a double-Plaisted encoding
 *
 * This function can be used to replace existential variables
 * which result from Tseitin transformation by splitting a gate's output
 * with two monotone outputs.<br>
 * See: Giunchiglia, Marin, Narizzano: <i>sQueezeBF: An Effective Preprocessor
 *      for QBFs based on Equivalence Reasoning</i>, SAT 2010.
 * \param g the gate defining the substituted variable; g must be a valid gate.
 */
void Formula::rewriteGate(Gate& g)
{
    val_assert(_gates.gateValid(g));

    if (g._type != GateType::AND_GATE && g._type != GateType::XOR_GATE && g._type != GateType::MUX_GATE) return;
    if (varDeleted(lit2var(g._output_literal))) return;

    ScopeTimer rew(getTimer(WhichTimer::REWRITING));
    VLOG(3) << __FUNCTION__ << " Rewriting gate " << g << '.';

    if (!_unit_stack.empty()) unitPropagation();

    Clause::ClauseData new_clause;
    std::vector<Clause::ClauseData> to_add;
    to_add.reserve(g._encoding_clauses.size() + 1);

    const Literal y = g._output_literal;
    const Literal y_neg = negate(y);

    val_assert(!varDeleted(lit2var(y)));

    // new variable representing the second output of the gate in the
    // modified Plaisted encoding
    const Literal z = var2lit(copyVar(lit2var(y)), false);
    const Literal z_neg = negate(z);

    VLOG(3) << __FUNCTION__ << " Variable " << lit2var(y) << " doubled as " << lit2var(z);

    for (const ClauseID c_nr: g._encoding_clauses) {
      removeClause(c_nr);
    }

    // Delete all optional clauses involving y or neg_y
    for (const ClauseID c_nr: _occ_list[y]) {
        if (clauseOptional(c_nr)) removeClause(c_nr);
    }
    for (const ClauseID c_nr: _occ_list[y_neg]) {
        if (clauseOptional(c_nr)) removeClause(c_nr);
    }

    // add efficiency clause
    new_clause.clear();
    new_clause.push_back(y_neg);
    new_clause.push_back(z_neg);
    to_add.push_back(std::move(new_clause));

    // Replace all occurrences of neg_y by z
    replaceLiteralMono(y_neg, z);

    if (g._type == GateType::AND_GATE) {
        // rebuild gate encoding

        // long clause: replace y by neg_z
        new_clause.clear();
        new_clause.push_back(z_neg);
        for (const Literal lit: g._input_literals) {
            // Skip literal if variable was deleted before
            val_assert(!varDeleted(lit2var(lit)));
            new_clause.push_back(negate(lit));
        }
        to_add.push_back(std::move(new_clause));

        // short clauses
        // keep unmodified
        for (const Literal lit: g._input_literals) {
            new_clause.clear();
            new_clause.push_back(y_neg);
            new_clause.push_back(lit);
            // Skip clause if variable was deleted before
            to_add.push_back(std::move(new_clause));
        }

    } else if (g._type == GateType::XOR_GATE) {

        // Rebuild the gate encoding
        val_assert(g._input_literals.size() == 2);
        const Literal a = g._input_literals[0];
        const Literal b = g._input_literals[1];

        new_clause.clear();
        new_clause.push_back(y_neg);
        new_clause.push_back(a);
        new_clause.push_back(b);
        to_add.push_back(std::move(new_clause));

        new_clause.clear();
        new_clause.push_back(y_neg);
        new_clause.push_back(negate(a));
        new_clause.push_back(negate(b));
        to_add.push_back(std::move(new_clause));

        new_clause.clear();
        new_clause.push_back(z_neg);
        new_clause.push_back(negate(a));
        new_clause.push_back(b);
        to_add.push_back(std::move(new_clause));

        new_clause.clear();
        new_clause.push_back(z_neg);
        new_clause.push_back(a);
        new_clause.push_back(negate(b));
        to_add.push_back(std::move(new_clause));

    } else if (g._type == GateType::MUX_GATE) {
        val_assert(g._input_literals.size() == 3);
        const Literal sel = g._input_literals[0];
        const Literal x1  = g._input_literals[1];
        const Literal x0  = g._input_literals[2];

        new_clause.clear();
        new_clause.push_back(z_neg);
        new_clause.push_back(negate(sel));
        new_clause.push_back(negate(x1));
        to_add.push_back(std::move(new_clause));

        new_clause.clear();
        new_clause.push_back(y_neg);
        new_clause.push_back(negate(sel));
        new_clause.push_back(x1);
        to_add.push_back(std::move(new_clause));

        new_clause.clear();
        new_clause.push_back(y_neg);
        new_clause.push_back(sel);
        new_clause.push_back(x0);
        to_add.push_back(std::move(new_clause));

        new_clause.clear();
        new_clause.push_back(z_neg);
        new_clause.push_back(sel);
        new_clause.push_back(negate(x0));
        to_add.push_back(std::move(new_clause));

    } else {
        val_assert_msg(false, "Gate type is not supported by rewriting!");
        return;
    }

    for (std::size_t i = 0; i < to_add.size(); ++i) {
        addClause(std::move(to_add[i]));
    }

    ++stat(Statistics::REWRITING);
}


/**
 * \brief Computes the cost of replacing a Tseitin variable by its gate definition.
 *
 * The costs are given as the number of literals by which the formula size increases.
 * Negative costs means that the formula becomes shorter. The value is only an
 * over-approximation as it does not take into account that substitution may lead
 * to tautological clauses or duplicate literals.
 *
 * \todo Better estimation of the substitution costs by taking into
 *       account that resolvents of gate-defining clauses are often
 *       tautologies.
 */
int Formula::computeSubstitutionCosts(const Gate& g) const
{
    int result = 0;
    int pos_cost = 0;
    int neg_cost = 0;

    for (const ClauseID c_nr: _occ_list[g._output_literal]) {
        if (!clauseOptional(c_nr)) pos_cost += static_cast<int>(_clauses[c_nr].size() - 1);
    }

    for (const ClauseID c_nr: _occ_list[negate(g._output_literal)]) {
        if (!clauseOptional(c_nr)) neg_cost += static_cast<int>(_clauses[c_nr].size() - 1);
    }

    for (const ClauseID c_nr: g._encoding_clauses) {
        if (_clauses[c_nr].containsLiteral(g._output_literal)) result += static_cast<int>(_clauses[c_nr].size() - 1) * neg_cost;
        else result += static_cast<int>(_clauses[c_nr].size() - 1) * pos_cost;
    }
    result -= (pos_cost + neg_cost + static_cast<int>(_occ_list[g._output_literal].size()) + static_cast<int>(_occ_list[negate(g._output_literal)].size()));

    return result;
}



/**
 * \brief Eliminates Tseitin variables by substitution.
 *
 * Substitution first determines the encoded gates. These gates
 * are handled in reverse topological order, i.e., from the circuit
 * outputs to the circuit inputs. This avoids that substitution
 * re-introduces variables which have already been substituted.
 * The substitution of a gates is only performed if it increases
 * the formula size by at most some Settings::max_substitution_cost literals.
 * \sa Formula::substitutionCosts(const Gate&)
 * \return true if the formula was modified
 */
bool Formula::applySubstitution()
{
    VLOG(1) << __FUNCTION__;

    ScopeTimer gate_sub(getTimer(WhichTimer::SUBSTITUTION));

    const std::size_t old_stat_substitute = stat(Statistics::SUBSTITUTION);
    const std::size_t old_stat_rewriting = stat(Statistics::REWRITING);

    determineGates(true, true, true, _settings.semantic_gates);

    val_assert(checkConsistency());

    for (auto it = _gates.rbegin(); it != _gates.rend(); ++it)
    {
        if (_interrupt) break;
        if (!_unit_stack.empty()) fastPreprocess();

        val_assert(_gates.checkConsistency());

        // Check if any of the gate's clauses has been
        // deleted (by subsumption checks during substitution
        // of other clauses)
        if (!_gates.gateValid(*it)) continue;

        if (_settings.preserve_gates) {
            if (!_gates.isGateInput(lit2var(it->_output_literal)) && substituteGate(*it)) {
                _gates.invalidateGate(*it);
            }
        } else {
            if (substituteGate(*it)) {
                _gates.invalidateGate(*it);
            } else if (_settings.rewrite && _prefix->type() != PrefixType::DQBF) {
                getTimer(WhichTimer::SUBSTITUTION).stop();
                rewriteGate(*it);
                _gates.invalidateGate(*it);
                getTimer(WhichTimer::SUBSTITUTION).start();
            }
        }
    }

    VLOG(2) << __FUNCTION__ << " substituted " << (stat(Statistics::SUBSTITUTION) - old_stat_substitute)
            << " and rewrote " << (stat(Statistics::REWRITING) - old_stat_rewriting) << " gates.";

    return (stat(Statistics::SUBSTITUTION) > old_stat_substitute) || (stat(Statistics::REWRITING) > old_stat_rewriting);
}


} // end namespace hqspre
