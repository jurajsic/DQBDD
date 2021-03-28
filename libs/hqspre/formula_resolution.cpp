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
#include <limits>
#include <unordered_set>
#include <utility>
#include <vector>

#include <easylogging++.hpp>

#include "auxil.hpp"
#include "clause.hpp"
#include "formula.hpp"
#include "gate.hpp"
#include "literal.hpp"
#include "prefix.hpp"
#include "process_limits.hpp"
#include "timer.hpp"
#include "varheap.hpp"

/**
 * \file formula_resolution.cpp
 * \brief Implementation of resolution-based operations on formulas
 * \author Ralf Wimmer
 * \author Sven Reimer
 * \date 2016-18
 */

namespace hqspre {

/**
 * \brief Checks if an existential variable can be eliminated by resolution.
 *
 * \return true iff the variable can be eliminated
 */
bool
Formula::isResolvable(const Variable var) const
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());

    return (isExistential(var) && _prefix->inRMB(var) && (!_settings.preserve_gates || !_gates.isGateOutput(var)));
}

/**
 * \brief Determines which variables may be eliminated using resolution.
 *
 * A variable \f$ x\f$ may be eliminated if
 * - it is a Tseitin variable (i.e., encodes a gate output)
 * - if the dependency sets of all literals the literal \f$ x\f$ appear with in
 * the same clause are a subset of \f$x\f$'s dependency set.
 * - the same as the second condition only for \f$\neg x\f$.
 * \return a vector of resolvable variables
 */
std::vector<Variable>
Formula::getResolvableVariables() const
{
    if (!_settings.preserve_gates) {
        const auto& rmb = _prefix->getRMB();
        if (!rmb.empty() && isExistential(*(rmb.begin()))) {
            return std::vector<Variable>(rmb.cbegin(), rmb.cend());
        } else {
            return std::vector<Variable>();
        }
    }

    std::vector<Variable> result;
    result.reserve(maxVarIndex() + 1 - numUVars());

    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (isResolvable(var)) result.push_back(var);
    }

    return result;
}

/**
 * \brief Computes the costs of eliminating an existential variable by
 * resolution. \param var the variable to be eliminated \return the number of
 * literals by which the formula length increases.
 */
int
Formula::computeResolutionCosts(const Variable var) const
{
    const Literal lit_pos = var2lit(var, false);
    const Literal lit_neg = var2lit(var, true);

    const int opv = static_cast<int>(_occ_list[lit_pos].size());
    const int onv = static_cast<int>(_occ_list[lit_neg].size());

    if (opv == 0 && onv == 0) return std::numeric_limits<int>::max();

    _process_limit.decreaseLimitBy(_occ_list[lit_pos].size());
    _process_limit.decreaseLimitBy(_occ_list[lit_neg].size());

    int spv = 0;
    for (const ClauseID c_nr : _occ_list[lit_pos]) spv += static_cast<int>(_clauses[c_nr].size());

    int snv = 0;
    for (const ClauseID c_nr : _occ_list[lit_neg]) snv += static_cast<int>(_clauses[c_nr].size());

    val_assert(opv + onv > 0);
    val_assert(spv >= 2 * opv);
    val_assert(snv >= 2 * onv);

    const int cost = opv * (snv - onv) + onv * (spv - opv) - spv - snv;
    return cost;
}

/**
 * \brief Eliminates an existential variable by resolution
 *
 * The caller has to check if elimination by resolution is allowed.
 * This is the case if (1) the eliminated variable is a gate output,
 * (2) the eliminated variable depends on at least the same variables
 * as all variable it appears with positively in the same clause,
 * (3) the same as 2, but in clauses with positive occurrence of the
 * eliminated variable.
 *
 * The function performs universal reduction on the resolvents before
 * adding them.
 * \param[in] var the variable that should be eliminated by resolution
 * \param[out] recalc_vars collection of literals whose costs have to be updated
 */
void
Formula::elimEVar(Variable var, std::unordered_set<Variable>* recalc_vars)
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    val_assert(isExistential(var));
    val_assert(isResolvable(var));
    val_assert(!recalc_vars || recalc_vars->empty());

    if (!_unit_stack.empty()) unitPropagation();

    const Literal lit_pos = var2lit(var, false);
    const Literal lit_neg = var2lit(var, true);

    for (const ClauseID clause_pos : _occ_list[lit_pos]) {
        // Ignore the optional clauses
        if (clauseDeleted(clause_pos)) continue;
        if (clauseOptional(clause_pos)) continue;

        _process_limit.decreaseLimitBy(3, _occ_list[lit_neg].size());

        for (const ClauseID clause_neg : _occ_list[lit_neg]) {
            // Ignore the optional clauses
            if (clauseDeleted(clause_pos)) break;
            if (clauseDeleted(clause_neg)) continue;
            if (clauseOptional(clause_neg)) continue;

            _process_limit.decreaseLimitBy(3, _clauses[clause_pos].size() + _clauses[clause_neg].size());
            auto resolvent = resolve(_clauses[clause_pos], _clauses[clause_neg], var);

            if (!resolvent.isTautology()) {
                // update costs for all literals in newly introduced clauses
                if (recalc_vars) {
                    auto pos = recalc_vars->begin();
                    for (Literal lit : resolvent) {
                        // We only need to recompute the costs of variables that are
                        // resolvable. However, this restriction only pays off when the
                        // check for resolvability is cheap.
                        if (isResolvable(lit2var(lit))) pos = recalc_vars->insert(pos, lit2var(lit));
                    }
                }
                addClause(std::move(resolvent));
            }
        }
    }

    _process_limit.decreaseLimitBy(2, _occ_list[lit_neg].size() + _occ_list[lit_pos].size());

    // Delete all clauses (incl. the optional ones) that contain the eliminated
    // variable
    while (!_occ_list[lit_pos].empty()) {
        // update costs for all literals in removed clauses
        const Clause& clause = _clauses[_occ_list[lit_pos].front()];

        if (recalc_vars) {
            auto pos = recalc_vars->begin();
            for (Literal lit : clause) {
                pos = recalc_vars->insert(pos, lit2var(lit));
            }
        }
        removeClause(_occ_list[lit_pos].front());
    }

    while (!_occ_list[lit_neg].empty()) {
        // update costs for all literals in removed clauses
        const Clause& clause = _clauses[_occ_list[lit_neg].front()];
        if (recalc_vars) {
            auto pos = recalc_vars->begin();
            for (Literal lit : clause) {
                pos = recalc_vars->insert(pos, lit2var(lit));
            }
        }
        removeClause(_occ_list[lit_neg].front());
    }
    removeVar(var);
    ++stat(Statistics::RESOLUTION);

    // Propagate possible new units
    unitPropagation();

    VLOG(3) << " - Length after resolution:  " << numLiterals();
}

/**
 * \brief Eliminates an existential variable by resolution if the cost is below
 * a limit
 *
 * The caller has to check if elimination by resolution is allowed.
 * This is the case if (1) the eliminated variable is a gate output,
 * (2) the eliminated variable depends on at least the same variables
 * as all variable it appears with positively in the same clause,
 * (3) the same as 2, but in clauses with positive occurrence of the
 * eliminated variable.
 *
 * The function performs universal reduction on the resolvents before
 * adding them.
 * \param[in] var the variable that should be eliminated by resolution
 * \param[in] max_cost the maximum cost up to which resolution should be
 * executed \param[out] recalc_vars collection of literals whose costs have to
 * be updated \return true if the variable has been resolved on, false otherwise
 */
bool
Formula::elimEVarLimit(const Variable var, const long int max_cost, std::unordered_set<Variable>* recalc_vars)
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    val_assert(isExistential(var));
    val_assert(isResolvable(var));
    val_assert(!recalc_vars || recalc_vars->empty());

    if (!_unit_stack.empty()) unitPropagation();

    const Literal lit_pos = var2lit(var, false);
    const Literal lit_neg = var2lit(var, true);

    long actual_cost = 0;

    // Subtract the cost of those clauses that will be deleted
    for (const ClauseID c_nr : _occ_list[lit_pos]) actual_cost -= _clauses[c_nr].size();
    for (const ClauseID c_nr : _occ_list[lit_neg]) actual_cost -= _clauses[c_nr].size();

    std::vector<Clause> to_add;
    to_add.reserve(_occ_list[lit_pos].size() * _occ_list[lit_neg].size());

    for (const ClauseID clause_pos : _occ_list[lit_pos]) {
        // Ignore the optional clauses
        if (clauseDeleted(clause_pos)) continue;
        if (clauseOptional(clause_pos)) continue;

        for (const ClauseID clause_neg : _occ_list[lit_neg]) {
            // Ignore the optional clauses
            if (clauseDeleted(clause_pos)) break;
            if (clauseDeleted(clause_neg)) continue;
            if (clauseOptional(clause_neg)) continue;

            auto resolvent = resolve(_clauses[clause_pos], _clauses[clause_neg], var);

            if (!resolvent.isTautology()) {
                universalReduction(resolvent, -1);
                actual_cost += resolvent.size();
                if (actual_cost > max_cost) return false;
                to_add.push_back(std::move(resolvent));
            }
        }
    }

    // Delete all clauses (incl. the optional ones) that contain the eliminated
    // variable
    while (!_occ_list[lit_pos].empty()) {
        // update costs for all literals in removed clauses
        const Clause& clause = _clauses[_occ_list[lit_pos].front()];

        if (recalc_vars) {
            auto pos = recalc_vars->begin();
            for (const Literal lit : clause) {
                pos = recalc_vars->insert(pos, lit2var(lit));
            }
        }
        removeClause(_occ_list[lit_pos].front());
    }

    while (!_occ_list[lit_neg].empty()) {
        // update costs for all literals in removed clauses
        const Clause& clause = _clauses[_occ_list[lit_neg].front()];
        if (recalc_vars) {
            auto pos = recalc_vars->begin();
            for (const Literal lit : clause) {
                pos = recalc_vars->insert(pos, lit2var(lit));
            }
        }
        removeClause(_occ_list[lit_neg].front());
    }

    removeVar(var);

    // add the new clauses
    for (Clause& resolvent : to_add) {
        if (recalc_vars) {
            auto pos = recalc_vars->begin();
            for (Literal lit : resolvent) {
                pos = recalc_vars->insert(pos, lit2var(lit));
            }
        }
        addClause(std::move(resolvent));
    }

    // Propagate possible new units
    if (!_unit_stack.empty()) unitPropagation();

    ++stat(Statistics::RESOLUTION);

    return true;
}

/**
 * \brief Tries to eliminate existential variables by resolution.
 *
 * First we check, which variables can be eliminated, see
 * Formula::getResolvableVariables(bool) const. Since this is
 * rather expensive if done completely, we first restrict ourselves
 * to gate detection. Only if all gate outputs (for which this
 * is feasible) have been eliminated, we switch to a complete
 * determination of resolvable variables. Variables are eliminated
 * as long as the cost of their eliminated is less than a constant
 * (currently: 0).
 * \return true if the formula was modified.
 */
bool
Formula::applyResolution()
{
    VLOG(1) << __FUNCTION__;

    ScopeTimer resolution(getTimer(WhichTimer::RESOLUTION));
    _process_limit.setLimit(PreproMethod::RESOLUTION);

    gateDependencies(DependencyOperation::ADD);

    std::size_t count = 0;

    std::unordered_set<Variable> recalc_vars;
    _variable_score.resize(maxVarIndex() + 1, 0);
    _candidates.clear();
    _candidates.resize(maxVarIndex() + 1);

    const auto resolvable_vars = getResolvableVariables();
    VLOG(3) << __FUNCTION__ << " has found " << resolvable_vars.size() << " resolvable variables.";

    for (Variable var : resolvable_vars) {
        _variable_score[var] = computeResolutionCosts(var);
        _candidates.insert(var);
    }

    while (!_candidates.empty()) {
        if (_interrupt) break;

        //        if (_process_limit.reachedLimit()) {
        //            VLOG(2) << "Terminate " << __FUNCTION__ << " due to process
        //            limit."; break;
        //        }

        const Variable next_var = _candidates.top();

        if (_variable_score[next_var] > _settings.max_resolution_cost + 300) {
            // Candidate list is sorted -> if we found a candidate which is too
            // expensive, every remaining is too expensive
            break;
        }

        if (varDeleted(next_var) || !isResolvable(next_var)) {
            // skip variables that are already deleted or not resolvable
            continue;
        }

        //        elimEVar(next_var, &recalc_vars); ++count;
        if (elimEVarLimit(next_var, _settings.max_resolution_cost, &recalc_vars)) {
            ++count;
        }

        if (stat(Statistics::RESOLUTION) % 10 == 0) fastPreprocess(false);

        if (_interrupt) break;

        // now update costs
        for (Variable var : recalc_vars) {
            if (var != next_var && isExistential(var)) {
                // If out-of-order resolution is enabled, check
                // if the variable can still be eliminated.
                const bool resolvable = isResolvable(var);
                if (!resolvable && _candidates.inHeap(var)) {
                    _candidates.remove(var);
                    continue;
                }

                if (resolvable) {
                    const int old_cost   = _variable_score[var];
                    _variable_score[var] = computeResolutionCosts(var);
                    // variable might be removed earlier from heap -> reinsert
                    if (!_candidates.inHeap(var)) {
                        _candidates.insert(var);
                    }
                    // Update heap only if new costs are lower
                    else if (old_cost > _variable_score[var]) {
                        _candidates.update(var);
                    }
                }
            }
        }
        recalc_vars.clear();
    }

    VLOG(2) << __FUNCTION__ << " eliminated " << count << " variables.";
    return count > 0;
}

/**
 * \brief Checks whether the variable `var` can be eliminated as part of an
 * implication chain
 *
 * This function eliminates the implication chain.
 * \return the replaced literal if the variable was eliminated, otherwise 0
 * \sa Formular::findImplicationChains()
 */
Literal
Formula::checkImplicationChain(const Literal lit)
{
    val_assert(minLitIndex() <= lit && lit <= maxLitIndex());

    const Variable var = lit2var(lit);

    if (!isExistential(var)) return 0;
    if (_settings.preserve_gates && _gates.isGateOutput(var)) return 0;
    if (_occ_list[lit].size() != 1) return 0;
    const Literal neg_lit = negate(lit);
    if (_implications[neg_lit].size() != 1) return 0;

    const Literal other_lit = (_implications[neg_lit].begin())->getLiteral();
    if (!dependenciesSubset(lit2var(other_lit), var)) return 0;

    removeClause(_occ_list[lit].front());
    replaceLiteral(neg_lit, other_lit);

    val_assert(_occ_list[lit].empty());
    val_assert(_occ_list[neg_lit].empty());

    ++stat(Statistics::IMPLICATION_CHAINS);
    VLOG(3) << __FUNCTION__ << "() replaced literal " << lit2dimacs(neg_lit) << " by literal " << lit2dimacs(other_lit);

    unitPropagation();
    return other_lit;
}

/**
 * \brief Finds and eliminates implication chains.
 *
 * This is a special case of variable elimination by resolution,
 * where the eliminated variable occurs positively (or negatively)
 * only in a single clause.
 * \return true if the formula was modified.
 */
bool
Formula::findImplicationChains()
{
    VLOG(1) << __FUNCTION__;

    ScopeTimer impl_chain(getTimer(WhichTimer::IMPLICATION_CHAINS));

    const std::size_t old_stat_impl_chains = stat(Statistics::IMPLICATION_CHAINS);

    for (Literal lit = minLitIndex(); lit <= maxLitIndex(); ++lit) {
        if (_interrupt) break;
        checkImplicationChain(lit);
    }

    VLOG(2) << __FUNCTION__ << "() replaced " << (stat(Statistics::IMPLICATION_CHAINS) - old_stat_impl_chains)
            << " variables.";

    return stat(Statistics::IMPLICATION_CHAINS) > old_stat_impl_chains;
}

}  // end namespace hqspre
