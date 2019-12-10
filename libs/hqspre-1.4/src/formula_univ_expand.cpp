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
#include <limits>
#include <set>
#include <stack>
#include <utility>
#include <vector>

#define ELPP_STL_LOGGING
#include <easylogging++.hpp>
#include <solver/antom.h>
#include <solver/antombase.h>

#include "aux.hpp"
#include "clause.hpp"
#include "exceptions.hpp"
#include "formula.hpp"
#include "literal.hpp"
#include "prefix.hpp"
#include "timer.hpp"
#include "varheap.hpp"


/**
 * \file formula_univ_expand.cpp
 * \brief Implementation of universal expansion
 * \author Ralf Wimmer
 * \author Sven Reimer
 * \date 2016-17
 */

namespace hqspre {

/**
 * \brief Tries to estimate the costs of universal expansion
 *
 * The function returns as the first entry of the return value the
 * number of existential variables that need to be doubled. The second
 * entry is the number of additional clauses after expansion.
 * These numbers are over-approximations as some of the variables might
 * be unit or clauses tautologies.
 */
std::pair<int, int> Formula::computeExpansionCosts(const Variable uvar) const
{
    int result_vars = 0; // number of additional existential variables
    int result_clauses = 0; // number of additional clauses

    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (isExistential(var) && depends(var, uvar)) ++result_vars;
    }

    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (clauseDeleted(c_nr) || clauseOptional(c_nr)) continue;
        bool contains_uvar = false;
        bool to_copy = false;

        for (Literal lit: _clauses[c_nr]) {
            const Variable var = lit2var(lit);
            if (var == uvar) contains_uvar = true;
            else if (isExistential(var) && depends(var, uvar)) to_copy = true;
        }
        if (to_copy && !contains_uvar) result_clauses += 1;
    }

    return std::make_pair(result_vars, result_clauses);
}

void Formula::markTransitiveUnits(std::stack<Literal>& units, std::vector<bool>& marked) const
{
    while (!units.empty()) {
        const Literal top = units.top();
        units.pop();
        for (const BinaryClause succ: _implications[top]) {
            const Literal lit = succ.getLiteral();
            if (!marked[lit]) {
                marked[lit] = true;
                units.push(lit);
            }
        }
    }
}

/**
 * \brief Tries to estimate the costs of universal expansion
 *
 * The function returns difference of total literals of the formula
 * This number is an over-approxmimation as further units and subsumptions
 * are not taken into account.
 */
long int Formula::computeExpansionCosts2(const Literal ulit, const std::set<Variable>& pseudo_deps)
{
    long int result_lits = 0; // number of additional literals

    const Variable uvar = lit2var(ulit);
    // mark new units
    std::vector<bool> new_units_pos(maxLitIndex() + 1, false);
    std::vector<bool> new_units_neg(maxLitIndex() + 1, false);

    std::stack<Literal> front;

    // Binary clause becomes unit after expansion, mark new units
    for (const BinaryClause bin_clause : _implications[ulit] ) {
        new_units_neg[bin_clause.getLiteral()] = true;
        front.push(bin_clause.getLiteral());
    }

    // Now mark also transitive units
    markTransitiveUnits(front, new_units_neg);

    for (const BinaryClause bin_clause : _implications[negate(ulit)] ) {
        new_units_pos[bin_clause.getLiteral()] = true;
        front.push(bin_clause.getLiteral());
    }

    // Now mark also transitive units
    markTransitiveUnits(front, new_units_pos);

    // Now determine and mark pseudo dependencies
    std::vector<bool> is_pseudo_dep(maxVarIndex() + 1, false);

    for (const Variable var : pseudo_deps) {
        is_pseudo_dep[var] = true;
    }

    bool to_copy = false;
    bool contains_uvar_pos = false;
    bool contains_uvar_neg = false;
    bool pos_satisfied = false;
    bool neg_satisfied = false;
    std::size_t pos_clause_size = 0;
    std::size_t neg_clause_size = 0;

    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (clauseDeleted(c_nr)) continue;
        if (clauseOptional(c_nr)) {
            result_lits -= _clauses[c_nr].size();
            continue;
        }

        to_copy = false;
        contains_uvar_pos = false;
        contains_uvar_neg = false;
        pos_satisfied = false;
        neg_satisfied = false;
        pos_clause_size = _clauses[c_nr].size();
        neg_clause_size = pos_clause_size;

        // Assume clause will be deleted
        result_lits -= pos_clause_size;

        _process_limit.decreaseLimitBy(6, _clauses[c_nr].size());

        for (const Literal lit: _clauses[c_nr]) {
            const Variable var = lit2var(lit);
            if (var == uvar) {
                contains_uvar_pos = isPositive(lit);
                contains_uvar_neg = !contains_uvar_pos;
                continue;
            }
            // we have to copy each clause where an existential literal exists which depends on "uvar"
            else if (isExistential(var) && depends(var, uvar) && !is_pseudo_dep[var]) {
                to_copy = true;
            }

            // lit will be unit in positive part of expansion
            // -> clause will be satisfied in positive expansion part
            if (new_units_pos[lit]) {
                pos_satisfied = true;
            }
            // ~lit will be unit in positive expansion
            // -> literal will be deleted in positive expansion part
            else if (new_units_pos[negate(lit)]) {
                --pos_clause_size;
            }

            // lit will be unit in negative part of expansion
            // -> clause will be satisfied in negative expansion part
            if (new_units_neg[lit]) {
                neg_satisfied = true;
            }
            // ~lit will be unit in negative expansion
            // -> literal will be deleted in negative expansion part
            else if (new_units_neg[negate(lit)]) {
                --neg_clause_size;
            }
        }

        // expanding universal var will be deleted in positive expansion
        if (contains_uvar_pos) {
            --pos_clause_size;
            // Binary clause containing expanding var will be deleted after expansion
            if (pos_clause_size == 1 || pos_satisfied) {
                pos_clause_size = 0;
            }
            result_lits += pos_clause_size;
        }
        // expanding universal var will be deleted in negative expansion
        else if (contains_uvar_neg) {
            --neg_clause_size;
            // Binary clause containing expanding var will be deleted after expansion
            if (neg_clause_size == 1 || neg_satisfied) {
                neg_clause_size = 0;
            }
            result_lits += neg_clause_size;
        }
        // clause does not contain any universal literal
        // -> will be possibly copied
        else {
            // Reset clause size in case clause will be satisfied
            if (pos_satisfied) pos_clause_size = 0;
            if (neg_satisfied) neg_clause_size = 0;

            result_lits += pos_clause_size;
            // We have to copy the clause
            if (to_copy && !contains_uvar_neg) {
                result_lits += neg_clause_size;
            }
        }
    }

    return result_lits;
}

/**
 * \brief Applies universal expansion to those universal variables which can cheaply be eliminated.
 * \return True if the formula was modified.
 */
bool Formula::applyUniversalExpansion()
{
    VLOG(1) << __FUNCTION__;

    ScopeTimer univ_expansion(getTimer(WhichTimer::UNIV_EXPANSION));

    std::size_t count = 0;
    if (numUVars() <= 6) {
        const std::size_t old_num_clauses = numClauses();
        for (Variable uvar = minVarIndex(); uvar <= maxVarIndex(); ++uvar) {
            if (_interrupt) break;

            if (isUniversal(uvar)) {
                elimUVar(var2lit(uvar));
                ++count;
            }
            if (numClauses() > 2 * old_num_clauses) break;
        }

        VLOG(2) << __FUNCTION__ << "() expanded " << count << " universal variables.";
        return count > 0;
    }

    if (_prefix->type() == PrefixType::QBF && count < 30) {
        const std::size_t old_size = numLiterals();
        std::set<Variable> pseudo_deps;
        for (int level = static_cast<int>(_qbf_prefix->getMaxLevel()); level >= 0; --level) {
            if (_interrupt) break;
            if (_qbf_prefix->getLevelQuantifier(level) != VariableStatus::UNIVERSAL) continue;

            const auto& block = _qbf_prefix->getVarBlock(level);
            if (block.size() > 20) continue;

            while (!block.empty()) {
                if (static_cast<double>(numLiterals()) > 1.5 * static_cast<double>(old_size)) break;
                int next_var = -1;
                int min_cost = std::numeric_limits<int>::max();
                int varcost = 0;
                for (Variable uvar: block) {
                    if (_interrupt) break;
                    auto cost = computeExpansionCosts(uvar);
                    sstdQuadDep(uvar, true, true, &pseudo_deps);
                    cost.first -= static_cast<int>(pseudo_deps.size());
                    pseudo_deps.clear();
                    varcost = cost.first;
                    if (cost.second < min_cost) { next_var = uvar; min_cost = cost.second; }
                }
                if (_interrupt) break;

                if (next_var > 0 && varcost < 20 && min_cost < 100) {
                    elimUVar(var2lit(next_var));
                    applyResolution();
                    ++count;
                } else break;
            }
        }

        VLOG(2) << __FUNCTION__ << "() expanded " << count << " universal variables.";
        return count > 0;
    }

    return false;
}

bool Formula::applyUniversalExpansion2()
{
    VLOG(1) << __FUNCTION__;

    ScopeTimer univ_expansion(getTimer(WhichTimer::UNIV_EXPANSION));
    if (_prefix->type() == PrefixType::QBF) {

        std::size_t count = 0;

        _candidates.resize(maxVarIndex() + 1);

        // All universal variables and its pseudo dependencies
        std::vector<std::set< Variable > > var_candidates(maxVarIndex()+1);

        // maximum formula size after expansion
        if (_settings.max_expansion_size == 0) {
            // set fixed value in first call.
            _settings.max_expansion_size = numLiterals() << 1;
        }


        while (true) {
            if (_interrupt) break;

            if (numUVars() <= 6) {
                const std::size_t old_num_clauses = numClauses();
                for (Variable uvar = minVarIndex(); uvar <= maxVarIndex(); ++uvar) {
                    if (isUniversal(uvar)) {
                        elimUVar(var2lit(uvar));
                        ++count;
                    }
                    if (numClauses() > 3 * old_num_clauses) break;
                }

                if (numUVars() == 0) {
                    return count > 0;
                }
            }

            std::size_t current_count = count;
            std::size_t literals = numLiterals();
            std::set<Variable> pseudo_deps;
            // Now try to remove a whole quantifier level
            for (int level = static_cast<int>(_qbf_prefix->getMaxLevel()); level >= 0; --level) {
                // Exit routine if we expand our maximal expansion size
                if (literals > _settings.max_expansion_size ) {
                    VLOG(2) << __FUNCTION__ << " expanded " << count << " universal variables.";
                    return count > 0;
                }

                if (_qbf_prefix->getLevelQuantifier(level) != VariableStatus::UNIVERSAL) continue;

                const auto& block = _qbf_prefix->getVarBlock(level);
                if (block.size() > 20) continue;
                VLOG(3) << __FUNCTION__ << " trying to remove quantifier level " << level << " with " << block.size() << " variables";

                const std::size_t old_size = literals;
                while (!block.empty()) {
                    int next_var = -1;
                    long int min_cost = std::numeric_limits<int>::max();
                    for (Variable uvar: block) {
                        // sstdQuadDep(uvar, true, true, &pseudo_deps);
                        const long int cost = computeExpansionCosts2(var2lit(uvar, false), pseudo_deps);
                        pseudo_deps.clear();
                        if (cost < min_cost) { next_var = uvar; min_cost = cost; }
                    }

                    elimUVar(var2lit(next_var, false));
                    if (_interrupt) break;
                    if (applyResolution()) {
                        fastPreprocess();
                    }

                    literals = numLiterals();
                    ++count;
                    if (static_cast<double>(literals) > 1.5 * static_cast<double>(old_size)) break;
                }

                if (block.empty()) {
                    updateVars();
                    level = static_cast<int>(_qbf_prefix->getMaxLevel());
                }
            }
            // We have done nothing -> break routine
            if (count == current_count) break;
            if (_interrupt) break;
        }
        VLOG(2) << __FUNCTION__ << " expanded " << count << " universal variables.";

        return count > 0;
    }
    return false;
}


std::vector<Variable> Formula::computeVarElimSet()
{
    antom::Antom antom;
    std::vector<Variable> var1_minus_var2;
    std::vector<Variable> var2_minus_var1;
    antom.SetMaxIndex(2 * maxVarIndex());
    Variable current = maxVarIndex() + 1;
    std::vector<Literal> bin_clause(2, 0);

    // Add hard constraints: eliminate all 2-cycles
    for (Variable var1 = minVarIndex(); var1 < maxVarIndex(); ++var1) {
        if (!isExistential(var1) || _prefix->inRMB(var1)) continue;
        for (Variable var2 = var1 + 1; var2 <= maxVarIndex(); ++var2) {
            if (!isExistential(var2) || _prefix->inRMB(var2)) continue;
            if (dependenciesSubset(var1, var2) || dependenciesSubset(var2, var1)) continue;
            var1_minus_var2.clear();
            var2_minus_var1.clear();
            const auto& dep1 = _dqbf_prefix->getDependencies(var1);
            const auto& dep2 = _dqbf_prefix->getDependencies(var2);
            two_sided_difference(
                                 dep1.cbegin(), dep1.cend(),
                                 dep2.cbegin(), dep2.cend(),
                                 std::back_inserter(var1_minus_var2), std::back_inserter(var2_minus_var1)
                                );
            val_assert(!var1_minus_var2.empty() && !var2_minus_var1.empty());
            const Variable edge1_2 = current++;
            const Variable edge2_1 = current++;

            // Eliminate either "edge1" or "edge2"
            bin_clause[0] = var2lit(edge1_2, false);
            bin_clause[1] = var2lit(edge2_1, false);
            antom.AddClause(bin_clause);

            // edge1
            for (const Variable uvar: var1_minus_var2) {
                bin_clause[0] = var2lit(edge1_2, true);
                bin_clause[1] = var2lit(uvar, false);
                antom.AddClause(bin_clause);
            }

            // edge2
            for (const Variable uvar: var2_minus_var1) {
                bin_clause[0] = var2lit(edge2_1, true);
                bin_clause[1] = var2lit(uvar, false);
                antom.AddClause(bin_clause);
            }
        }
    }

    // Add soft clauses
    bin_clause.resize(1);
    for (Variable uvar = minVarIndex(); uvar <= maxVarIndex(); ++uvar) {
        if (isUniversal(uvar)) {
            bin_clause[0] = var2lit(uvar, true);
            antom.AddSoftClause(bin_clause);
            bin_clause.resize(1);
        }
    }

    // Solve and determine optimal solution:
    std::int64_t optimum = 0;
    const auto result = antom.MaxSolve(optimum);
    if (result != ANTOM_SAT) {
        throw ElimSetException("Antom could not find a solution for universal expansion.");
    }

    const auto& model = antom.Model();
    std::vector<Variable> to_eliminate;
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (isUniversal(var) && isPositive(model[var])) {
            to_eliminate.push_back(var);
        }
    }

    return to_eliminate;
}


bool Formula::applyUniversalExpansionDQBF()
{
    val_assert(_prefix->type() == PrefixType::DQBF);
    val_assert(_dqbf_prefix);

    VLOG(1) << __FUNCTION__;

    const bool force_eliminate = (_settings.convert_to_qbf == DQBFtoQBFmethod::VAR_ELIM);

    ScopeTimer univ_expansion(getTimer(WhichTimer::UNIV_EXPANSION));

    std::vector<Variable> to_eliminate = computeVarElimSet();
    VLOG(2) << __FUNCTION__ << ": We need to eliminate " << to_eliminate.size() << " univ vars out of " << numUVars() << ": " << to_eliminate;

    const auto old_size = numLiterals();
    std::set<Variable> pseudo_deps;
    std::size_t count = 0;

    while (!to_eliminate.empty()) {
        int next_pos = -1;
        long int min_cost = std::numeric_limits<long int>::max();
        for (std::size_t pos = 0; pos < to_eliminate.size(); ++pos) {
            const Variable uvar = to_eliminate[pos];
            if (!isUniversal(uvar)) continue;
            sstdQuadDep(uvar, true, true, &pseudo_deps);
            const long int cost = computeExpansionCosts2(var2lit(uvar, false), pseudo_deps);
            pseudo_deps.clear();
            if (cost < min_cost) { next_pos = static_cast<int>(pos); min_cost = cost; }
        }

        if (next_pos >= 0 && (force_eliminate || to_eliminate.size() <= 3 || static_cast<double>(min_cost) < 0.1 * static_cast<double>(old_size))) {
            std::swap(to_eliminate[next_pos], to_eliminate.back());
            const Variable uvar = to_eliminate.back();
            to_eliminate.pop_back();

            elimUVar(var2lit(uvar));
            ++count;
            if (_interrupt) break;
            applyResolution();
            if (_interrupt) break;
            if (static_cast<double>(numLiterals()) > 1.1 * static_cast<double>(old_size)) break;
        } else break;
    }

    if (_prefix->type() == PrefixType::DQBF && isQBF()) {
        convertToQBF();
        VLOG(2) << __FUNCTION__ << " formula is now a QBF after expanding " << stat(Statistics::UNIV_EXPANSION) << " variables.";
    }

    return count > 0;
}


/**
 * \brief Eliminates a universal variable by expansion.
 * \param lit a literal of the universal variable to be eliminated. The polarity of `lit` determines which variable copy appears in which cofactor.
 * \note This function creates new variables which serve as
 *       copies of those existential variables which depend
 *       upon `lit2var(lit)`.
 * \todo If \f$x\f$ is the eliminated universal variable,
 *       \f$y\f$ an existential one, and there is an implication
 *       like \f$x\rightarrow y\f$, one can replace \f$y\f$ in
 *       the 1-cofactor by 1. The remaining three cases
 *       are similar. We could look for such implications using
 *       a SAT-solver or check if the implication is a blocked
 *       (binary clause).
 */
void Formula::elimUVar(const Literal lit)
{
    enum Occurrence {
        NOT_AT_ALL = 0,
        POSITIVE = 1,
        NEGATIVE = 2
    };

    const Variable var = lit2var(lit);

    val_assert(_prefix);
    val_assert_msg(isUniversal(var),    "You can only expand universal variables!");

    VLOG(2) << __FUNCTION__ << "(" << var << ") -> size before: " << numLiterals();

    if (!_unit_stack.empty()) unitPropagation();

    // Remove all pseudo-dependencies of the expanded variable.
    std::set<Variable> pseudo_deps;
    sstdQuadDep(var, true, true, &pseudo_deps);

    const Literal u_pos = lit;
    const Literal u_neg = negate(lit);

    // Clauses that will be deleted
    std::vector<ClauseID> to_delete;
    // Clauses that will be added
    std::vector<Clause::ClauseData> to_add;
    // Mapping the dependent literals to their copies
    std::vector<Literal> lit_map(maxLitIndex() + 1, 0);

    // Create a copy of each dependent var and store it in 'var_map'
    const Variable old_maxVarIndex = maxVarIndex();
    for (Variable e_var = minVarIndex(); e_var <= old_maxVarIndex; ++e_var) {
        if (!isExistential(e_var) || !depends(e_var, var) || pseudo_deps.find(e_var) != pseudo_deps.cend()) {
            lit_map[ var2lit(e_var, false) ] = var2lit(e_var, false);
            lit_map[ var2lit(e_var, true) ] = var2lit(e_var, true);
        } else {
            const Variable e_var_copy = copyVar(e_var);
            val_assert(isExistential(e_var_copy));
            lit_map[var2lit(e_var, false)] = var2lit(e_var_copy, false);
            lit_map[var2lit(e_var, true) ] = var2lit(e_var_copy, true);
        }

        if (_interrupt) return;
    }

    VLOG_IF(!pseudo_deps.empty(), 2) << "Dependency scheme removed " << pseudo_deps.size() << " out of " << (lit_map.size() / 2) + pseudo_deps.size() << " dependencies.";

    Clause::ClauseData new_clause;
    new_clause.reserve(50);

    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (clauseDeleted(c_nr)) continue;

        Occurrence contains_u = NOT_AT_ALL;
        bool contains_dependent_var = false;

        // Check if the clause contains the eliminated universal variables
        // or any variable which needs to be copied.
        for (const Literal cl_lit: _clauses[c_nr]) {
            if (cl_lit == u_pos) contains_u = POSITIVE;
            else if (cl_lit == u_neg) contains_u = NEGATIVE;
            else if (isExistential(lit2var(cl_lit)) && lit_map[cl_lit] != cl_lit) {
                contains_dependent_var = true;
            }
        }

        if (contains_u == NOT_AT_ALL && !contains_dependent_var) {
            // clause remains unchanged
            continue;
        } else if (clauseOptional(c_nr)) {
            // delete optional clauses that would not remain unchanged
            to_delete.push_back(c_nr);
            continue;
        } else if (contains_u == POSITIVE) {
            // clause[c_nr] contains 'var' positively
            // -> only remove u_pos
            removeLiteral(c_nr, u_pos);
        } else if (contains_u == NEGATIVE) {
            // clause contains 'var' negatively
            // -> remove ~var and replace the dependent variables by their copies
            for (const Literal cl_lit: _clauses[c_nr]) {
                if (cl_lit == u_neg) continue;
                new_clause.push_back(lit_map[cl_lit]);
            }

            to_add.push_back(std::move(new_clause));
            new_clause.clear();
            to_delete.push_back(c_nr);
        } else if (contains_u == NOT_AT_ALL && contains_dependent_var) {
            // clause does not contain 'var', but exist. vars that depend upon it
            // -> Add a copy of the clause with replaced variables.
            for (const Literal cl_lit: _clauses[c_nr]) {
                new_clause.push_back(lit_map[cl_lit]);
            }

            to_add.push_back(std::move(new_clause));
            new_clause.clear();
        }
    } // end for c_nr


    // Remove old and add new clauses
    for (const ClauseID c_nr: to_delete) removeClause(c_nr);
    for (auto& clause: to_add) addClause(std::move(clause));

    // Delete the expanded variable
    removeVar(var);
    ++stat(Statistics::UNIV_EXPANSION);

    fastPreprocess();

    VLOG(2) << "Size after univ. expansion: " << numLiterals();

    if (numClauses() == 0) throw SATException("No clauses left after universal expansion");
}


} // end namespace hqspre
