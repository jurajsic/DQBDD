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

#include "formula.hpp"
#include "formula.ipp"

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <set>
#include <stack>
#include <utility>
#include <vector>

#include <easylogging++.hpp>
#include "auxil.hpp"
#include "clause.hpp"
#include "gate.hpp"
#include "literal.hpp"
#include "prefix.hpp"
#include "timer.hpp"
#include "varheap.hpp"

#define TEST




/**
 * \file formula.cpp
 * \brief Implementation of basic functions for variable, clause, and dependency
 * management \author Ralf Wimmer \date 12/2015-03/2016
 */

namespace hqspre {


/**
 * \brief Performs preprocessing on the formula.
 */
void
Formula::preprocess()
{
    ScopeTimer prepro(getTimer(WhichTimer::TOTAL_TIME));

    val_assert(checkConsistency());

    if (_interrupt) {
        return;
    }

    _blocked_candidates.resize(_clauses.size());

    // Do initial unit and pure propagation
    fastPreprocess();

	/*
	// try to perform forkSplitting
	if (!isQBF() && _settings.enableFork) {
	std::cout << "entering forkext" << std::endl;
		_settings.forkExt = forkExtension();
	}

	if (_settings.forkExt) {
	std::cout << "formula is disjunct" << std::endl;
		_settings.findHidden = resolveAfterForkSplit();
	}
	
	if (_settings.findHidden) {
		findHiddenEquivalences();
		solveSAT();
	}
	*/


    val_assert(checkConsistency());

    if (_interrupt) {
        return;
    }

    determineGates(true, true, true, _settings.semantic_gates);
    gateDependencies(DependencyOperation::ADD);

    val_assert(checkConsistency());

    if (_interrupt) {
        return;
    }

    if (_settings.equiv_gates) {
        if (findEquivalentGates()) {
            fastPreprocess();
        }
    }

    val_assert(checkConsistency());

    bool again       = true;
    auto temp_hidden = _settings.hidden;

    while (again && stat(Statistics::PREPRO_LOOPS) < _settings.max_loops) {
        if (_interrupt) {
            return;
        }
        ++stat(Statistics::PREPRO_LOOPS);

        VLOG(1) << "Preprocessing loop #" << stat(Statistics::PREPRO_LOOPS) << " started.";

        again = false;
        val_assert(checkConsistency());

        // Substitution
        if (_settings.substitution && stat(Statistics::PREPRO_LOOPS) <= _settings.max_substitution_loops
            && applySubstitution()) {
            again = true;
            fastPreprocess();
            val_assert(checkConsistency());
        }

        // Subsumption checks
        if (_settings.subsumption && removeSubsumedClauses()) {
            again = true;
            fastPreprocess();
            val_assert(checkConsistency());
            if (_interrupt) {
                return;
            }
        }

        // Always skip hidden literals in first run
        if (stat(Statistics::PREPRO_LOOPS) == 1) {
            _settings.hidden = 0u;
        }

        // Blocked clause and friends (tautology/subsumption) elimination
        val_assert(checkConsistency());
        if ((_settings.bce > 0 || _settings.hse) && removeBlockedAndFriends()) {
            val_assert(checkConsistency());
            again = true;
            fastPreprocess();
            if (_interrupt) {
                return;
            }
        }
        if (stat(Statistics::PREPRO_LOOPS) == 1) {
            _settings.hidden = temp_hidden;
        }

        if (_settings.hec && findHiddenEquivAndContraDefinitions()) {
            again = true;
            fastPreprocess();
            val_assert(checkConsistency());
            if (_interrupt) {
                return;
            }
        } else {
            _settings.hec = false;
        }

        // Self-subsuming resolution
        if (_settings.self_subsumption && selfSubsumingResolution()) {
            again = true;
            fastPreprocess();
            val_assert(checkConsistency());
            if (_interrupt) {
                return;
            }
        }

        // Use vivification (Piette et al.) to shorten clauses
        if (_settings.vivify) {
            while (applyVivification() && _settings.vivify_fp) {
                again = true;
                fastPreprocess();
                val_assert(checkConsistency());
                if (_interrupt) {
                    return;
                }
            }
        }

        // Unit propagation lookahead
        if (_settings.upla) {
            if (!upla()) {
                _settings.upla = false;
            }
            val_assert(checkConsistency());
            if (_interrupt) {
                return;
            }
        }

        // Elimination of implication chains
        if (_settings.impl_chains && findImplicationChains()) {
            again = true;
            fastPreprocess();
            val_assert(checkConsistency());
            if (_interrupt) {
                return;
            }
        }

        // Variable elimination by resolution
        if (_settings.resolution && applyResolution()) {
            again = true;
            fastPreprocess();
            val_assert(checkConsistency());
            if (_interrupt) {
                return;
            }
        }

        // Subsumption checks
        if (_settings.subsumption && removeSubsumedClauses()) {
            again = true;
            fastPreprocess();
            val_assert(checkConsistency());
            if (_interrupt) {
                return;
            }
        }

        // Detecting contradictions
        if (_settings.contradictions && findContradictions()) {
            again = true;
            fastPreprocess();
            val_assert(checkConsistency());
            if (_interrupt) {
                return;
            }
        }

        // Finding constant variables (a.k.a. backbones) by SAT calls
        if (stat(Statistics::PREPRO_LOOPS) == 1) {
            if (_settings.sat_const && findConstantsBySAT()) {
                again = true;
                fastPreprocess();
                val_assert(checkConsistency());
                if (_interrupt) {
                    return;
                }
            }
            if (_settings.bia && addBlockedImplications()) {
                again = true;
                fastPreprocess();
                val_assert(checkConsistency());
                if (_interrupt) {
                    return;
                }
            }

            // Finding implications (= implied binary clauses) by SAT calls
            if (_settings.sat_impl && findImplicationsBySAT()) {
                again = true;
                fastPreprocess();
                val_assert(checkConsistency());
                if (_interrupt) {
                    return;
                }
            }

            determineGates(true, true, true, false);
            removeOptionalClauses();
            if (_interrupt) {
                return;
            }
        } else if (static_cast<unsigned int>(stat(Statistics::PREPRO_LOOPS)) <= _settings.max_substitution_loops) {
            determineGates(true, true, true, false);
            removeOptionalClauses();
        }

        // Incomplete SAT/UNSAT checks
        if (stat(Statistics::PREPRO_LOOPS) == 1 && _settings.sat_incomplete) {
            checkSAT();
            checkUNSAT(_settings.num_random_patterns);
        }

        // Apply universal expansion
        if (_settings.univ_expand == 1 && _qbf_prefix != nullptr && applyUniversalExpansion()) {
            again = true;
            val_assert(checkConsistency());
            if (_interrupt) {
                return;
            }
        } else if (_settings.univ_expand == 2 && _qbf_prefix != nullptr && applyUniversalExpansion2()) {
            again = true;
            val_assert(checkConsistency());
            if (_interrupt) {
                return;
            }
        } else if (_settings.univ_expand > 0) {
            if (_dqbf_prefix != nullptr && applyUniversalExpansionDQBF()) {
                again = true;
                val_assert(checkConsistency());
                if (_interrupt) {
                    return;
                }
            }
        }

        if (numUVars() == 0 && _settings.sat_decide) {
            checkSAT();
        }

        // Subsumption checks
        if (_settings.subsumption && removeSubsumedClauses()) {
            again = true;
            fastPreprocess();
            val_assert(checkConsistency());
            if (_interrupt) {
                return;
            }
        }

        // Incomplete SAT/UNSAT checks
        if (stat(Statistics::PREPRO_LOOPS) == 1 && _settings.sat_incomplete) {
            checkSAT();
            if (_interrupt) {
                return;
            }
            checkUNSAT(_settings.num_random_patterns);
            if (_interrupt) {
                return;
            }
        }

        updateVars();
        if (_interrupt) {
            return;
        }
        VLOG(1) << "Preprocessing loop #" << stat(Statistics::PREPRO_LOOPS) << " finished.";
        

    }  // end main prepro loop

    // If the formula is propositional, but could not be solved
    // by the SAT solver within the given time limit, abort the
    // preprocessing phase.
    if (numUVars() == 0 && _settings.sat_decide) {
        checkSAT();
        if (_interrupt) {
            return;
        }
    }

    val_assert(checkConsistency());

}

/**
 * \brief Performs preprocessing on the formula.
 * Does not touch any gate defining variables and clauses
 * \return The set of found gates
 */
std::vector<Gate>
Formula::preprocessGates()
{
    _settings.preserve_gates = true;
    preprocess();
    determineGates(true, true, true, false);
    return _gates.getGates();
}

/**
 * \brief Applies unit propagation, equivalance checking, pure literal detection
 * and implication chains
 *
 * Should be called after an expensive preprocessor method for quick and cheap
 * detection of new units. \param until_fixedpoint If true, then the operations
 * are called until no more changes happen.
 */
void
Formula::fastPreprocess(const bool until_fixedpoint)
{
    const std::size_t old_stat_unit        = stat(Statistics::UNIT);
    const std::size_t old_stat_pure        = stat(Statistics::PURE);
    const std::size_t old_stat_equiv       = stat(Statistics::EQUIV_LITS);
    const std::size_t old_stat_impl_chains = stat(Statistics::IMPLICATION_CHAINS);
    bool              changes              = false;
    std::size_t       iter                 = 0;

    do {
        changes = false;
        ++iter;
        if (unitPropagation()) {
            changes = true;
        }
        if (findEquivalences()) {
            changes = true;
        }
        if (findPure()) {
            changes = true;
        }
    } while (changes && until_fixedpoint && iter < _settings.max_fast_loops);

    VLOG_IF(stat(Statistics::UNIT) > old_stat_unit || stat(Statistics::PURE) > old_stat_pure
                || stat(Statistics::EQUIV_LITS) > old_stat_equiv
                || stat(Statistics::IMPLICATION_CHAINS) > old_stat_impl_chains,
            2)
        << __FUNCTION__ << " detected " << (stat(Statistics::UNIT) - old_stat_unit) << " units, "
        << (stat(Statistics::PURE) - old_stat_pure) << " pures, " << (stat(Statistics::EQUIV_LITS) - old_stat_equiv)
        << " equivalent literals, and " << (stat(Statistics::IMPLICATION_CHAINS) - old_stat_impl_chains)
        << " implication chains.";
}

/**
 * \brief Checks which variables still occur in the matrix.
 *
 * Marks all variables which do not occur in the matrix as deleted.
 * For QBFs, subsequent blocks with the same quantifier are merged.
 */
void
Formula::updateVars()
{
    val_assert(_prefix);
    val_assert(_unit_stack.empty());

    // Determine which variables actually occur in the formula.
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (!varDeleted(var) && _occ_list[var2lit(var, false)].empty() && _occ_list[var2lit(var, true)].empty()) {
            removeVar(var);
        }
    }

    _prefix->updateVars();
}

/**
 * \brief Returns the accumulated number of dependencies of all existential
 * variables.
 */
std::size_t
Formula::numDependencies() const noexcept
{
    val_assert(_prefix);

    return _prefix->numDependencies();
}

/**
 * \brief Checks if `var1`'s dependencies are a subset of `var2`'s dependencies.
 *
 * If `var1` is universal, this function returns whether `var2` depends on
 * `var1`. If both variables are existential, this function checks whether the
 * dependency set of `var1` is a subset of `var2`'s dependencies.
 *
 * \note `var2` has to be existential.
 */
bool
Formula::dependenciesSubset(const Variable var1, const Variable var2) const
{
    val_assert(_prefix);

    return _prefix->dependenciesSubset(var1, var2);
}

/**
 * \brief add both implications for a binary clause "lit1 lit2"
 *
 * \note performs subsumption check: if the binary clause already exists, the
 * clause with ID "c_nr" will be removed and the original one will be kept
 */
void
Formula::addImplications(const Literal lit1, const Literal lit2, const ClauseID c_nr)
{
    const int id = hasImplication(negate(lit1), lit2);
    if (id >= 0) {
        // found implication has same id -> nothing to do
        if (id == static_cast<int>(c_nr)) {
            return;
        }

        // the binary already exists -> remove new one completely
        // due to current datastructure also the original binary clause will be
        // deleted
        // => re-add implications of duplicated clause already existing in database
        removeClause(c_nr);
        addImplication(negate(lit1), lit2, static_cast<ClauseID>(id));
        addImplication(negate(lit2), lit1, static_cast<ClauseID>(id));
    } else {
        // new implications -> just add
        addImplication(negate(lit1), lit2, c_nr);
        addImplication(negate(lit2), lit1, c_nr);
    }
    _implications_added = true;
}


/**
 * \brief Inserts a clause into the occurrence lists and implication lists.
 *
 * \note Only for internal use.
 * \param c_nr the number of the newly created clause
 * \param check_subsumption if true, it is checked whether the clause subsumes
 * other clauses in database \return the clause ID if the clause was actually
 * added; -2 if the clause was a tautology; -1 if the clause was unit.
 */
int
Formula::addClauseToLists(const ClauseID c_nr, const bool check_subsumption)
{
    val_assert(c_nr <= maxClauseIndex());
    val_assert(!clauseDeleted(c_nr));

    auto& clause = _clauses[c_nr];

#ifndef NDEBUG
    for (const Literal lit : clause) {
        LOG_IF(varDeleted(lit2var(lit)), ERROR)
            << "Variable " << lit2var(lit) << " is deleted in clause " << clause << '\n';
        val_assert_msg(!varDeleted(lit2var(lit)), "Trying to add a clause with a deleted variable");
    }
#endif

    if (clause.isTautology()) {
        VLOG(3) << __FUNCTION__ << "(" << c_nr << ") Ignoring tautology " << clause;

        clause.clear();
        clause.setStatus(ClauseStatus::DELETED);
        _deleted_clause_numbers.push_back(c_nr);
        return -2;
    }

    // Perform universal reduction
    if (_settings.univ_reduction) {
        universalReduction(clause);
    }

    _clause_sizes[c_nr] = clause.size();

    if (clause.empty()) {
        throw UNSATException("Trying to add an empty clause!");
        return -2;
    }

    if (clause.size() == 1) {
        pushUnit(clause[0], PureStatus::UNIT);
        clause.clear();
        clause.setStatus(ClauseStatus::DELETED);
        _deleted_clause_numbers.push_back(c_nr);
        return -1;
    }

    // fill occurrence lists
    for (const Literal lit : clause) {
        _occ_list[lit].push_back(c_nr);
    }

    if (check_subsumption) {
        const auto num_subsumed = isBackwardSubsuming(_clauses[c_nr], static_cast<int>(c_nr));
        if (num_subsumed > 0) {
            _clauses[c_nr].setStatus(ClauseStatus::MANDATORY);
            VLOG(2) << __FUNCTION__ << "(" << _clauses[c_nr] << ") removed " << num_subsumed << " subsumed clauses.";
        }

    }

    // For binary clauses, add implications.
    // Attention! first check subsumptions, then add implications!
    // Otherwise we might remove the implication in subsumption check.
    if (clause.size() == 2) {
        // binary clauses are additionally stored as implications.
        addImplication(negate(clause[0]), clause[1], c_nr);
        addImplication(negate(clause[1]), clause[0], c_nr);
    }

    // mark clause as modified
    _gates.resizeClauses(c_nr + 1);
    _gates.touchClause(c_nr);

    val_assert(_clauses[c_nr].checkConsistency());
    return static_cast<int>(c_nr);
}

/**
 * \brief Adds a clause to the formula.
 *
 * \note Calling this function is only allowed as long as none of
 *       the functions has been called which modify the formula
 *       (like unitPropagation(), findPure(), ...)
 * \param clause vector with the literals of the clause
 * \return the clause ID if the clause was actually added; -2 if the clause was
 * a tautology; -1 if the clause was unit. \throw UNSATException if an empty
 * clause is added.
 */
int
Formula::addClause(const Clause& clause)
{
    if (clause.empty()) {
        throw UNSATException("Trying to add an empty clause");
    }

    if (clause.isTautology()) {
        return -2;
    }

    if (clause.size() == 1) {
        pushUnit(clause[0], PureStatus::UNIT);
        return -1;
    }

    ClauseID c_nr = 0;
    val_assert(_clauses.size() == _clause_sizes.size());

    // Try to recycle clause numbers where possible.
    if (_deleted_clause_numbers.empty()) {
        _clause_sizes.push_back(clause.size());
        _clauses.push_back(clause);
        c_nr = _clauses.size() - 1;
    } else {
        c_nr = _deleted_clause_numbers.back();
        _deleted_clause_numbers.pop_back();
        _clause_sizes[c_nr] = clause.size();
        _clauses[c_nr]      = clause;
    }

    return addClauseToLists(c_nr);
}

/**
 * \brief Adds a clause to the formula.
 *
 * \note Calling this function is only allowed as long as none of
 *       the functions has been called which modify the formula
 *       (like unitPropagation(), findPure(), ...)
 * \param clause vector with the literals of the clause
 * \return the clause ID if the clause was actually added; -2 if the clause was
 * a tautology; -1 if the clause was unit. \throw UNSATException if an empty
 * clause is added.
 */
int
Formula::addClause(Clause&& clause)
{
    if (clause.empty()) {
        throw UNSATException("Trying to add an empty clause");
    }

    if (clause.isTautology()) {
        return -2;
    }

    if (clause.size() == 1) {
        pushUnit(clause[0], PureStatus::UNIT);
        return -1;
    }

    ClauseID c_nr = 0;
    val_assert(_clauses.size() == _clause_sizes.size());

    // Try to recycle clause numbers where possible.
    if (_deleted_clause_numbers.empty()) {
        _clause_sizes.push_back(clause.size());
        _clauses.push_back(std::move(clause));
        c_nr = _clauses.size() - 1;
    } else {
        c_nr = _deleted_clause_numbers.back();
        _deleted_clause_numbers.pop_back();
        _clause_sizes[c_nr] = clause.size();
        _clauses[c_nr]      = std::move(clause);
    }

    return addClauseToLists(c_nr);
}

/**
 * \brief Adds a clause to the formula.
 *
 * \note Calling this function is only allowed as long as none of
 *       the functions has been called which modify the formula
 *       (like unitPropagation(), findPure(), ...)
 * \param clause vector with the literals of the clause
 * \param needs_sorting if true, the clause is sorted and duplicate literals are
 * removed. \pre If needs_sorting is false, the caller has to guarantee that the
 * literals are ordered in increasing order and that no literal appears more
 * than once. \param check_subsumption if true, it is checked whether the clause
 * subsumes other clauses in database \param status Specifies if the clause is
 * mandatory or optional (deleted does not make much sense here!) \return the
 * clause ID if the clause was actually added; -1 if the clause was a tautology.
 * \throw UNSATException if an empty clause is added.
 */
int
Formula::addClause(Clause::ClauseData&& clause, bool needs_sorting, bool check_subsumption, ClauseStatus status)
{
    if (clause.empty()) {
        throw UNSATException("Trying to add an empty clause");
    }

    if (clause.size() == 1) {
        pushUnit(clause[0], PureStatus::UNIT);
        return -1;
    }

    ClauseID c_nr = 0;
    val_assert(_clauses.size() == _clause_sizes.size());

    // Try to recycle clause numbers where possible.
    if (_deleted_clause_numbers.empty()) {
        _clause_sizes.push_back(clause.size());
        _clauses.emplace_back(std::move(clause), needs_sorting, status);
        c_nr = _clauses.size() - 1;
    } else {
        c_nr = _deleted_clause_numbers.back();
        _deleted_clause_numbers.pop_back();
        _clause_sizes[c_nr] = clause.size();
        _clauses[c_nr]      = Clause(std::move(clause), needs_sorting, status);
    }

    return addClauseToLists(c_nr, check_subsumption);
}

/**
 * \brief Removes a clause from the formula.
 *
 * \param c_nr the ID of the clause
 * \note This function modifies Formula::_occ_list and Formula::_implications.
 *       Therefore Formula::removeClause may not be called while iterating over
 *       one of these lists.
 * \note The removed clause may not be accessed anymore.
 */
void
Formula::removeClause(const ClauseID c_nr)
{
    val_assert(c_nr < _clauses.size());
    if (clauseDeleted(c_nr)) {
        val_assert(_clauses[c_nr].empty());
        return;
    }

    for (const Literal lit : _clauses[c_nr]) {
        // remove c_nr from _occ_list[lit]:
        removeFromOccList(lit, c_nr);
    }

    if (_clauses[c_nr].size() == 2) {
        // remove both implications
        removeImplication(negate(_clauses[c_nr][0]), _clauses[c_nr][1]);
        removeImplication(negate(_clauses[c_nr][1]), _clauses[c_nr][0]);
    }

    // clear the clause and free its number
    _clauses[c_nr].clear();
    _clauses[c_nr].setStatus(ClauseStatus::DELETED);
    _clause_sizes[c_nr] = 0;

    _gates.touchClause(c_nr);
    _deleted_clause_numbers.push_back(c_nr);
}

/**
 * \brief Removes all clauses which are marked as optional.
 * \return true if the formula was modified
 * \sa Formula::removeClause(ClauseID)
 */
bool
Formula::removeOptionalClauses()
{
    std::size_t count = 0;

    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (clauseOptional(c_nr)) {
            removeClause(c_nr);
            ++count;
        }
    }

    return count > 0;
}

/**
 * \brief Removes the given literal from a clause.
 * \param c_nr the ID of the clause from this the literal is to be deleted
 * \param lit the literal to be deleted
 * \pre The literal has to be contained in the clause. Otherwise the
 *      result is undefined.
 * \return Returns true if clause was further reduced by universal reduction
 */
bool
Formula::removeLiteral(const ClauseID c_nr, const Literal lit)
{
    val_assert(!clauseDeleted(c_nr));
    val_assert(std::find(_clauses[c_nr].cbegin(), _clauses[c_nr].cend(), lit) != _clauses[c_nr].cend());
    auto& clause = _clauses[c_nr];
    bool  result = false;

    // If the clause is currently binary, it becomes unit.
    if (clause.size() == 2) {
        val_assert(clause[0] == lit || clause[1] == lit);
        if (clause[0] == lit) {
            pushUnit(clause[1], PureStatus::UNIT);
        } else {
            pushUnit(clause[0], PureStatus::UNIT);
        }
        removeClause(c_nr);
        return true;
    }

    // Remove the clause from the literal's occurrence list.
    removeFromOccList(lit, c_nr);

    // Remove the literal from the clause.
    clause.erase(std::remove(clause.begin(), clause.end(), lit), clause.end());
    clause.computeSignature();

    // Perform universal reduction
    if (_settings.univ_reduction) {
        result = universalReduction(clause, static_cast<int>(c_nr));
    }

    // If the clause is now binary, add the implications.
    if (clause.size() == 2) {
        addImplications(clause[0], clause[1], c_nr);
    }

    _gates.touchClause(c_nr);

    return result;
}

/**
 * \brief Replace, in the whole formula, a literal by another one.
 *
 * The negated literal is replaced by the negated replacement, i.e.,
 * afterwards the formula does not contain var(lit) anymore.
 * \pre var(lit) != var(replacement)
 * \param lit the literal that is to be replaced
 * \param replacement the replacement of the literal
 */

void
Formula::replaceLiteral(const Literal lit, const Literal replacement)
{
    replaceLiteralMono(lit, replacement);
    replaceLiteralMono(negate(lit), negate(replacement));
    removeVar(lit2var(lit));
    unitPropagation();
}

/**
 * \brief Replace a literal by another one in all clauses.
 *
 * \note The negated literal is not replaced!
 * \pre \f$var(lit) \neq var(replacement)\f$
 * \pre \f$var(replacement) \not\in clause\f$
 * \param lit the literal that is to be replaced
 * \param replacement the replacement literal of the literal `lit`
 */
void
Formula::replaceLiteralMono(const Literal lit, const Literal replacement)
{
    val_assert(!varDeleted(lit2var(lit)));
    val_assert(lit2var(lit) != lit2var(replacement));

    const Literal neg_lit         = negate(lit);
    const Literal neg_replacement = negate(replacement);

    // Remove all implications of ~lit
    for (const auto impl : _implications[neg_lit]) {
        removeImplication(negate(impl.getLiteral()), lit);
    }
    _implications[neg_lit].clear();

    while (!_occ_list[lit].empty()) {
        const ClauseID c_nr = _occ_list[lit].back();
        _gates.touchClause(c_nr);
        _occ_list[lit].pop_back();
        auto& clause = _clauses[c_nr];
        bool  skip   = false;

        if (lit < replacement) {
            unsigned int l_pos = 0;
            // find the position of 'lit'
            while (clause[l_pos] != lit) {
                ++l_pos;
            }
            // copy all literals between lit and replacement one position to the left
            while (l_pos + 1 < clause.size() && clause[l_pos + 1] <= replacement) {
                clause[l_pos] = clause[l_pos + 1];
                ++l_pos;
            }
            if (clause[l_pos] == replacement) {
                // The clause already contains the replacement
                // -> move all literals right of l_pos one position left and
                //    shorten the clause by one
                while (l_pos + 1 < clause.size()) {
                    clause[l_pos] = clause[l_pos + 1];
                    ++l_pos;
                }

                // shrink the clause length by 1
                clause.resize(clause.size() - 1);

                // if the clause has become binary, add to implication lists
                if (clause.size() == 2) {
                    addImplications(clause[0], clause[1], c_nr);
                }

                // if the clause has become unit, propagate it.
                if (clause.size() == 1) {
                    pushUnit(clause[0], PureStatus::UNIT);
                    skip = true;
                }
            } else {
                clause[l_pos] = replacement;

                // check for tautology
                if (l_pos > 0 && clause[l_pos - 1] == neg_replacement) {
                    skip = true;
                } else if (l_pos < clause.size() - 1 && clause[l_pos + 1] == neg_replacement) {
                    skip = true;
                }

                if (!skip) {
                    _occ_list[replacement].push_back(c_nr);
                    if (clause.size() == 2) {
                        addImplications(clause[0], clause[1], c_nr);
                    }
                }
            }
        } else {
            // lit > replacement
            auto l_pos = static_cast<int>(clause.size() - 1);
            while (clause[l_pos] > lit) {
                --l_pos;
            }
            while (l_pos > 0 && clause[l_pos - 1] >= replacement) {
                clause[l_pos] = clause[l_pos - 1];
                --l_pos;
            }

            if (clause[l_pos] == replacement) {
                // The clause already contained the replacement -> shorten it
                while (l_pos < static_cast<int>(clause.size()) - 1) {
                    clause[l_pos] = clause[l_pos + 1];
                    ++l_pos;
                }

                // shrink the clause length by 1
                clause.resize(clause.size() - 1);

                // if the clause has become binary, add to implication lists
                if (clause.size() == 2) {
                    addImplications(clause[0], clause[1], c_nr);
                }

                // if the clause has become unit, propagate it.
                if (clause.size() == 1) {
                    pushUnit(clause[0], PureStatus::UNIT);
                    skip = true;
                }
            } else {
                clause[l_pos] = replacement;

                // check for tautology
                if (l_pos > 0 && clause[l_pos - 1] == neg_replacement) {
                    skip = true;
                } else if (l_pos < static_cast<int>(clause.size()) - 1 && clause[l_pos + 1] == neg_replacement) {
                    skip = true;
                }

                if (!skip) {
                    _occ_list[replacement].push_back(c_nr);
                    if (clause.size() == 2) {
                        addImplications(clause[0], clause[1], c_nr);
                    }
                }
            }
        }

        if (skip) {
            removeClause(c_nr);
            _clause_sizes[c_nr] = 0;
        } else {
            if (_settings.univ_reduction) {
                universalReduction(clause, static_cast<int>(c_nr));
            }
            clause.computeSignature();
            _clause_sizes[c_nr] = clause.size();
            // clause can be removed/identified as unit in universalReduction()
            if (!clause.empty()) {
                const auto num_subsumed = isBackwardSubsuming(clause, static_cast<int>(c_nr));
                VLOG_IF(num_subsumed > 0, 3) << __FUNCTION__ << " removed " << num_subsumed << " subsumed clauses.";
            }
        }
    }

    val_assert(_occ_list[lit].empty());
    val_assert(_implications[negate(lit)].empty());
}

/**
 * \brief Resets the formula and clears all the data.
 *
 * Afterwards, the list of variables and clauses as well as
 * all auxiliary data structures (implications, occurrence lists)
 * are empty. The statistics variables are reset to 0.
 */
void
Formula::reset()
{
    val_assert(_prefix);
    _prefix->clear();
    _clauses.clear();
    _gates.clear();
    _occ_list.clear();
    _implications.clear();
    _deleted_clause_numbers.clear();
    _deleted_var_numbers.clear();
    _unit_stack.clear();
    _assignment.clear();
    _dont_touch.clear();
    _variable_score.clear();
    _clause_sizes.clear();
    _candidates.clear();
    _blocked_candidates.clear();
    _seen.clear();
    _seen2.clear();
#ifdef SKOLEM
    _skolem_data.clear();
#endif

    for (auto& timer : _timers) {
        timer.reset();
    }
    std::fill(_statistics.begin(), _statistics.end(), 0);
}

/**
 * \brief Creates a copy of a given variable.
 *
 * The copy has the same dependencies as the copied variable.
 * \return the index of the copy
 */
Variable
Formula::copyVar(const Variable var)
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    val_assert(!varDeleted(var));

    const Variable new_var = nextVar();

    val_assert(minVarIndex() <= new_var && new_var <= maxVarIndex());
    val_assert(varDeleted(new_var));

    _prefix->makeCopy(new_var, var);
    return new_var;
}

/**
 * \brief Checks if start_lit implies target_var or NOT(target_var) taking
 * transitivity into account.
 *
 * \return a pair (a,b) such that a = true means start_lit implies target_var
 * and b = true that start_lit implies NOT(target_var).
 */
std::pair<bool, bool>
Formula::hasImplicationTransitive(const Literal start_lit, const Variable target_var) const
{
    val_assert(start_lit < _implications.size());

    const auto reached = impliedLiterals(start_lit);
    return std::make_pair(reached[var2lit(target_var, false)], reached[var2lit(target_var, true)]);
}

std::vector<bool>
Formula::impliedLiterals(const Literal start_lit) const
{
    std::stack<Literal> front;
    front.push(start_lit);

    std::vector<bool> reached(_implications.size(), false);
    reached[start_lit] = true;

    while (!front.empty()) {
        const Literal top = front.top();
        front.pop();
        for (const BinaryClause succ : _implications[top]) {
            const Literal lit = succ.getLiteral();
            if (!reached[lit]) {
                reached[lit] = true;
                front.push(lit);
            }
        }
    }

    return reached;
}

/**
 * \brief Checks if the DQBF is actually a QBF with a linearly ordered
 * quantifier prefix
 *
 * This is done by testing whether for the dependency sets \f$D_y\f$ and
 * \f$D_{y'}\f$ of all pairs \f$y,y'\f$ of existential variables either
 * \f$D_y\subseteq D_{y'}\f$ or \f$D_{y'}\subseteq D_y\f$ holds.
 *
 * \todo Check if the signature approach (cf. subsumption checks) makes this
 *       procedure more efficient. A signature can be defined by hashing the
 *       the variables in the dependency set of an existential variable to
 *       the range 0...63 and setting the corresponding bit position in the
 *       signature to 1. If sig(prefix[i]) & ~sig(prefix[j]) != 0, then
 * prefix[i] cannot be a subset of prefix[j]. However, since the dependency sets
 * are typically larger than an average clause, this might not help much.
 */
bool
Formula::isQBF() const noexcept
{
    val_assert(_prefix);
    val_assert(_prefix->type() == PrefixType::QBF || _dqbf_prefix);

    if (_prefix->type() == PrefixType::QBF && _qbf_prefix != nullptr && _dqbf_prefix == nullptr) {
        return true;
    }

    return _dqbf_prefix->isQBF();
}

/**
 * \brief If the current formula has an equivalent QBF prefix, this function
 * constructs such a QBF prefix. \pre This function requires that the DQBF has
 * an equivalent QBF prefix.
 */
void
Formula::convertToQBF()
{
    val_assert(_prefix);
    val_assert(isQBF());

    if (_prefix->type() == PrefixType::QBF && _qbf_prefix != nullptr) {
        return;
    }

    QBFPrefix* new_prefix = _dqbf_prefix->convertToQBF();
    delete _prefix;
    _qbf_prefix  = new_prefix;
    _prefix      = new_prefix;
    _dqbf_prefix = nullptr;
}

void
Formula::initCandidateLists()
{
    _blocked_candidates.clear();
    _blocked_candidates.resize(_clauses.size());
    for (unsigned c_nr = 0; c_nr != _clauses.size(); ++c_nr) {
        if (clauseDeleted(c_nr)) {
            continue;
        }
        if (!_blocked_candidates.inHeap(c_nr)) {
            _blocked_candidates.insert(c_nr);
        }
    }
}

/**
 * \brief Counts occurences of all literals and stores it into `_variable_score`
 */
void
Formula::countOccurences()
{
    _variable_score.resize(maxLitIndex() + 1, 0);

    for (Literal l = minLitIndex(); l <= maxLitIndex(); ++l) {
        _variable_score[l] = static_cast<int>(_occ_list[l].size());
    }
}

/**
 * \brief Counts occurences of all variables and stores it into
 * `_variable_score`
 */
void
Formula::countOccurencesVar()
{
    _variable_score.resize(maxVarIndex() + 1, 0);

    for (Variable v = minVarIndex(); v <= maxVarIndex(); ++v) {
        _variable_score[v] = static_cast<int>(_occ_list[var2lit(v)].size() + _occ_list[negate(var2lit(v))].size());
    }
}

/**
 * \brief Counts implications of all variables and stores it into
 * `_variable_score`
 */
void
Formula::countImplications()
{
    val_assert(_variable_score.size() > maxLitIndex());

    for (Literal l = minLitIndex(); l <= maxLitIndex(); ++l) {
        _variable_score[l] = static_cast<int>(_implications[l].size());
    }
}

/**
 * \brief Prints the current configuration of the preprocessor.
 */
void
Formula::printSettings() const
{
    LOG(INFO) << _settings << '\n';
}

/*
void Formula::setVerbosity(const unsigned short val) noexcept
{
    _settings.verbosity = val;
    el::Loggers::setVerboseLevel(val);
}
*/

/**
 * \brief Prints some statistics about the effectiveness of the preprocessor.
 */
void
Formula::printStatistics() const
{
    LOG(INFO) << "STATISTICS ON THE PREPROCESSING PROCESS:\n"
              << "c num_unit .....................= " << getStatistics(Statistics::UNIT) << '\n'
              << "c num_pure .....................= " << getStatistics(Statistics::PURE) << '\n'
              << "c num_univ_reduction............= " << getStatistics(Statistics::UNIV_REDUCTION) << '\n'
              << "c num_univ_expansion............= " << getStatistics(Statistics::UNIV_EXPANSION) << '\n'
              << "c num_equiv.....................= " << getStatistics(Statistics::EQUIV_LITS) << '\n'
              << "c num_hidden_equiv..............= " << getStatistics(Statistics::HIDDEN_EQUIV_LITS) << '\n'
              << "c num_equiv_gates...............= " << getStatistics(Statistics::EQUIV_GATES) << '\n'
              << "c num_gates.....................= "
              << getStatistics(Statistics::GATES_AND) + getStatistics(Statistics::GATES_XOR)
                     + getStatistics(Statistics::GATES_MUX) + getStatistics(Statistics::GATES_SEMANTIC)
              << '\n'
              << "c   AND gates...................= " << getStatistics(Statistics::GATES_AND) << '\n'
              << "c   XOR gates...................= " << getStatistics(Statistics::GATES_XOR) << '\n'
              << "c   MUX gates...................= " << getStatistics(Statistics::GATES_MUX) << '\n'
              << "c   semantic gates..............= " << getStatistics(Statistics::GATES_SEMANTIC) << '\n'
              << "c num_blocked_clauses...........= " << getStatistics(Statistics::BCE) << '\n'
              << "c num_blocked_univ_literals.....= " << getStatistics(Statistics::BLE) << '\n'
              << "c num_added_blocked_literals....= " << getStatistics(Statistics::BLA) << '\n'
              << "c num_added_blocked_implications= " << getStatistics(Statistics::BIA) << '\n'
              << "c num_hidden_tautologies........= " << getStatistics(Statistics::HTE) << '\n'
              << "c num_hidden_subsumptions.......= " << getStatistics(Statistics::HSE) << '\n'
              << "c num_impl_chains...............= " << getStatistics(Statistics::IMPLICATION_CHAINS) << '\n'
              << "c num_contradictions............= " << getStatistics(Statistics::CONTRADICTIONS) << '\n'
              << "c num_hidden_equiv..............= " << getStatistics(Statistics::HIDDEN_EQUIV_LITS) << '\n'
              << "c num_hidden_unit...............= " << getStatistics(Statistics::HIDDEN_UNIT) << '\n'
              << "c num_substitution..............= " << getStatistics(Statistics::SUBSTITUTION) << '\n'
              << "c num_rewritings................= " << getStatistics(Statistics::REWRITING) << '\n'
              << "c num_subsumption...............= " << getStatistics(Statistics::SUBSUMPTION) << '\n'
              << "c num_selfsubsuming_resolution..= " << getStatistics(Statistics::SELF_SUBSUMPTION) << '\n'
              << "c num_resolution................= " << getStatistics(Statistics::RESOLUTION) << '\n'
              << "c num_unit_by_sat...............= " << getStatistics(Statistics::UNIT_SAT) << '\n'
              << "c num_impl_by_sat...............= " << getStatistics(Statistics::IMPL_SAT) << '\n'
              << "c num_sat_calls.................= " << getStatistics(Statistics::SAT_CALLS) << '\n'
              << "c num_vivified_literals.........= " << getStatistics(Statistics::VIVIFIED_LITERALS) << '\n'
              << "c num_vivified_clauses..........= " << getStatistics(Statistics::VIVIFIED_CLAUSES) << '\n'
              << "c num_add_dep_schemes...........= " << getStatistics(Statistics::ADD_DEPENDENCY_SCHEMES) << '\n'
              << "c num_rem_dep_schemes...........= " << getStatistics(Statistics::REM_DEPENDENCY_SCHEMES) << '\n'
              << "c num_preprocessor_loops........= " << getStatistics(Statistics::PREPRO_LOOPS) << '\n'
              << std::setprecision(2) << std::fixed
              << "c time_univ_reduction...........= " << getTime(WhichTimer::UNIV_REDUCTION) << "s\n"
              << "c time_univ_expansion...........= " << getTime(WhichTimer::UNIV_EXPANSION) << "s\n"
              << "c time_equiv....................= " << getTime(WhichTimer::EQUIV_LITS) << "s\n"
              << "c time_equiv_gates..............= " << getTime(WhichTimer::EQUIV_GATES) << "s\n"
              << "c time_hidden_equiv_contra......= " << getTime(WhichTimer::HIDDEN_EQUIV_CONTRA) << "s\n"
              << "c time_blocked_clauses..........= " << getTime(WhichTimer::BCE) << "s\n"
              << "c time_hidden_subsumptions......= " << getTime(WhichTimer::HSE) << "s\n"
              << "c time_impl_chains..............= " << getTime(WhichTimer::IMPLICATION_CHAINS) << "s\n"
              << "c time_contradictions...........= " << getTime(WhichTimer::CONTRADICTIONS) << "s\n"
              << "c time_substitution.............= " << getTime(WhichTimer::SUBSTITUTION) << "s\n"
              << "c time_rewriting................= " << getTime(WhichTimer::REWRITING) << "s\n"
              << "c time_subsumption..............= " << getTime(WhichTimer::SUBSUMPTION) << "s\n"
              << "c time_selfsubsuming_resolution.= " << getTime(WhichTimer::SELF_SUBSUMPTION) << "s\n"
              << "c time_resolution...............= " << getTime(WhichTimer::RESOLUTION) << "s\n"
              << "c time_unit_by_sat..............= " << getTime(WhichTimer::UNIT_SAT) << "s\n"
              << "c time_impl_by_sat..............= " << getTime(WhichTimer::IMPLICATIONS_SAT) << "s\n"
              << "c time_incomplete_sat_checks....= " << getTime(WhichTimer::INCOMPLETE_SAT) << "s\n"
              << "c time_dep_schemes..............= " << getTime(WhichTimer::DEPENDENCY_SCHEMES) << "s\n"
              << "c time_upla.....................= " << getTime(WhichTimer::UPLA) << "s\n"
              << "c time_vivification.............= " << getTime(WhichTimer::VIVIFICATION) << "s\n"
              << "c time_gate_detection...........= " << getTime(WhichTimer::GATE_DETECTION) << "s\n"
              << "c   time for AND detection......= " << getTime(WhichTimer::GATE_AND_DETECTION) << "s\n"
              << "c   time for XOR detection......= " << getTime(WhichTimer::GATE_XOR_DETECTION) << "s\n"
              << "c   time for MUX detection......= " << getTime(WhichTimer::GATE_MUX_DETECTION) << "s\n"
              << "c   time for semantic detection.= " << getTime(WhichTimer::GATE_SEMANTIC_DETECTION) << "s\n"
              << "c time_preprocessing............= " << getTime(WhichTimer::TOTAL_TIME) << "s\n"
              << std::flush;
}



Formula& Formula::operator=(Formula& other)
{
    _clauses = other._clauses;
    _gates = other._gates;
    _occ_list = other._occ_list;
    _implications = other._implications;
    _implications_added = other._implications_added;
    _deleted_clause_numbers = other._deleted_clause_numbers;
    _deleted_var_numbers = other._deleted_var_numbers;
    _unit_stack = other._unit_stack;
    _assignment = other._assignment;
    _dont_touch = other._dont_touch;
    _variable_score = other._variable_score;
    _clause_sizes = other._clause_sizes;
    _candidates = other._candidates;
    _blocked_candidates = other._blocked_candidates;
    _removed_lits = other._removed_lits;
    _seen = other._seen;
    _seen2 = other._seen2;
#ifdef SKOLEM
    _skolem_data = other._skolem_data;
#endif
    _settings = other._settings;
    _process_limit = other._process_limit;
    _interrupt = other._interrupt;
    _statistics = other._statistics;
    _timers = other._timers;

    if (other._dqbf_prefix) {
        _dqbf_prefix = other._dqbf_prefix;
        _prefix = _dqbf_prefix;
    } else if (other._qbf_prefix) {
        _qbf_prefix = other._qbf_prefix;
        _prefix = _qbf_prefix;
    }
	
	other._prefix = nullptr;
    other._dqbf_prefix = nullptr;
    other._qbf_prefix = nullptr;
	
	return *this;

}



}  // end namespace hqspre
