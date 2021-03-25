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

#include <cstddef>
#include <cstdint>
#include <set>
#include <vector>

#include <easylogging++.hpp>

#include "aux.hpp"
#include "clause.hpp"
#include "formula.hpp"
#include "gate.hpp"
#include "literal.hpp"
#include "prefix.hpp"
#include "process_limits.hpp"
#include "timer.hpp"
#include "varheap.hpp"

/**
 * \file formula_blocked_et_al.cpp
 * \brief Implementation of hidden tautology elimination, hidden subsumption
 * elimination \author Sven Reimer \date 06/2016
 */

namespace hqspre {

/**
 * \brief Checks if the resolvent of two clauses is a tautology.
 *
 * The resolvent is only a valid tautology if the dependency set of the
 * variable that makes the resolvent a tautology is a subset of the pivot var's
 * dependency set. This is the requirement, e.g., in blocked clause detection.
 *
 * \pre The literals of the other clause must be marked in the "_seen"
 * datastructure \param clause the clause containing the pivot variable \param
 * pivot_var the pivot variable \return true if the resolvent of c1 and c2
 * w.r.t. pivot_var is a valid tautology
 */
template <typename Container>
bool
Formula::checkResolventTautology(const Container& clause, const Variable pivot_var) const
{
    for (const Literal lit : clause) {
        // skip pivot var
        if (lit2var(lit) == pivot_var) {
            continue;
        }
        // found tautology
        if (_seen[negate(lit)]) {
            if (dependenciesSubset(lit2var(lit), pivot_var)) {
                return true;
            }
        }
    }
    return false;
}

/* \brief Tries to add blocking universal literals.
 * \param clause a copy of the clause as an ordered container into which the
 * blocked literals are inserted. \note The clause will not be extended; it is
 * only checked whether universal literals can be added which are blocked
 * literals \return true if clause can be extended by universal blocked literals
 */
template <typename Container>
bool
Formula::addBlockingLiterals(const Container& clause)
{
    // Blocking literals are only sound for QBF, but not for DQBF.
    if (!isQBF()) {
        return false;
    }

    for (const Literal lit : clause) {
        if (isExistential(lit2var(lit))) {
            continue;
        }

        if (clauseBlockedByLit(negate(lit))) {
            ++stat(Statistics::BLA);
            clearSeen(clause);
            return true;
        }
    }
    return false;
}

/**
 * \brief Tries to add further implications to the formula which are blocked.
 * \return true if the formula was modified.
 */
bool
Formula::addBlockedImplications()
{
    VLOG(1) << __FUNCTION__;

    std::size_t        count = 0;
    Clause::ClauseData bin_clause(2, 0);

    for (Literal lit1 = minLitIndex(); lit1 < maxLitIndex(); ++lit1) {
        if (_interrupt) {
            break;
        }
        const Variable var1 = lit2var(lit1);
        if (varDeleted(var1)) {
            continue;
        }
        for (Literal lit2 = lit1 + 1; lit2 <= maxLitIndex(); ++lit2) {
            if (_interrupt) {
                break;
            }
            const Variable var2 = lit2var(lit2);
            if (varDeleted(var2) || var1 == var2) {
                continue;
            }
            if (hasImplication(negate(lit1), lit2) != -1) {
                continue;
            }

            bin_clause[0] = lit1;
            bin_clause[1] = lit2;
            _seen[lit1]   = 1u;
            _seen[lit2]   = 1u;
            //            addHiddenLiteralsBinary(-1, bin_clause, sign);
            if (clauseBlocked(bin_clause) != 0u) {
                //                bin_clause.resize(2);
                //                bin_clause[0] = lit1;
                //                bin_clause[1] = lit2;
                addClause(bin_clause);  //, false, true, ClauseStatus::OPTIONAL);
                ++count;
                VLOG(3) << "Adding blocked binary clause (" << lit2dimacs(lit1) << ", " << lit2dimacs(lit2) << ')';
            }
            clearSeen();
        }
    }

    stat(Statistics::BIA) += count;
    VLOG(2) << __FUNCTION__ << "(): added " << count << " blocked binary clauses.";
    return count > 0;
}

/**
 * \brief Removes the clause "c_nr" due to hidden tautology or blocked clause
 * elimination and updates internal candidate lists
 */
void
Formula::removeClauseAndUpdateCandidates(const ClauseID c_nr)
{
    checkPure(_clauses[c_nr], 0);

    // Recreate original clause
    // Since it will be deleted in the next step and the upcoming
    // Implication chain check should be applied without this clause
    const auto org_clause = _clauses[c_nr].getLiterals();

    // Delete clause and its occurences
    removeClause(c_nr);

    unitPropagation();

    // No update of candidates needed if implication chains and blocked clauses
    // are not active
    if (!_settings.impl_chains && _settings.bce == 0u) {
        return;
    }

    // Check whether new implication chains are produced
    // And update candidates
    for (const Literal lit : org_clause) {
        unsigned int check_lit = lit;
        if (_settings.impl_chains) {
            check_lit = checkImplicationChain(lit);
            if (check_lit != 0) {
                if (_settings.bce > 0) {
                    // add new candidates if bce is applied
                    for (const ClauseID occ_nr : _occ_list[check_lit]) {
                        if (clauseDeleted(occ_nr)) {
                            continue;
                        }
                        if (!_blocked_candidates.inHeap(occ_nr)) {
                            _blocked_candidates.insert(occ_nr);
                        }
                    }
                }
            } else {
                check_lit = lit;
            }
        }

        // Update candidates
        if (_settings.bce > 0) {
            // add new candidates if bce is applied
            for (const ClauseID occ_nr : _occ_list[negate(check_lit)]) {
                val_assert(!clauseDeleted(occ_nr));
                if (!_blocked_candidates.inHeap(occ_nr)) {
                    _blocked_candidates.insert(occ_nr);
                }
            }
        }
    }  // end for (lit)
}

/**
 * \brief Finds and removes hidden tautologies, blocked clauses, subsumptions
 * from the formula.
 *
 * \return true iff the formula was modified.
 */
bool
Formula::removeBlockedAndFriends()
{
    VLOG(1) << __FUNCTION__;

    bool       modified = false;
    ScopeTimer bce(getTimer(WhichTimer::BCE));
    _process_limit.setLimit(PreproMethod::BCE);

    const std::size_t old_stat_hte         = stat(Statistics::HTE);
    const std::size_t old_stat_bce         = stat(Statistics::BCE);
    const std::size_t old_stat_hse         = stat(Statistics::HSE);
    const std::size_t old_stat_bla         = stat(Statistics::BLA);
    const std::size_t old_stat_impl_chains = stat(Statistics::IMPLICATION_CHAINS);

    val_assert(_seen.size() > maxLitIndex());
    Clause::ClauseData clause;

    initCandidateLists();

    while (!_blocked_candidates.empty()) {
        const ClauseID c_nr = _blocked_candidates.top();

        if (_interrupt) {
            break;
        }
        if (clauseDeleted(c_nr)) {
            continue;
        }
        if (_settings.preserve_gates && _gates.isGateClause(c_nr)) {
            continue;
        }
        if (_clauses[c_nr].size() < _settings.min_bce_size) {
            continue;
        }

        if (_process_limit.reachedLimit()) {
            VLOG(3) << "terminate " << __FUNCTION__ << " due to process limit.\n";
            break;
        }

        const Clause& c    = _clauses[c_nr];
        std::uint64_t sign = 0;

        clause.resize(c.size());

        for (std::size_t j = 0; j != c.size(); ++j) {
            clause[j] = c[j];
            val_assert(!_seen[c[j]]);
            _seen[c[j]] = 1u;
        }

        bool hidden_tautology = false;

        if (_settings.hidden != 0u || _settings.covered || _settings.bla) {
            // use incomplete hidden literal addition
            if (_settings.hidden == 1) {
                hidden_tautology = addHiddenLiteralsBinary(static_cast<int>(c_nr), clause, sign);
            }
            // use complete hidden literal addition
            else if (_settings.hidden == 2) {
                hidden_tautology = addHiddenLiterals(static_cast<int>(c_nr), clause, sign);
            }

            // add covered literals (if clause in not already a hidden tautology)
            if (!hidden_tautology && _settings.covered) {
                hidden_tautology = addCoveredLiterals(clause, sign);
            }

            // adds blocked universal literals
            if (!hidden_tautology && _settings.bla && _prefix->type() == PrefixType::QBF) {
                hidden_tautology = addBlockingLiterals(clause);
            }
        }

        if (hidden_tautology) {
            VLOG(3) << __FUNCTION__ << "(): clause " << c << " is a hidden tautology.";
            ++stat(Statistics::HTE);

            removeClauseAndUpdateCandidates(c_nr);
            modified = true;
        }
        // check for blocked clauses
        else {
            bool clause_blocked = false;

            if (_settings.bce > 0) {
                if (clauseBlocked(clause) != 0u) {
                    VLOG(3) << __FUNCTION__ << "(): clause " << _clauses[c_nr] << " is (hidden) blocked.";
                    ++stat(Statistics::BCE);

                    removeClauseAndUpdateCandidates(c_nr);
                    modified       = true;
                    clause_blocked = true;
                }
            }

            // check for hidden subsumptions if clause is not a tautology or blocked
            // also skip if we applied complete hidden literal addition.
            // In this case the hidden subsumption check was done implicitely during
            // the hidden literal addtion
            if (_settings.hidden != 2 && _settings.hse && !clause_blocked) {
                getTimer(WhichTimer::BCE).stop();
                getTimer(WhichTimer::HSE).start();

                // Check of subsumption, excluded current clause id "i"
                if (isForwardSubsumed(clause, sign, static_cast<int>(c_nr))) {
                    VLOG(3) << __FUNCTION__ << "(): clause " << c << " is (hidden) subsumed.";
                    ++stat(Statistics::HSE);
                    removeClauseAndUpdateCandidates(c_nr);
                    modified = true;
                }
                getTimer(WhichTimer::HSE).stop();
                getTimer(WhichTimer::BCE).start();
            }
        }
        clearSeen(clause);
    }

    VLOG(2) << __FUNCTION__ << " found " << (stat(Statistics::HTE) - old_stat_hte) << " hidden tautologies, "
            << (stat(Statistics::BCE) - old_stat_bce) << " hidden blocked clauses with "
            << (stat(Statistics::BLA) - old_stat_bla) << " added blocked literals and "
            << (stat(Statistics::HSE) - old_stat_hse) << " hidden subsumptions. Additionally "
            << (stat(Statistics::IMPLICATION_CHAINS) - old_stat_impl_chains) << " implication chains are found.";

    if (numClauses() == 0) {
        throw SATException("No clauses left after BCE");
    }
    return modified;
}

/**
* \brief Adds hidden literals obtained by only checking binary clauses

* Also checks if the given clause is a hidden tautology.
* \return true if `clause` is a hidden tautology
* \pre Assumes that "_seen" vector was initialized with "clause"
* \param clause the clause to check
* \param sign the signature of `clause`
* \param c_nr the clause ID of `clause`
*/
bool
Formula::addHiddenLiteralsBinary(const int c_nr, Clause::ClauseData& clause, std::uint64_t& sign) const
{
    std::size_t csize = clause.size();

    for (std::size_t i = 0; i != csize; ++i) {
        const Literal l0 = clause[i];
        addSignatureLit(sign, l0);

        _process_limit.decreaseLimitBy(_implications[negate(l0)].size());

        for (const BinaryClause bin : _implications[negate(l0)]) {
            // Skip currently considered binary clause
            if (static_cast<int>(bin.getClauseID()) == c_nr) {
                continue;
            }

            const Literal impl = bin.getLiteral();
            // We found a hidden tautology!
            if (_seen[impl] == 1u) {
                return true;
            } else if (_seen[negate(impl)] == 0u) {
                // otherwise add hidden literal
                clause.push_back(negate(impl));
                ++csize;
                _seen[negate(impl)] = 1u;
            }
        }
    }

    return false;
}

/**
 * \brief Adds hidden literals to a given clause.
 * also checks for hidden tautologies
 * \param c_nr the ID of the clause which is to be extended.
 * \param clause a copy of the clause as an ordered container into which the
 * hidden literals are inserted. \note Only the copy of the clause is modified,
 * the original clause in clauses[c_nr] stays untouched. \pre Assumes that
 * "_seen" vector was initialized with "clause" \return true if clause is a
 * hidden tautology
 */
bool
Formula::addHiddenLiterals(const int c_nr, Clause::ClauseData& clause, std::uint64_t& sign) const
{
    const bool hidden_binary = addHiddenLiteralsBinary(c_nr, clause, sign);
    if (hidden_binary) {
        return true;
    }  // true = clause is hidden tautology and may be deleted.

    std::size_t csize = clause.size();

    for (std::size_t c = 0; c != csize; ++c) {
        const Literal cur_lit = clause[c];
        addSignatureLit(sign, cur_lit);

        // proceed all clauses in which the current literal occurs
        for (std::size_t j = 0; j != _occ_list[cur_lit].size(); ++j) {
            const ClauseID other_c_nr = _occ_list[cur_lit][j];
            if (clauseDeleted(other_c_nr) || c_nr == static_cast<int>(other_c_nr)) {
                continue;
            }  // skip the clause c
            const Clause& other = _clauses[other_c_nr];

            // We can skip all binary clauses; they have already been handled.
            // Clause "other" can never contain a hidden literal if it is larger or
            // has the same size as current clause
            if (other.size() <= 2 || other.size() > clause.size() + 1) {
                continue;
            }

            // Either this clause was already successfully used -> always skip
            // Or clause was used unsuccessfully -> skip this time until it might get
            // candidate again
            if (other.isMarked()) {
                continue;
            }

            bool    candidate_found = false;
            Literal candidate       = 0;

            _process_limit.decreaseLimitBy(3, clause.size() + other.size());

            // Now check whether we can add a hidden literal by "other"
            for (const Literal lit : other) {
                // Literal not in clause -> candidate for adding hidden literal
                if (_seen[lit] == 0u) {
                    // We have already found a candidate? There can never be more than one
                    // -> clause can not used to add a hidden literal
                    if (candidate != 0) {
                        candidate_found = false;
                        break;
                    }
                    candidate       = lit;
                    candidate_found = true;
                }
            }

            // Mark current clause as touched, so that we do not check it again
            // until it becomes necessary
            other.mark(1);

            // We did not find any candidate literal?
            // That means all literals in "other" are already in "clause"
            // -> "other" subsumes the current clause and we have found a hidden
            // subsumption
            if (candidate == 0) {
                // reset clause marking
                for (ClauseID my_c_nr = 0; my_c_nr <= maxClauseIndex(); ++my_c_nr) {
                    _clauses[my_c_nr].unMark();
                }
                return true;
            }
            // We found a candidate for adding a hidden literal
            // -> add negated literal to clause
            else if (candidate_found) {
                val_assert(candidate != 0 && _seen[candidate] == 0u);
                // If negated candidate is already in clause, we would add a duplicated
                // literal
                // -> skip in this case
                if (_seen[negate(candidate)] == 0u) {
                    ++csize;
                    const Literal hidden_lit = negate(candidate);
                    _seen[hidden_lit]        = 1u;
                    clause.push_back(hidden_lit);

                    // Unmark occurence list of added literal
                    // These clauses are candidates again for producing a new hidden
                    // literal
                    for (std::size_t k = 0; k != _occ_list[hidden_lit].size(); ++k) {
                        _clauses[_occ_list[hidden_lit][k]].unMark(1);
                    }
                }
                // Mark current clause as used, so that we do not check it again for
                // "clause" Since "other" produced a hidden literal, it can never be a
                // candidate for adding another hidden literal
                other.mark(2);
            }
        }
    }

    // reset clause marking
    for (ClauseID my_c_nr = 0; my_c_nr <= maxClauseIndex(); ++my_c_nr) {
        _clauses[my_c_nr].unMark();
    }

    return false;
}

/**
 * \brief Adds covered literals to a given clause.
 * \param clause a copy of the clause as an ordered container into which the
 * covered literals are inserted. \param sign the signature of the clause passed
 * to this function \note The passed clause is modified together with its
 * signature \pre Assumes that "_seen" vector was initialized with "clause"
 * \return true if clause is a hidden tautology
 */
bool
Formula::addCoveredLiterals(Clause::ClauseData& clause, std::uint64_t& sign) const
{
    std::size_t          csize = clause.size();
    std::vector<Literal> candidate_literals;
    for (std::size_t i = 0; i != csize; ++i) {
        const Literal pivot = clause[i];
        addSignatureLit(sign, pivot);
        if (isUniversal(lit2var(pivot))) {
            continue;
        }

        candidate_literals.clear();
        bool first = true;

        const std::vector<ClauseID>& resolution_partners = _occ_list[negate(pivot)];
        for (const ClauseID partner_index : resolution_partners) {
            // Check if clause (resolved)_pivot _clauses[partner_index] is a tautology
            if (clauseDeleted(partner_index)) {
                continue;
            }
            _process_limit.decreaseLimitBy(3, _clauses[partner_index].size());
            if (checkResolventTautology(_clauses[partner_index], lit2var(pivot))) {
                continue;
            }

            if (first) {
                first = false;
                _process_limit.decreaseLimitBy(_clauses[partner_index].size());
                for (const Literal literal : _clauses[partner_index]) {
                    if (literal != negate(pivot) && dependenciesSubset(lit2var(literal), lit2var(pivot))) {
                        val_assert(!_seen2[literal]);
                        _seen2[literal] = 1u;
                        candidate_literals.push_back(literal);
                    }
                }
            } else {
                for (const Literal partner_lit : _clauses[partner_index]) {
                    if (!dependenciesSubset(lit2var(partner_lit), lit2var(pivot))) {
                        continue;
                        // "partner_lit" is part of candidates -> add to intersection
                    } else if (_seen2[partner_lit] != 0u) {
                        // mark "partner_lit" as part of the intersection
                        _seen2[partner_lit] = 'i';
                    }
                }

                // Update candidates by intersection values
                // Update also the _seen2 vector
                unsigned int j = 0;
                for (const auto lit : candidate_literals) {
                    val_assert(_seen2[lit] != 0u);

                    if (_seen2[lit] == 'i') {
                        // literal is part of the intersection
                        candidate_literals[j] = lit;
                        _seen2[lit]           = 1u;
                        ++j;
                    } else {
                        // literal not part of intersection -> unmark literal
                        _seen2[lit] = 0u;
                    }
                }
                candidate_literals.resize(j);

                // Intersection empty -> We are done for the current pivot element
                // At this point all literals in the _seen2 vector are cleaned -> no
                // update needed
                if (candidate_literals.empty()) {
                    break;
                }
            }
        }

        // clear all remaining marked literals
        clearSeen2(candidate_literals);

        for (const Literal literal : candidate_literals) {
            // Found hidden tautology
            if (_seen[negate(literal)] != 0u) {
                return true;
            } else if (_seen[literal] == 0u) {
                clause.push_back(literal);
                ++csize;
                _seen[literal] = 1u;
            }
        }
    }

    return false;
}

}  // end namespace hqspre
