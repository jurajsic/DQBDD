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

#define ELPP_STL_LOGGING
#include <easylogging++.hpp>
#include "aux.hpp"
#include "clause.hpp"
#include "formula.hpp"
#include "literal.hpp"
#include "timer.hpp"
#include "varheap.hpp"

/**
 * \file formula_subsumption.cpp
 * \brief Implementation of subsumption checking
 * \author Ralf Wimmer
 * \date 02/2016
 */

//#define COMPLETE_FORWARDSUB

namespace hqspre {

/**
 * \brief Checks whether a clause is subsumed by an existing binary clause in
 * the database
 *
 * If the formula contains clauses \f$C\f$ and \f$C'\f$ with \f$C\neq C'\f$ such that
 * \f$C\subset C'\f$, then \f$C'\f$ can be deleted. The result is
 * a logically equivalent formula.
 *
 * This function checks if `clause` is subsumed by any of the formulas binary
 * clauses. \param clause the clause for which we check if it subsumed \param
 * except ID of the clause `clause` which should be skipped (as each clause
 * subsumes itself) \return true if the clause is subsumed
 */
template<typename Container>
bool
Formula::isForwardSubsumedByBinary(const Container& clause, int except)
{
    for (const Literal check_lit : clause) {
        const Literal neg_lit = negate(check_lit);
        _process_limit.decreaseLimitBy(_implications[neg_lit].size());

        for (const BinaryClause& bin_clause : _implications[neg_lit]) {
            // Skip explicit stated exception clause
            const ClauseID clauseid = bin_clause.getClauseID();
            if (static_cast<int>(clauseid) == except) continue;

            const Literal sec_lit = bin_clause.getLiteral();
            // If literal-ID is already larger than largest ID of current clause
            // the binary can never subsume the current clause
            if (sec_lit > *clause.rbegin()) continue;

            // Found subsumption
            if (_seen[sec_lit]) {
                VLOG(3) << __FUNCTION__ << " clause " << lit2dimacs(check_lit) << " " << lit2dimacs(sec_lit)
                        << " subsumes " << clause << '.';

                val_assert(!clauseDeleted(clauseid));
                val_assert(_clauses[clauseid].size() == 2);
                _clauses[clauseid].setStatus(ClauseStatus::MANDATORY);
                return true;
            }
        }
    }

    return false;
}

/**
 * \brief Checks whether a clause is subsumed by an existing non-binary clause
 * in the database
 *
 * If the formula contains clauses \f$C\f$ and \f$C'\f$, \f$C\neq C'\f$ such that
 * \f$C\subset C'\f$, then \f$C'\f$ can be deleted. The result is
 * a logically equivalent formula.
 *
 * This function checks if `clause` is subsumed by any of the formulas
 * non-binary clauses. \param clause the clause for which we check if it
 * subsumed \param sign the signatur of `clause`, due to efficency reasons it is
 * assumed that this value is computated in advance \param except ID of the
 * clause `clause` which should be skipped (as each clause subsumes itself)
 * \return true if the clause is subsumed
 */
template<typename Container>
bool
Formula::isForwardSubsumed(const Container& clause, const std::uint64_t sign, int except)
{
#ifndef NDEBUG
    for (Literal lit : clause) {
        val_assert(_seen[lit]);
    }
#endif

    if (isForwardSubsumedByBinary(clause, except)) return true;

    const Literal check_lit = getMinOccLit(clause);

    for (const ClauseID other_c_nr : _occ_list[check_lit]) {
        // Skip explicit stated exception clause
        if (except == static_cast<int>(other_c_nr)) continue;

        const Clause& other_clause = _clauses[other_c_nr];

        // Skip to short clauses
        if (other_clause.size() > clause.size()) continue;
        // Check matching signatures
        if ((other_clause.getSignature() & ~sign) != 0) continue;

        bool is_subsumed = true;

        for (const Literal lit : other_clause) {
            // If there is a literal in "other_clause" which is not contained in
            // "clause" "other_clause" cannot subsume "clause"
            if (!_seen[lit]) {
                is_subsumed = false;
                break;
            }
        }

        if (is_subsumed) {
            VLOG(3) << __FUNCTION__ << " clause " << other_clause << " subsumes " << clause << '.';
            return true;
        }
    }
    return false;
}

// explicit instantiation of template arguments for this function
template bool Formula::isForwardSubsumed(const std::set<Literal>& clause, const std::uint64_t sign, int except);
template bool Formula::isForwardSubsumed(const Clause::ClauseData& clause, const std::uint64_t sign, int except);

bool
Formula::isForwardSubsumed(const Clause& clause, int except)
{
    return isForwardSubsumed(clause.getLiterals(), clause.getSignature(), except);
}

/**
 * \brief Finds and removes subsumed clauses.
 *
 * If the formula contains clauses \f$C\f$ and \f$C'\f$, \f$C\neq C'\f$ such that
 * \f$C\subset C'\f$, then \f$C'\f$ can be deleted. The result is
 * a logically equivalent formula.
 * \return true iff the formula was modified.
 * \sa Formula::binarySubsumption()
 * \sa Formula::narySubsumption()
 */
bool
Formula::removeSubsumedClauses()
{
    VLOG(1) << __FUNCTION__;

    ScopeTimer subsumption(getTimer(WhichTimer::SUBSUMPTION));

    const std::size_t old_stat_subsumption = stat(Statistics::SUBSUMPTION);

    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (clauseDeleted(c_nr)) continue;
        const auto num_del = isBackwardSubsuming(_clauses[c_nr], static_cast<int>(c_nr), true);
        if (num_del > 0) {
            _clauses[c_nr].setStatus(ClauseStatus::MANDATORY);
        }
    }

    VLOG(2) << __FUNCTION__ << " found " << (stat(Statistics::SUBSUMPTION) - old_stat_subsumption)
            << " subsumed clauses.";
    return (stat(Statistics::SUBSUMPTION) > old_stat_subsumption);
}

/**
 * \brief Returns true if 'clause' contains the literal 'other_lit'
 *
 * The clause needs to be in increasing order.
 * \sa Formula::selfSubsumingResolution()
 */
static bool
almostSubsumedByBinary(const Clause& clause, const Literal other_lit)
{
    for (const Literal lit : clause) {
        if (lit == other_lit) {
            return true;
        } else if (lit > other_lit) {
            return false;
        }
    }

    return false;
}

/**
 * \brief Returns true if the clause c2 is -- with the exception of almostLit --
 * a subset of c.
 *
 * In this case resolution of c and c2 w.r.t. almostLit returns a subset of c.
 * \sa Formula::selfSubsumingResolution()
 */
static bool
almostSubsumedByNary(const Clause& c, const Clause& c2, const Literal almostLit)
{
    val_assert(c.size() >= c2.size());

    unsigned int i = 0;
    unsigned int j = 0;

    const auto c_size  = c.size();
    const auto c2_size = c2.size();

    while (i != c_size && j != c2_size) {
        if (c[i] == c2[j]) {
            ++i;
            ++j;
        } else if (c[i] == almostLit) {
            ++i;
        } else if (c2[j] == negate(almostLit)) {
            ++j;
        } else if (c[i] < c2[j]) {
            ++i;
        } else {
            return false;
        }
    }

    if (j != c2_size) return false;
    return true;
}

/**
 * \brief Tries to make clauses shorter by resolution.
 *
 * If c1 and c2 are clauses such that c1 contains literal l and
 * c2 !l and the resolvent of c1 and c2 w.r.t. l is a subset of c1,
 * then c1 can be replaced by the resolvent. This is known as
 * self-subsuming resolution.
 * This method tries to make clauses as short as possible by applying
 * self-subsuming resolution.
 * \return true if the formula was modified.
 */
bool
Formula::selfSubsumingResolution()
{
    // For all n-nary clauses c = (l1,l2,...,ln)
    // check if there exists a clause c' = (-l1,l2,...,lm) with m<=n
    // -> replace c by c'' = (l2,...ln)

    VLOG(1) << __FUNCTION__;

    ScopeTimer self_sub_timer(getTimer(WhichTimer::SELF_SUBSUMPTION));
    _process_limit.setLimit(PreproMethod::SELF_SUBSUMPTION);

    const std::size_t old_stat_selfsubsuming = stat(Statistics::SELF_SUBSUMPTION);
    const std::size_t old_stat_subsumption   = stat(Statistics::SUBSUMPTION);
    std::size_t       subsumed_clauses       = 0;

    // collect all candidates
    _candidates.clear();
    _candidates.resize(_clauses.size());
    _variable_score.clear();
    _variable_score.resize(_clauses.size(), 0);

    for (ClauseID cl_nr = 0; cl_nr != _clauses.size(); ++cl_nr) {
        if (clauseDeleted(cl_nr)) continue;
        if (clauseOptional(cl_nr)) continue;
        // ignore all clauses with size 2 -> they cannot be candidates for
        // selfsubsumption
        if (_clauses[cl_nr].size() < 3) continue;
        _variable_score[cl_nr] = static_cast<int>(_clauses[cl_nr].size());
        _candidates.insert(static_cast<Variable>(cl_nr));
    }

    // Iterate over all candidates
    while (!_candidates.empty()) {
        if (_process_limit.reachedLimit()) {
            VLOG(2) << "Terminate " << __FUNCTION__ << " due to process limit.";
            break;
        }
        if (_interrupt) break;

        bool           subsumed = false;
        const ClauseID c_nr     = _candidates.top();

        for (std::size_t i = 0; i != _clauses[c_nr].size(); ++i) {
            // Clause could be potentially deleted due to universal reduction after
            // removing a literal due to self-subsumption
            if (clauseDeleted(c_nr)) break;
            const auto& c = _clauses[c_nr];

            // Skip binary clauses, they can never be self-subsuming
            // Clause can get binary due selfsubsumption
            if (c.size() < 3) break;

            const Literal x = c[i];

            // self subsumption with binary clause c'
            for (BinaryClause clause : _implications[x]) {
                const Literal impl = clause.getLiteral();

                std::uint32_t binsign = 0;
                addSignatureLit(binsign, lit2var(x));
                addSignatureLit(binsign, lit2var(impl));

                // Skip if signatures do not match
                if ((binsign & ~c.getVarSignature()) != 0) {
                    continue;
                }

                _process_limit.decreaseLimitBy(c.size());
                if (almostSubsumedByBinary(c, impl)) {
                    ++stat(Statistics::SELF_SUBSUMPTION);
                    VLOG(3) << __FUNCTION__ << " removes " << lit2dimacs(x) << " from clause " << c << " due to binary "
                            << lit2dimacs(negate(x)) << " " << lit2dimacs(impl);

                    // Reset index if literal was further reduced by universal reduction
                    // Otherwise
                    if (removeLiteral(c_nr, x))
                        i = -1;
                    else
                        --i;
                    subsumed                = true;
                    const ClauseID bin_c_nr = clause.getClauseID();
                    _clauses[bin_c_nr].setStatus(ClauseStatus::MANDATORY);
                    break;
                }
            }  // end for (clause)

            // self subsumption with n-nary clause c'
            if (!subsumed) {
                const Literal x_neg = negate(x);
                for (const ClauseID c_nr2 : _occ_list[x_neg]) {
                    const auto& c2 = _clauses[c_nr2];
                    if (c2.size() > c.size()) continue;
                    // We already have checked binary clauses
                    if (c2.size() < 3) continue;
                    // Skip if signatures do not match
                    if ((c2.getVarSignature() & ~c.getVarSignature()) != 0) {
                        continue;
                    }

                    _process_limit.decreaseLimitBy(c.size() + c2.size());
                    if (almostSubsumedByNary(c, c2, x)) {
                        ++stat(Statistics::SELF_SUBSUMPTION);
                        VLOG(3) << __FUNCTION__ << " removes " << lit2dimacs(x) << " from clause " << c
                                << " due to clause " << c2;

                        // Reset index if literal was further reduced by universal reduction
                        // Otherwise
                        if (removeLiteral(c_nr, x))
                            i = -1;
                        else
                            --i;
                        subsumed = true;
                        _clauses[c_nr2].setStatus(ClauseStatus::MANDATORY);
                        break;
                    }
                }
            }
        }

        if (subsumed) {
            // Check whether reduced clause now subsumes other clauses
            // Clause size can be 0, if the clause is reduced to a unit due to
            // universal reduction In this case we have nothing to do here
            if (_clauses[c_nr].size() > 0) {
                isBackwardSubsuming(_clauses[c_nr], static_cast<int>(c_nr), true);
                // Add new candidates which can be selfsubsumed with the newly generated
                // clause
                for (Literal lit : _clauses[c_nr]) {
                    for (const ClauseID occ_nr : _occ_list[negate(lit)]) {
                        if (!_candidates.inHeap(static_cast<Variable>(occ_nr))) {
                            _candidates.insert(static_cast<Variable>(occ_nr));
                        }
                    }
                }
            }
            ++subsumed_clauses;
        }
    }

    VLOG(2) << __FUNCTION__ << " found " << (stat(Statistics::SELF_SUBSUMPTION) - old_stat_selfsubsuming)
            << " self-subsumed literals in " << subsumed_clauses << " clauses and "
            << (stat(Statistics::SUBSUMPTION) - old_stat_subsumption) << " subsumptions.";

    return (subsumed_clauses > 0);
}

}  // end namespace hqspre
