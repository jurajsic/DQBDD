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

#include <iterator>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <easylogging++.hpp>
#include "auxil.hpp"
#include "clause.hpp"
#include "formula.hpp"
#include "literal.hpp"

namespace hqspre {

bool
Formula::pairBlockedClauses()
{
    VLOG(1) << __FUNCTION__;

#ifdef SKOLEM
    LOG(WARNING) << __FUNCTION__ << " does not support Skolem function computation yet.\n";
    return false;
#endif

    val_assert(_unit_stack.empty());

    unsigned int count = 0;

    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (clauseDeleted(c_nr)) continue;
        Clause::ClauseData current_clause(_clauses[c_nr].size());
        for (std::size_t i = 0; i != _clauses[c_nr].size(); ++i) {
            current_clause[i] = _clauses[c_nr][i];
            val_assert(!_seen[current_clause[i]]);
            _seen[current_clause[i]] = true;
        }
        auto sig = _clauses[c_nr].getSignature();
        addHiddenLiteralsBinary(static_cast<int>(c_nr), current_clause, sig);
        clearSeen(current_clause);
        bool blocked = true;
        bool found1  = false;
        bool found2  = false;

        for (auto lit1_pos = current_clause.cbegin(); lit1_pos != current_clause.cend(); ++lit1_pos) {
            const Literal  lit1 = *lit1_pos;
            const Variable var1 = lit2var(lit1);
            if (!isExistential(var1)) continue;
            found1 = true;

            auto lit2_pos = lit1_pos;
            ++lit2_pos;
            for (; lit2_pos != current_clause.cend(); ++lit2_pos) {
                const Literal  lit2 = *lit2_pos;
                const Variable var2 = lit2var(lit2);
                if (!isExistential(var2)) continue;
                found2 = true;

                std::set<ClauseID> res_candidates;
                res_candidates.insert(_occ_list[negate(lit1)].cbegin(), _occ_list[negate(lit1)].cend());
                res_candidates.insert(_occ_list[negate(lit2)].cbegin(), _occ_list[negate(lit2)].cend());

                for (const ClauseID other_c_nr : res_candidates) {
                    const Clause&     other_clause = _clauses[other_c_nr];
                    std::set<Literal> resolvent;
                    resolvent.insert(negate(lit1));
                    resolvent.insert(negate(lit2));
                    bool lit1_found = false;
                    bool lit2_found = false;
                    for (Literal ell : current_clause) {
                        if (ell == lit1 || ell == lit2)
                            continue;
                        else
                            resolvent.insert(ell);
                    }
                    for (Literal ell : other_clause) {
                        if (ell == negate(lit1))
                            lit1_found = true;
                        else if (ell == negate(lit2))
                            lit2_found = true;
                        resolvent.insert(ell);
                    }
                    auto it1 = resolvent.cbegin();
                    auto it2 = it1;
                    ++it2;
                    bool tautology = false;
                    while (it2 != resolvent.cend()) {
                        if (lit2var(*it1) == lit2var(*it2) && *it1 != *it2) {
                            // tautology found;
                            const Variable tvar = lit2var(*it1);
                            if ((!lit1_found || dependenciesSubset(tvar, lit2var(lit1)))
                                || (!lit2_found || dependenciesSubset(tvar, lit2var(lit2)))) {
                                tautology = true;
                                break;
                            }
                        }
                        ++it1;
                        ++it2;
                    }
                    if (!tautology) {
                        blocked = false;
                        break;
                    }
                }
                if (blocked) break;
            }
            if (blocked) break;
        }

        if (blocked && found1 && found2) {
            ++count;
            removeClause(c_nr);
        }
    }

    VLOG(2) << __FUNCTION__ << " removed " << count << " pair-blocked clauses.";

    return count > 0;
}

}  // end namespace hqspre
