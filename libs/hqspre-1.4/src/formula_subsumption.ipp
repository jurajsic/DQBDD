#ifndef HQSPRE_FORMULA_SUBSUMPTION_IPP_
#define HQSPRE_FORMULA_SUBSUMPTION_IPP_

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


namespace hqspre {

/**
 * \brief Finds and removes all clauses that are subsumed by `short_clause`
 *
 * If the formula contains clauses \f$C\f$ and \f$C'\f$ such that
 * \f$C\subseteq C'\f$, then \f$C'\f$ can be deleted. The result is
 * a logically equivalent formula.
 *
 * This function checks if `short_clause` subsumes any of the formula's clauses.
 * All subsumed clauses are deleted. If `short_clause` is already part of the
 * formula, we have to take care that we do not delete it as each clause is a
 * subset of itself. To avoid this, the parameter `c_nr` is the clause ID if
 * `short_clause` is contained in the formula (otherwise set `c_nr` to a negative value).
 *
 * \param short_clause the clause for which we check if it subsumes any other clause
 * \param c_nr the ID of `short_clause` if it is part of the formula, -1 otherwise
 * \return the number of deleted clauses
 * \sa Formula::binarySubsumption()
 * \sa Formula::narySubsumption()
 */
inline std::size_t Formula::isBackwardSubsuming(const Clause& short_clause, const int c_nr, bool delete_subsumed)
{
    return isBackwardSubsuming(short_clause.getLiterals(), short_clause.getSignature(), c_nr, delete_subsumed);
}


template <typename Container>
inline std::size_t Formula::isBackwardSubsuming(const Container& short_clause, const std::uint64_t signature, const int c_nr, bool delete_subsumed)
{
    val_assert(!short_clause.empty());

    // List of subsumed clauses.
    std::vector<ClauseID> subsumed;

    // determine a literal which appears in the fewest clauses
    const Literal min_lit = getMinOccLit(short_clause);

    // Mark every literal in short_clause
    setSeen2(short_clause);

    const auto short_clause_size = short_clause.size();

    for (const ClauseID other_c_nr: _occ_list[min_lit]) {
        if (c_nr == static_cast<int>(other_c_nr)) continue;

        const Clause& other_clause = _clauses[other_c_nr];

        // If short clause is longer, it cannot subsume the other one
        if (short_clause.size() > other_clause.size()) continue;
        // Check matching signatures
        if ((signature & ~other_clause.getSignature()) != 0) continue;

        // Count literals in "other_clause" which also appears in "short_clause"
        // If we reached the size of the short clause we are done
        std::size_t count = 0;
        for (const Literal lit: other_clause) {
            if (_seen2[lit]) {
                ++count;
                if (count == short_clause_size) {
                    subsumed.push_back(other_c_nr);
                    break;
                }
            }
        }
    }

    if (delete_subsumed) {
        for (const ClauseID del_c_nr: subsumed) {
            if (short_clause_size == _clauses[del_c_nr].size() && short_clause_size == 2 && c_nr >= 0) {
                // We have a duplicate binary clause in the database. We must take
                // care that the implications are not deleted when deleting the
                // subsumed clause.
                removeClause(del_c_nr);
                addImplication(negate(short_clause[0]), short_clause[1], c_nr);
                addImplication(negate(short_clause[1]), short_clause[0], c_nr);
            } else {
                removeClause(del_c_nr);
            }
        }
        stat(Statistics::SUBSUMPTION) += subsumed.size();
    }

    clearSeen2(short_clause);

    return subsumed.size();
}

} // end namespace hqspre


#endif
