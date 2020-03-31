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

#include <algorithm>
#include <cassert>
#include <ostream>

#include <easylogging++.hpp>
#include "auxil.hpp"
#include "clause.hpp"
#include "formula.hpp"
#include "literal.hpp"
#include "prefix.hpp"

/**
 * \file formula_debug.cpp
 * \brief Functions for debugging formula operations
 * \author Ralf Wimmer
 * \date 2016
 */

namespace hqspre {

/**
 * \brief Check if the formula is in a consistent state.
 * \return true if no inconsistencies could be detected.
 *
 * Executed checks:
 * - the cached number of existential and universal variables is correct
 * - all deleted clauses are empty
 * - all available clause numbers (Formula::_deleted_clause_numbers)
 *   correspond to deleted clauses.
 * - if \f$x\in D_y\f$ then \f$y\in D_x\f$ and vice versa.
 * - if \f$(a,b)\f$ is a binary clause, we have the implications
 *   \f$\neg a\rightarrow b\f$ and \f$\neg b\rightarrow a\f$.
 * - If c_nr occurs in an occurrence list, then clause c_nr is not
 *   deleted and contains the corresponding literal.
 */
bool
Formula::checkConsistency() const
{
    if (!_settings.consistency_check) return true;

    // Check the prefix
    if (!_prefix) {
        LOG(ERROR) << "Prefix is null.";
        return false;
    }

    if (_prefix->type() == PrefixType::QBF && !_qbf_prefix) {
        LOG(ERROR) << "Prefix is a QBF prefix, but qbf_prefix is null.";
        return false;
    }

    if (_prefix->type() == PrefixType::DQBF && !_dqbf_prefix) {
        LOG(ERROR) << "Prefix is a DQBF prefix, but dqbf_prefix is null.";
        return false;
    }

    if (_qbf_prefix && _dqbf_prefix) {
        LOG(ERROR) << "Both qbf_prefix and dqbf_prefix are in use.";
        return false;
    }

    if (!_prefix->checkConsistency()) return false;

    if (!_gates.checkConsistency()) return false;

    for (const Variable var : _deleted_var_numbers) {
        if (!varDeleted(var)) {
            LOG(ERROR) << "Variable marked as deleted by is not.";
            return false;
        }
    }

    // Check all clauses for consistency
    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (clauseDeleted(c_nr)) continue;
        if (!_clauses[c_nr].checkConsistency()) return false;
    }

    // Clauses do not contain deleted variables
    for (const auto& clause : _clauses) {
        for (const Literal lit : clause) {
            if (varDeleted(lit2var(lit))) {
                LOG(ERROR) << "Clause contains a deleted variable (var = " << lit2var(lit) << ").";
                return false;
            }
        }
    }

    // Implications do not contain deleted variables
    for (Literal lit = minLitIndex(); lit <= maxLitIndex(); ++lit) {
        const Variable var = lit2var(lit);
        if (varDeleted(var) && !_implications[lit].empty()) {
            LOG(ERROR) << "Implication list of deleted literal " << lit2dimacs(lit) << " is not empty.";
            return false;
        }

        for (BinaryClause impl : _implications[lit]) {
            Literal literal = impl.getLiteral();
            if (varDeleted(lit2var(literal))) {
                LOG(ERROR) << "Literal " << lit2dimacs(literal) << " is implied, but deleted.";
                return false;
            }
        }
    }

    // All deleted clauses are empty
    for (ClauseID c_nr : _deleted_clause_numbers) {
        if (!_clauses[c_nr].empty()) {
            LOG(ERROR) << "Clause " << c_nr << " is deleted, but not empty.";
            return false;
        }
        if (!clauseDeleted(c_nr)) {
            LOG(ERROR) << "Clause number " << c_nr << " is available but clause not deleted.";
            return false;
        }
    }

    // Consistency of implications
    for (Literal a = minLitIndex(); a <= maxLitIndex(); ++a) {
        for (BinaryClause bclause : _implications[a]) {
            const Literal b = bclause.getLiteral();
            if (_implications[negate(b)].find(BinaryClause(negate(a))) == _implications[negate(b)].cend()) {
                LOG(ERROR) << lit2dimacs(a) << " implies " << lit2dimacs(b) << ", but " << lit2dimacs(negate(b))
                           << " does not imply " << lit2dimacs(negate(a)) << '.';
                return false;
            }
        }
    }

    // All binary clauses are stored in the implication lists
    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (clauseDeleted(c_nr) || _clauses[c_nr].size() != 2) continue;

        if (hasImplication(negate(_clauses[c_nr][0]), _clauses[c_nr][1]) == -1) {
            LOG(ERROR) << "We have the binary clause " << _clauses[c_nr] << ", but the implication "
                       << lit2dimacs(negate(_clauses[c_nr][0])) << " => " << lit2dimacs(_clauses[c_nr][1])
                       << " is missing.";
            return false;
        }

        if (hasImplication(negate(_clauses[c_nr][1]), _clauses[c_nr][0]) == -1) {
            LOG(ERROR) << "We have the binary clause " << _clauses[c_nr] << ", but the implication "
                       << lit2dimacs(negate(_clauses[c_nr][1])) << " => " << lit2dimacs(_clauses[c_nr][0])
                       << " is missing.";
            return false;
        }
    }

    // Consistency of occurrence lists
    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (clauseDeleted(c_nr)) continue;

        for (Literal lit : _clauses[c_nr]) {
            if (std::find(_occ_list[lit].cbegin(), _occ_list[lit].cend(), c_nr) == _occ_list[lit].cend()) {
                LOG(ERROR) << "Literal " << lit2dimacs(lit) << " appears in clause #" << c_nr << ", but " << c_nr
                           << " not in _occ_list[" << lit << "].";
                return false;
            }
        }
    }
    for (Literal lit = minLitIndex(); lit <= maxLitIndex(); ++lit) {
        for (ClauseID c_nr : _occ_list[lit]) {
            if (clauseDeleted(c_nr)) {
                LOG(ERROR) << "Deleted clause #" << c_nr << " appears in the occurrence list of literal "
                           << lit2dimacs(lit) << '.';
                return false;
            }
            if (std::find(_clauses[c_nr].cbegin(), _clauses[c_nr].cend(), lit) == _clauses[c_nr].cend()) {
                LOG(ERROR) << "Clause #" << c_nr << " appears in the occurrence list of literal " << lit2dimacs(lit)
                           << ", but the literal not in the clause.";
                return false;
            }
        }
    }

    // The _seen vector has to be cleared
    for (bool b : _seen) {
        if (b) {
            LOG(ERROR) << "Seen vector is not cleared.";
            return false;
        }
    }
    return true;
}

bool
Formula::checkSeen() const
{
    bool error = false;
    for (Literal lit = minLitIndex(); lit <= maxLitIndex(); ++lit) {
        if (_seen[lit]) {
            LOG(ERROR) << lit2dimacs(lit) << " is already seen.";
            error = true;
        }
    }
    val_assert(!error);
    return error;
}

Literal
Formula::getAssignment(const Literal lit) const
{
    if (_assignment[lit2var(lit)] == TruthValue::TRUE) {
        return lit;
    } else if (_assignment[lit2var(lit)] == TruthValue::FALSE) {
        return negate(lit);
    } else {
        return 0;
    }
}

void
Formula::printFullOccurenceList(const Variable var, std::ostream& stream) const
{
    assert(var <= maxVarIndex());
    printOccurenceList(var2lit(var), stream);
    printOccurenceList(negate(var2lit(var)), stream);
}

void
Formula::printOccurenceList(const Literal lit, std::ostream& stream) const
{
    assert(lit <= maxLitIndex());
    stream << "occurrences of " << lit2dimacs(lit) << " (" << _occ_list[lit].size() << " entries) :\n";
    for (ClauseID c_nr : _occ_list[lit]) {
        stream << "#" << c_nr << ": " << _clauses[c_nr] << '\n';
    }
}

void
Formula::printImplications(const Literal lit, std::ostream& stream) const
{
    assert(lit <= maxLitIndex());
    stream << "implications of " << lit2dimacs(lit) << " (" << _implications[lit].size() << ") :\n";
    stream << lit2dimacs(lit) << " => ";
    for (BinaryClause cl : _implications[lit]) {
        stream << lit2dimacs(cl.getLiteral()) << " ";
    }
    stream << '\n';
}

void
Formula::printAllClauses(std::ostream& stream, const bool print_implications) const
{
    stream << "clause database:\n";
    for (unsigned int i = 0; i != _clauses.size(); ++i) {
        if (clauseDeleted(i)) continue;
        stream << (i + 1) << ": " << _clauses[i] << '\n';
    }

    if (print_implications) {
        stream << "clause implications:\n";
        for (Literal l = minLitIndex(); l <= maxLitIndex(); ++l) {
            printImplications(l, stream);
        }
    }
}

}  // end namespace hqspre
