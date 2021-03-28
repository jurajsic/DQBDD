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

#include "propagator.hpp"

#include <random>
#include <stack>

#define ELPP_STL_LOGGING
#include <easylogging++.hpp>

#include "auxil.hpp"
#include "clause.hpp"
#include "literal.hpp"
#include "prefix.hpp"

/**
 * \file propagator.cpp
 * \brief Implementation of all BCP-related functions
 * \author Ralf Wimmer
 * \date 11/2017
 */

namespace hqspre {

/**
 * \brief Checks if a literal is unsatisfied by a given assignment.
 * \return true iff the literal is assigned false, false if it is true or
 * unknown.
 */
inline static bool
isUnsatisfied(const Literal lit, const std::vector<TruthValue>& assignment)
{
    return isUnsatisfied(lit, assignment[lit2var(lit)]);
}

/**
 * \brief Checks if a literal is satisfied by a given assignment.
 * \return true iff the literal is assigned true, false if it is false or
 * unknown.
 */
inline static bool
isSatisfied(const Literal lit, const std::vector<TruthValue>& assignment)
{
    return isSatisfied(lit, assignment[lit2var(lit)]);
}

/**
 * \brief Inserts a literal (and all its direct implications) into the queue of
 * literals that are to be propagated.
 *
 * This function also handles the binary clauses (i.e., implications).
 *
 * \throws ConflictException if a conflict occurs (literal was already assigned
 * in the opposite way) \param lit the literal to be enqueued \param the ID of
 * the clause which should be ignored (used e.g. in vivification)
 */
void
SatPropagator::enqueue(const Literal lit, const int except_clause)
{
    const Variable var         = lit2var(lit);
    const bool     already_set = (_assignment[var] != TruthValue::UNKNOWN);

    if (isPositive(lit)) {
        if (_assignment[var] == TruthValue::FALSE)
            throw ConflictException();
        else
            _assignment[var] = TruthValue::TRUE;
    } else {
        // literal is negative
        if (_assignment[var] == TruthValue::TRUE)
            throw ConflictException();
        else
            _assignment[var] = TruthValue::FALSE;
    }

    if (already_set) return;

    _unit_stack.push_back(lit);
    for (std::size_t ptr = _unit_stack.size() - 1; ptr != _unit_stack.size(); ++ptr) {
        const Literal current_lit = _unit_stack[ptr];
        for (const auto bin : _implications[current_lit]) {
            if (static_cast<int>(bin.getClauseID()) == except_clause) continue;
            const Literal  implied_lit = bin.getLiteral();
            const Variable implied_var = lit2var(implied_lit);
            if (_assignment[implied_var] == TruthValue::UNKNOWN) {
                _assignment[implied_var] = makeSatisfied(implied_lit);
                _unit_stack.push_back(implied_lit);
            } else if (isUnsatisfied(implied_lit, _assignment))
                throw ConflictException();
        }
    }
}

/**
 * \brief Creates watch literals for all clauses with at least 3 literals.
 *
 * Binary clauses are handled separately using the implication datastructure.
 * \sa SatPropagator::createWatches(ClauseID)
 */
void
SatPropagator::createWatches()
{
    _watches.clear();
    _watches.resize(_prefix->maxLitIndex() + 1);
    _watchpos.assign(_clauses.size(), std::make_pair(0, 0));

    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        const Clause& clause = _clauses[c_nr];

        if (clause.getStatus() != ClauseStatus::DELETED && clause.size() > 2) {
            createWatches(c_nr);
        }
    }

    val_assert(checkConsistency());
}

/**
 * \brief Creates watch literals for the given clause.
 *
 * The clause must contain at least three literals.
 * \sa SatPropagator::createWatches()
 */
void
SatPropagator::createWatches(const ClauseID c_nr)
{
    const Clause&     clause      = _clauses[c_nr];
    const std::size_t clause_size = clause.size();

    if (clause_size <= 2) return;

    const std::size_t pos1 = (std::rand() % clause_size);
    std::size_t       pos2 = (std::rand() % (clause_size - 1));
    if (pos2 >= pos1) ++pos2;

    _watchpos[c_nr] = std::make_pair(pos1, pos2);
    _watches[clause[pos1]].push_back(c_nr);
    _watches[clause[pos2]].push_back(c_nr);
}

/**
 * \brief Performs SAT-based BCP for a literal.
 *
 * This function searches for all unit literals which result from assigning the
 * given literal. The found implications are stored in `_units`, their
 * satisfying assignment in `_assignment`. If a conflict occurs, i.e., two
 * contradictory literals are found, a ConflictException is raised. The
 * parameter `except_clause` can be used to ignore a certain clause of the
 * database; this is necessary, e.g., for vivification where the vivified clause
 * is ignored.
 *
 * \param assumption the literals which should be propagated
 * \param except_clause the clause which should be ignored during BCP
 * \throw ConflictException if a conflict was detected
 */
void
SatPropagator::bcp(const Literal assumption, const int except_clause)
{
    val_assert(_prefix);
    val_assert(_watches.size() > _prefix->maxLitIndex());
    val_assert(_watchpos.size() >= _clauses.size());
    val_assert(_assignment.size() > _prefix->maxVarIndex());
    val_assert(checkConsistency());

    enqueue(assumption, except_clause);

    while (!_unit_stack.empty()) {
        const Literal neg_lit1 = _unit_stack.back();
        _unit_stack.pop_back();
        _units.push_back(neg_lit1);
        const Literal lit1 = negate(neg_lit1);

        val_assert(isUnsatisfied(lit1, _assignment));

        std::size_t next_free = 0;

        const std::size_t num_watches = _watches[lit1].size();
        for (std::size_t current_pos = 0; current_pos < num_watches; ++current_pos) {
            const ClauseID c_nr = _watches[lit1][current_pos];

            if (static_cast<int>(c_nr) == except_clause) {
                _watches[lit1][next_free] = c_nr;
                ++next_free;
                continue;
            }

            const Clause&     clause      = _clauses[c_nr];
            const std::size_t clause_size = clause.size();
            val_assert(clause_size > 2);

            // Take care that lit1 is the literal over whose watch list we are
            // running; pos1 is its position in the current clause. The other watch
            // literal is lit2 with position pos2.
            std::size_t pos1 = 0;
            std::size_t pos2 = 0;
            Literal     lit2 = 0;
            if (lit1 == clause[_watchpos[c_nr].first]) {
                pos1 = _watchpos[c_nr].first;
                pos2 = _watchpos[c_nr].second;
                lit2 = clause[pos2];
            } else {
                val_assert(lit1 == clause[_watchpos[c_nr].second]);
                pos1 = _watchpos[c_nr].second;
                pos2 = _watchpos[c_nr].first;
                lit2 = clause[pos2];
            }
            val_assert(clause[pos1] == lit1 && clause[pos2] == lit2);
            val_assert(isUnsatisfied(lit1, _assignment));

            if (isSatisfied(lit2, _assignment)) {
                // We have to keep the watch entry -> copy it from current_pos to
                // next_free
                _watches[lit1][next_free] = c_nr;
                ++next_free;
                continue;
            }

            // We have to find a replacement for lit1
            std::size_t new_pos1 = pos1 + 1;
            if (new_pos1 >= clause_size) new_pos1 = 0;
            while (new_pos1 != pos1) {
                if (new_pos1 != pos2 && !isUnsatisfied(clause[new_pos1], _assignment)) {
                    break;
                }
                ++new_pos1;
                if (new_pos1 >= clause_size) new_pos1 = 0;
            }

            if (new_pos1 == pos1) {
                // We could not find a replacement for the first watch -> second watch
                // is unit (or unsatisfied as well).

                try {
                    enqueue(lit2, except_clause);
                } catch (const ConflictException& e) {
                    // We need to leave the _watches structure in a consistent state.
                    // That means we have to copy the remaining entries from 'current_pos'
                    // on to 'next_free' and resize _watches[lit1]. Then re-throw the
                    // exception.
                    if (next_free != current_pos) {
                        while (current_pos < _watches[lit1].size()) {
                            _watches[lit1][next_free++] = _watches[lit1][current_pos++];
                        }
                        _watches[lit1].resize(next_free);
                    }
                    throw;
                }

                // We have to keep the current watch entry -> copy it from current_pos
                // to next_free
                _watches[lit1][next_free] = c_nr;
                ++next_free;
                continue;
            } else {
                // We have to update the first watch literal entry; the entry in
                // _watches[lit1] is automatically deleted.
                val_assert(new_pos1 != pos1);
                val_assert(new_pos1 != pos2);
                val_assert(clause[new_pos1] != lit1 && clause[new_pos1] != lit2);
                val_assert(!isUnsatisfied(clause[new_pos1], _assignment));

                _watchpos[c_nr] = std::make_pair(new_pos1, pos2);
                _watches[clause[new_pos1]].push_back(c_nr);
            }
        }

        _watches[lit1].resize(next_free);
    }  // end while !_unit_stack.empty()

    val_assert(checkConsistency());
}  // end bcp()

/**
 * \brief Clears the unit literals found in previous calls to
 * SatPropagator::bcp(Literal, int).
 *
 * Afterward, all literals are unassigned again such that new BCP operations can
 * be executed.
 */
void
SatPropagator::clear()
{
    std::fill(_assignment.begin(), _assignment.end(), TruthValue::UNKNOWN);
    _units.clear();
    _unit_stack.clear();
}

/**
 * \brief Clears all found units and recreates the watch datastructure.
 */
void
SatPropagator::reset()
{
    clear();
    createWatches();
}

/**
 * \brief Adds a clause to the watched literal data structure.
 *
 * This function must be called <i>after</i> Formula::addClause() with the
 * clause's ID.
 */
void
SatPropagator::addClause(const ClauseID c_nr)
{
    if (_watchpos.size() <= c_nr) _watchpos.resize(c_nr + 1);
    createWatches(c_nr);
    val_assert(checkConsistency());
}

/**
 * \brief Removes a clause from the watched literal data structure.
 *
 * This function must be called <i>before</i> Formula::removeClause(ClauseID).
 */
void
SatPropagator::removeClause(const ClauseID c_nr)
{
    val_assert(_clauses[c_nr].size() > 2);

    const auto    w_pos = _watchpos[c_nr];
    const Literal lit1  = _clauses[c_nr][w_pos.first];
    const Literal lit2  = _clauses[c_nr][w_pos.second];

    val_assert(!_watches[lit1].empty());
    val_assert(!_watches[lit2].empty());

    auto pos = std::find(_watches[lit1].begin(), _watches[lit1].end(), c_nr);
    val_assert(pos != _watches[lit1].end());
    *pos = _watches[lit1].back();
    _watches[lit1].pop_back();

    pos = std::find(_watches[lit2].begin(), _watches[lit2].end(), c_nr);
    val_assert(pos != _watches[lit2].end());
    *pos = _watches[lit2].back();
    _watches[lit2].pop_back();

    val_assert(std::find(_watches[lit1].cbegin(), _watches[lit1].cend(), c_nr) == _watches[lit1].cend());
    val_assert(std::find(_watches[lit2].cbegin(), _watches[lit2].cend(), c_nr) == _watches[lit2].cend());
}

/**
 * \brief Checks if the watch literal datastructure is consistent.
 * \return true if the watch literal datastructure is consistent; false if an
 * error was detected.
 */
bool
SatPropagator::checkConsistency() const
{
    if (!_prefix) {
        LOG(ERROR) << "Prefix is null.";
        return false;
    }

    for (Literal lit = _prefix->minLitIndex(); lit <= _prefix->maxLitIndex(); ++lit) {
        for (ClauseID c_nr : _watches[lit]) {
            if (_clauses[c_nr].getStatus() == ClauseStatus::DELETED) {
                LOG(ERROR) << "Deleted clause #" << c_nr << " in watch list of literal" << lit2dimacs(lit) << '.';
                return false;
            }

            const auto& w = _watchpos[c_nr];
            if (_clauses[c_nr][w.first] != lit && _clauses[c_nr][w.second] != lit) {
                LOG(ERROR) << "Watch points to the wrong literal.";
                LOG(ERROR) << "  _clause = " << _clauses[c_nr];
                LOG(ERROR) << "  lit = " << lit2dimacs(lit);
                LOG(ERROR) << "  pos1 = " << w.first;
                LOG(ERROR) << "  pos2 = " << w.second;
                return false;
            }
        }
    }

    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (_clauses[c_nr].getStatus() == ClauseStatus::DELETED || _clauses[c_nr].size() <= 2) continue;

        const Clause& clause = _clauses[c_nr];

        const auto    w_pos = _watchpos[c_nr];
        const Literal lit1  = clause[w_pos.first];
        const Literal lit2  = clause[w_pos.second];

        if (std::find(_watches[lit1].cbegin(), _watches[lit1].cend(), c_nr) == _watches[lit1].cend()) {
            LOG(ERROR) << "Clause no. " << c_nr << " not in watch list of its first watched literal "
                       << lit2dimacs(lit1) << '.';
            LOG(ERROR) << "_watchpos[" << c_nr << "] = (" << _watchpos[c_nr].first << ", " << _watchpos[c_nr].second
                       << ").";
            LOG(ERROR) << "_clauses[" << c_nr << "] = " << clause;
            LOG(ERROR) << "_watches[" << lit2dimacs(lit1) << "] = ";
            for (auto c : _watches[lit1]) LOG(ERROR) << c;
            return false;
        }

        if (std::find(_watches[lit2].cbegin(), _watches[lit2].cend(), c_nr) == _watches[lit2].cend()) {
            LOG(ERROR) << "Clause no. " << c_nr << " not in watch list of its second watched literal.";
            return false;
        }
    }

    return true;
}

/**
 * \brief Inserts a literal (and all its direct implications) into the queue of
 * literals that are to be propagated.
 *
 * This function also handles the binary clauses (i.e., implications) and
 * applies universal reduction to the binary clauses.
 *
 * \throws ConflictException if a conflict occurs (literal was already assigned
 * in the opposite way) \param lit the literal to be enqueued \param the ID of
 * the clause which should be ignored (used e.g. in vivification)
 */
void
QbfPropagator::enqueue(const Literal lit, const int except_clause)
{
    const Variable var         = lit2var(lit);
    const bool     already_set = (_assignment[var] != TruthValue::UNKNOWN);

    if (isUnsatisfied(lit, _assignment))
        throw ConflictException();
    else
        _assignment[var] = makeSatisfied(lit);

    if (!already_set) {
        _unit_stack.push_back(lit);
        for (std::size_t ptr = _unit_stack.size() - 1; ptr != _unit_stack.size(); ++ptr) {
            const Literal current_lit = _unit_stack[ptr];
            for (const auto bin : _implications[current_lit]) {
                if (static_cast<int>(bin.getClauseID()) == except_clause) continue;

                const Literal  implied_lit = bin.getLiteral();
                const Variable implied_var = lit2var(implied_lit);

                if (_assignment[implied_var] == TruthValue::UNKNOWN) {
                    if (_prefix->isUniversal(implied_var))
                        throw ConflictException();
                    else {
                        _assignment[implied_var] = makeSatisfied(implied_lit);
                        _unit_stack.push_back(implied_lit);
                    }
                } else if (isUnsatisfied(implied_lit, _assignment)) {
                    throw ConflictException();
                }
            }
        }
    }
}

/**
 * \brief Creates watch literals for all clauses with at least 3 literals.
 *
 * Binary clauses are handled separately using the implication datastructure.
 * \sa SatPropagator::createWatches(ClauseID)
 * \sa QbfPropagator::createUnivWatches(ClauseID)
 */
void
QbfPropagator::createWatches()
{
    _watches.clear();
    _watches.resize(_prefix->maxLitIndex() + 1);
    _watchpos.assign(_clauses.size(), std::make_pair(0, 0));

    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        const Clause& clause = _clauses[c_nr];
        if (clause.getStatus() != ClauseStatus::DELETED && clause.size() > 2) {
            bool only_exist = true;
            for (const auto lit : clause) {
                if (_prefix->isUniversal(lit2var(lit))) {
                    only_exist = false;
                    break;
                }
            }

            if (only_exist)
                SatPropagator::createWatches(c_nr);
            else
                createUnivWatches(c_nr);
        }
    }

    val_assert(checkConsistency());
}

void
QbfPropagator::createUnivWatches(const ClauseID /* c_nr */)
{}

}  // end namespace hqspre
