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

#ifndef HQSPRE_PROPAGATOR_HPP_
#define HQSPRE_PROPAGATOR_HPP_

#include <exception>
#include <unordered_set>
#include <vector>

#include "clause.hpp"
#include "prefix.hpp"


/**
 * \file propagator.hpp
 * \brief Contains the definitons of all BCP-related classes and functions
 * \author Ralf Wimmer
 * \date 11/2017
 */


namespace hqspre {


/**
 * \brief Exception that is thrown when a conflict occurs during BCP
 */
class ConflictException: public std::exception
{
public:

    /**
     * \brief Default constructor
     */
    ConflictException() = default;

    /**
     * \brief Default destructor
     */
    virtual ~ConflictException() noexcept = default;
};


class QbfPropagator;

/**
 * \brief SAT-based Boolean constraint propagation (BCP)
 *
 * This class implements the two-watch-literal scheme for SAT for efficient unit literal detection.
 * In contrast to QBF-based BCP, SAT-based BCP does not perform universal reduction to detect more
 * unit literals.
 */
class SatPropagator
{
public:
    SatPropagator(const Prefix* const prefix,
                  const std::vector<Clause>& clauses,
                  const std::vector<std::unordered_set<BinaryClause>>& implications) noexcept:
        _clauses(clauses),
        _implications(implications),
        _prefix(prefix),
        _units(),
        _assignment(prefix->maxLitIndex() + 1, TruthValue::UNKNOWN),
        _unit_stack(),
        _watches(),
        _watchpos()
    {
        createWatches();
        _unit_stack.reserve(_prefix->maxVarIndex() + 1);
    }

    SatPropagator(const SatPropagator& other) = default;
    SatPropagator(SatPropagator&& other) = default;
    ~SatPropagator() noexcept = default;
    SatPropagator& operator=(const SatPropagator& other) = default;
    SatPropagator& operator=(SatPropagator&& other) = default;

    void bcp(Literal assumption, int except_clause = -1);
    void removeClause(ClauseID c_nr);
    void addClause(ClauseID c_nr);
    void clear();
    void reset();

    /**
     * \brief Returns a reference to the list of found unit clauses (including the initial assumptions)
     */
    const std::vector<Literal>&    getUnits()      const noexcept { return _units; }
    std::vector<Literal>&          getUnits()      noexcept { return _units; }
    const std::vector<TruthValue>& getAssignment() const noexcept { return _assignment; }

private:

    friend class QbfPropagator;

    void createWatches();
    void createWatches(ClauseID c_nr);
    void enqueue(Literal lit, int except_clause = -1);
    bool checkConsistency() const;

    /** \brief Reference to the clause database of the formula */
    const std::vector<Clause>& _clauses;

    /** \brief Reference to the implication datastructure of the formula */
    const std::vector<std::unordered_set<BinaryClause>>& _implications;

    /** \brief Pointer to the formula's prefix */
    const Prefix* const _prefix;

    /** \brief Target vector for storing the found units. */
    std::vector<Literal> _units;

    /** \brief Assignment for the found unit literals (provided there was no conflict) */
    std::vector<TruthValue> _assignment;

    /** \brief Stores the unprocessed unit literals */
    std::vector<Literal> _unit_stack;

    /** \brief Watch literal datastructure. _watches[lit] contains the IDs of all clauses in which `lit` is watched. */
    std::vector<std::vector<ClauseID>> _watches;

    /** \brief Positions of the watched literals in each clause. */
    std::vector<std::pair<std::size_t, std::size_t>> _watchpos;
};


class QbfPropagator: private SatPropagator
{
public:
    QbfPropagator(const Prefix* const prefix,
                  const std::vector<Clause>& clauses,
                  const std::vector<std::unordered_set<BinaryClause>>& implications) noexcept:
        SatPropagator(prefix, clauses, implications)
    {
    }

    ~QbfPropagator() = default;

    using SatPropagator::clear;

private:
    void enqueue(Literal lit, int except_clause = -1);
    void createWatches();
    void createUnivWatches(ClauseID c_nr);
};


} // end namespace hqspre
#endif

