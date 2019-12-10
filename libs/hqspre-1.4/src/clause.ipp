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

#ifndef HQSPRE_CLAUSE_IPP_
#define HQSPRE_CLAUSE_IPP_

/**
 * \file clause.ipp
 * \brief Inline implementation of operations on clauses
 * \author Ralf Wimmer
 * \date 05-06/2016
 */

namespace hqspre {

/**
 * \brief Checks if the given sorted range contains a duplicate element
 */
template<typename Iterator>
inline bool is_unique(Iterator first, Iterator last)
{
    return std::adjacent_find(first, last) == last;
}


/**
 * \brief Creates a clause from a vector of literals
 *
 * \param literals the vector of literals
 * \param needs_check if true the literals get sorted and duplicate literals are removed
 * \param status status of the clause (mandatory, optional, deleted)
 * \pre If needs_check is false, the user has to guarantee that the clause is
 *      sorted and does not contain duplicate literals.
 */
template<typename Allocator>
inline Clause::Clause(const std::vector<Literal, Allocator>& literals, const bool needs_check, const ClauseStatus status):
    _literals(literals.cbegin(), literals.cend()),
    _status(status)
{
    if (needs_check) {
        // sort the literals and remove duplicates
        std::sort(_literals.begin(), _literals.end());
        _literals.erase(std::unique(_literals.begin(), _literals.end()), _literals.end());
    }

    val_assert(std::is_sorted(_literals.cbegin(), _literals.cend()));
    val_assert(is_unique(_literals.cbegin(), _literals.cend()));

    computeSignature();
}


/**
 * \brief Creates a clause from a vector of literals
 *
 * \param literals the vector of literals
 * \param needs_check if true the literals get sorted and duplicate literals are removed
 * \param status status of the clause (mandatory, optional, deleted)
 * \pre If needs_check is false, the user has to guarantee that the clause is
 *      sorted and does not contain duplicate literals.
 */
template<>
inline Clause::Clause(const Clause::ClauseData& literals, const bool needs_check, const ClauseStatus status):
    _literals(literals),
    _status(status)
{
    if (needs_check) {
        // sort the literals and remove duplicates
        std::sort(_literals.begin(), _literals.end());
        _literals.erase(std::unique(_literals.begin(), _literals.end()), _literals.end());
    }

    val_assert(std::is_sorted(_literals.cbegin(), _literals.cend()));
    val_assert(is_unique(_literals.cbegin(), _literals.cend()));

    computeSignature();
}


/**
 * \brief Returns the number of literals in the clause
 */
inline std::size_t Clause::size() const noexcept
{
    return _literals.size();
}

/**
 * \brief Checks if the clause is empty.
 */
inline bool Clause::empty() const noexcept
{
    return _literals.empty();
}


/**
 * \brief Returns a read-only iterator pointing to the first literal and providing read access.
 */
inline Clause::const_iterator Clause::cbegin() const noexcept
{
    return _literals.cbegin();
}

/**
 * \brief Returns a read-only iterator pointing to the first literal and providing read access.
 */
inline Clause::const_iterator Clause::begin() const noexcept
{
    return _literals.begin();
}

/**
 * \brief Returns a read-only iterator pointing beyond the last literal and providing read access.
 */
inline Clause::const_iterator Clause::cend() const noexcept
{
    return _literals.cend();
}

/**
 * \brief Returns a read-only iterator pointing beyond the last literal and providing read access.
 */
inline Clause::const_iterator Clause::end() const noexcept
{
    return _literals.end();
}


/**
 * \brief Returns a read-only reverse iterator pointing to the last literal of the clause.
 */
inline Clause::const_reverse_iterator Clause::crbegin() const noexcept
{
    return _literals.crbegin();
}


/**
 * \brief Returns a read-only reverse iterator pointing before the first literal of the clause.
 */
inline Clause::const_reverse_iterator Clause::crend() const noexcept
{
    return _literals.crend();
}


/**
 * \brief Returns a read-only reverse iterator pointing to the last literal of the clause.
 */
inline Clause::const_reverse_iterator Clause::rbegin() const noexcept
{
    return _literals.rbegin();
}


/**
 * \brief Returns a read-only reverse iterator pointing before the first literal of the clause.
 */
inline Clause::const_reverse_iterator Clause::rend() const noexcept
{
    return _literals.rend();
}


/**
 * \brief Returns an iterator pointing to the first literal and providing read and write access.
 */
inline Clause::iterator Clause::begin() noexcept
{
    return _literals.begin();
}


/**
 * \brief Returns an iterator pointing beyond the last literal and providing read and write access.
 */
inline Clause::iterator Clause::end() noexcept
{
    return _literals.end();
}


/**
 * \brief Returns a reverse iterator pointing to the last literal and providing read and write access.
 */
inline Clause::reverse_iterator Clause::rbegin() noexcept
{
    return _literals.rbegin();
}


/**
 * \brief Returns a reverse iterator pointing before the first literal and providing read and write access.
 */
inline Clause::reverse_iterator Clause::rend() noexcept
{
    return _literals.rend();
}

/**
 * \brief Returns the first literal in the clause.
 * \pre The clause must not be empty.
 */
inline Literal Clause::front() const noexcept
{
    return _literals.front();
}


/**
 * \brief Returns the last literal in the clause.
 * \pre The clause must not be empty.
 */
inline Literal Clause::back() const noexcept
{
    return _literals.back();
}


inline const Clause::ClauseData& Clause::getLiterals() const noexcept
{
    return _literals;
}

/**
 * \brief Provides read access to the literals.
 *
 * The parameter `pos` must be in the right range. Otherwise the behavior
 * is undefined (similar to access to a std::vector).
 * \return the literal at the specified position
 */
inline Literal Clause::operator[](const std::size_t pos) const noexcept
{
    return _literals[pos];
}

/**
 * \brief Compares two clauses for equality
 *
 * Two clauses are equal if they contain the same literals.
 */
inline bool Clause::operator==(const Clause& other) const noexcept
{
    return this->_signature == other._signature && this->_literals == other._literals;
}


/**
 * \brief Compares two clauses for disequality
 *
 * Two clauses are not equal if they contain different sets of literals.
 */
inline bool Clause::operator!=(const Clause& other) const noexcept
{
    return this->_signature != other._signature || this->_literals != other._literals;
}


/**
 * \brief Returns the status of a clause -- MANDATORY, OPTIONAL, or DELETED.
 */
inline ClauseStatus Clause::getStatus() const noexcept
{
    return _status;
}


/**
 * \brief Sets the status of a clause -- MANDATORY, OPTIONAL, or DELETED.
 */
inline void Clause::setStatus(const ClauseStatus new_status) noexcept
{
    _status = new_status;
}


/**
 * \brief Returns the signature of a clause.
 *
 * This signature is used for efficiently excluding candidates
 * of subsumption.
 */
inline std::uint64_t Clause::getSignature() const noexcept
{
    return _signature;
}


/**
 * \brief Returns the signature of a clause w.r.t. the variables appearing in the clause.
 *
 * This signature is used for efficiently excluding candidates
 * of subsumption.
 */
inline std::uint32_t Clause::getVarSignature() const noexcept
{
    return _var_signature;
}

/**
 * \brief Returns whether the clause was marked
 *
 */
inline bool Clause::isMarked() const noexcept
{
    return _marked != 0;
}

/**
 * \brief Returns whether a certain position of the clause was marked
 * \note The position of the marker must be in the range from 0 to 7.
 */
inline bool Clause::isMarked(const unsigned int pos) const noexcept
{
    val_assert(pos < 32);
    return static_cast<bool>(_marked & (1 << pos));
}

/**
 * \brief Mark a position of the clause
 * \note The position of the marker must be in the range from 0 to 7.
 */
inline void Clause::mark(const unsigned int pos) const noexcept
{
    val_assert(pos < 8);
    _marked |= (1u << pos);
}

/**
 * \brief Unmark a clause as touched
 *
 */
inline void Clause::unMark() const noexcept
{
    _marked = 0u;
}

/**
 * \brief Unmark a certain position of the clause
 * \note The position of the marker must be in the range from 0 to 7.
 */
inline void Clause::unMark(const unsigned int pos) const noexcept
{
    val_assert(pos < 8);
    _marked &= ~ static_cast<std::uint8_t>(1u << pos);
}


/**
 * \brief Constructs a new BinaryClause object.
 * \param lit the implied literal
 * \param c_nr the ID of the corresponding regular clause
 */
inline BinaryClause::BinaryClause(const Literal lit, const ClauseID c_nr) noexcept :
  _literal(lit),
  _clause_id(c_nr)
{}


/**
 * \brief Checks two BinaryClause objects for equality.
 *
 * Two such objects are equal if they contain the same implied literal,
 * independent of the corresponding clause ID.
 */
inline bool BinaryClause::operator==(const BinaryClause& other) const noexcept
{
    return this->_literal == other._literal;
}

/**
 * \brief Checks two BinaryClause objects for disequality.
 *
 * Two such objects are equal if they contain the same implied literal,
 * independent of the corresponding clause ID.
 */
inline bool BinaryClause::operator!=(const BinaryClause& other) const noexcept
{
    return this->_literal != other._literal;
}

/**
 * \brief Provides a linear order on BinaryClause objects.
 *
 * Only the implied literals are compared, not the clause IDs.
 */
inline bool BinaryClause::operator<(const BinaryClause& other) const noexcept
{
    return this->_literal < other._literal;
}


/**
 * \brief Returns the implied literal of the BinaryClause object.
 */
inline Literal BinaryClause::getLiteral() const noexcept
{
    return _literal;
}


/**
 * \brief Returns the ID of the corresponding regular clause.
 */
inline ClauseID BinaryClause::getClauseID() const noexcept
{
    return _clause_id;
}


} // end namespace hqspre


#endif
