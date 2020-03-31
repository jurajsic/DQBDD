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

#ifndef HQSPRE_CLAUSE_HPP_
#define HQSPRE_CLAUSE_HPP_

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <type_traits>
#include <vector>

// #include <boost/pool/pool_alloc.hpp>

#include "auxil.hpp"
#include "literal.hpp"

/**
 * \file clause.hpp
 * \author Ralf Wimmer
 * \date 05/2016
 * \brief Data structures related to clauses
 */

namespace hqspre {

/**
 * \brief Status of a clause (mandatory/optional/deleted)
 */
enum class ClauseStatus : short
{
    /**
     * \brief Actual problem clause
     *
     * Mandatory clauses may not be deleted.
     */
    MANDATORY,

    /**
     * \brief Redundant clause which may be deleted if not used.
     *
     * Optional clauses must be implications from the mandatory clauses.
     * They can be deleted without changing the set of satisfying assignments
     * of the (mandatory) matrix.
     */
    OPTIONAL,

    /**
     * \brief A clause which has been deleted.
     *
     * Deleted clauses may not be accessed anymore and will
     * be recycled in future.
     */
    DELETED
};

/**
 * \brief Add the signature bit of a single literal to signature
 */
template<typename T>
inline void
addSignatureLit(T& signature, const Literal literal) noexcept
{
    static_assert(std::is_integral<T>::value, "addSignatureLit(T&, Literal) only works for integral types T");

    constexpr std::size_t num_bits = sizeof(signature) * CHAR_BIT;
    signature |= static_cast<T>(1) << (literal % num_bits);
}

// Forward declaration
class Formula;

/**
 * \brief Representation of a clause, i.e., a disjunction of literals.
 */
class Clause
{
   public:
    ///\brief Data structure holding the literals
    using ClauseData = std::vector<Literal>;
    //    using ClauseData = std::vector<Literal,
    //    boost::fast_pool_allocator<Literal>>;

    ///\brief Type of a forward (read/write) iterator
    using iterator = ClauseData::iterator;

    ///\brief Type of a forward read-only iterator
    using const_iterator = ClauseData::const_iterator;

    ///\brief Type of a reverse (read/write) iterator
    using reverse_iterator = ClauseData::reverse_iterator;

    ///\brief Type of a reverse read-only iterator
    using const_reverse_iterator = ClauseData::const_reverse_iterator;

    /**
     * \brief Creates a copy of a given clause.
     */
    Clause(const Clause&) = default;

    /**
     * \brief Moves a given clause into a new one.
     */
    Clause(Clause&&) noexcept = default;

    template<typename Allocator>
    explicit Clause(const std::vector<Literal, Allocator>& literals, bool needs_check = true,
                    ClauseStatus status = ClauseStatus::MANDATORY);
    explicit Clause(ClauseData&& literals, bool needs_check = true,
                    ClauseStatus status = ClauseStatus::MANDATORY) noexcept;

    /// \brief Destructor with default behavior
    ~Clause() noexcept = default;

    /// \brief Assignment operator with default behavior
    Clause& operator=(const Clause&) = default;

    /// \brief Move assignment operator with default behavior
    Clause& operator=(Clause&&) noexcept = default;

    std::size_t size() const noexcept;
    bool        empty() const noexcept;
    bool        isTautology() const noexcept;
    bool        containsLiteral(Literal lit) const;

    template<typename Container>
    bool subsetOf(const Container& other, std::uint64_t other_signature) const;

    //@{
    /**
     * \name Getting iterators for traversing the clause
     */
    const_iterator cbegin() const noexcept;
    const_iterator cend() const noexcept;
    const_iterator begin() const noexcept;
    const_iterator end() const noexcept;

    const_reverse_iterator crbegin() const noexcept;
    const_reverse_iterator crend() const noexcept;
    const_reverse_iterator rbegin() const noexcept;
    const_reverse_iterator rend() const noexcept;
    //@}

    Literal           front() const noexcept;
    Literal           back() const noexcept;
    Literal           operator[](std::size_t pos) const noexcept;
    const ClauseData& getLiterals() const noexcept;

    bool operator==(const Clause& other) const noexcept;
    bool operator!=(const Clause& other) const noexcept;

    ClauseStatus  getStatus() const noexcept;
    void          setStatus(ClauseStatus new_status) noexcept;
    std::uint64_t getSignature() const noexcept;
    std::uint32_t getVarSignature() const noexcept;

    bool isMarked() const noexcept;
    bool isMarked(unsigned int pos) const noexcept;
    void mark(unsigned int pos) const noexcept;
    void unMark() const noexcept;
    void unMark(unsigned int pos) const noexcept;

    bool checkConsistency() const;

   private:
    // Private methods; all methods which provide write access to the
    // clause are private.

    iterator         begin() noexcept;
    iterator         end() noexcept;
    reverse_iterator rbegin() noexcept;
    reverse_iterator rend() noexcept;

    ///\brief Gives (read and write) access to the literal at position `pos`.
    Literal& operator[](const std::size_t pos) noexcept { return _literals[pos]; }

    ///\brief Makes the clause empty.
    void clear()
    {
        _literals.clear();
        _status        = ClauseStatus::DELETED;
        _signature     = 0;
        _var_signature = 0;
        _marked        = 0;
    }

    ///\brief Resizes the underlying vector of literals to the given size.
    void resize(const std::size_t new_size) { _literals.resize(new_size, 0); }

    ///\brief Reserves memory for (at least) the given number of literals.
    void reserve(const std::size_t new_capacity) { _literals.reserve(new_capacity); }

    ///\brief Erases a range of literals from the clause.
    iterator erase(iterator begin, iterator end) { return _literals.erase(begin, end); }

    void computeSignature();

    // private data
    ClauseData   _literals;  ///< Vector of literals
    ClauseStatus _status;    ///< Status of the clause (mandatory/optional/deleted)

    /**
     * \brief The signature of the clause.
     * signature is the clause signature produced by the IDs of the literals
     * var_signature is produced by the IDs of the variables
     * \sa Clause::computeSignature()
     */
    std::uint64_t _signature     = 0;
    std::uint32_t _var_signature = 0;

    /**
     * \brief Marking for the clause
     * Used by various methods to mark whether a clause was touched/edited during
     * the process
     */
    mutable unsigned int _marked = 0;

    // friend declarations
    friend class Formula;
};

Clause        resolve(const Clause& c_pos, const Clause& c_neg, Variable var);
std::ostream& operator<<(std::ostream& stream, const Clause& clause);

/**
 * \brief Representation of a binary clause, i.e., the implied literal and its
 * clause id
 *
 * The BinaryClause objects are used to store the implication graph
 * of the binary clauses. For a binary clause \f$\{a,b\}\f$ with ID
 * \f$i\f$ we have two BinaryClause objects: One is stored in `implications`
 * at position \f$\neg a\f$ and contains the implied literal \f$b\f$
 * together with the ID \f$i\f$. The other one is stored at \f$\neg b\f$
 * and contains \f$a\f$ and the ID \f$i\f$.
 */
class BinaryClause
{
   public:
    explicit BinaryClause(Literal lit, ClauseID c_nr = 0) noexcept;
    bool operator==(const BinaryClause& other) const noexcept;
    bool operator!=(const BinaryClause& other) const noexcept;
    bool operator<(const BinaryClause& other) const noexcept;

    Literal  getLiteral() const noexcept;
    ClauseID getClauseID() const noexcept;

   private:
    Literal  _literal;    ///< The implied literal
    ClauseID _clause_id;  ///< The ID of the corresponding regular clause
};

}  // end namespace hqspre

namespace std {

/**
 * \brief Hash functor for binary clauses
 */
template<>
struct hash<hqspre::BinaryClause>
{
    std::size_t operator()(const hqspre::BinaryClause& c) const { return std::hash<hqspre::Literal>()(c.getLiteral()); }
};
}  // namespace std

#include "clause.ipp"

#endif
