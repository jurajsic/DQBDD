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

#ifndef HQSPRE_BOOL_VECTOR_HPP_
#define HQSPRE_BOOL_VECTOR_HPP_

#include <climits>
#include <cstddef>
#include <cstdint>
#include <iosfwd>

/**
 * \file bool_vector.hpp
 * \brief Header file for BoolVector
 * \author Florian Pigorsch, Ralf Wimmer
 * \date 2007-2016
 */

namespace hqspre {

/**
 * \brief Data structure for compactly storing bit vectors
 *
 * It stores bools with 1 bit per value on average and provides
 * dedicated operators like bit-wise logical operations.
 */
class BoolVector
{
   public:
    using BinType = std::uint64_t;

    BoolVector();
    explicit BoolVector(std::size_t size, bool initial = true);
    BoolVector(const BoolVector& v);
    BoolVector(BoolVector&& v) noexcept;
    ~BoolVector() noexcept;

    BoolVector& operator=(const BoolVector& v);
    BoolVector& operator=(BoolVector&& v) noexcept;

    std::size_t size() const noexcept;
    std::size_t countTrue() const;
    std::size_t getFirstTrue() const;

    bool uninitialized() const noexcept;
    void initialize(std::size_t size, bool initial = true);

    bool get(std::size_t index) const;
    void set(std::size_t index, bool value);
    void flip(std::size_t index);

    void setBin(std::size_t index, BinType b);

    void set(bool value);

    bool allTrue() const;
    bool allFalse() const;
    bool operator==(const BoolVector& b) const;

    BoolVector operator~() const;
    void       flip();

    BoolVector& operator&=(const BoolVector& v);
    BoolVector& operator|=(const BoolVector& v);
    BoolVector& operator^=(const BoolVector& v);

    static bool intersectAndCheckEmpty(BoolVector& v1, const BoolVector& v2);

   protected:
    /**
     * Compile-time constants used in BoolVector.
     */
    constexpr static std::size_t BinSize = CHAR_BIT * sizeof(BinType); /*!< Bit-size of one bin */

    std::size_t _size;
    std::size_t _binCount;
    BinType     _lastBinMask;

    BinType* _bins;

    friend std::ostream& operator<<(std::ostream& os, const BoolVector& b);
};

std::ostream& operator<<(std::ostream& os, const BoolVector& b);

}  // end namespace hqspre

#include "bool_vector.ipp"

#endif /* LRABS_BOOLVECTOR_HH */
