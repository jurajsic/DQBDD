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

#ifndef HQSPRE_BOOL_VECTOR_IPP_
#define HQSPRE_BOOL_VECTOR_IPP_

/**
 * \file bool_vector.ipp
 * \brief Implementation of inline functions for the BoolVector datastructure
 * \author Florian Pigorsch, Ralf Wimmer
 * \date 2007, 2017
 */

namespace hqspre {

inline std::size_t BoolVector::size() const noexcept
{
    return _size;
}

inline bool BoolVector::uninitialized() const noexcept
{
    return _bins == nullptr;
}

inline bool BoolVector::get( const std::size_t index ) const
{
    return ( _bins[ ( index / BinSize ) ] >> ( index % BinSize ) ) & 1ul;
}

inline void BoolVector::setBin( const std::size_t index, BoolVector::BinType b )
{
    _bins[ index ] = b;
}

} // end namespace hqspre

#endif
