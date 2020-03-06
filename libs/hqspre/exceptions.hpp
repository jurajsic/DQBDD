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

#ifndef HQSPRE_EXCEPTIONS_HPP_
#define HQSPRE_EXCEPTIONS_HPP_

#include <exception>
#include <string>

namespace hqspre {

class FileFormatException : public std::runtime_error
{
   public:
    explicit FileFormatException() noexcept : std::runtime_error("") {}

    explicit FileFormatException(const char* msg) noexcept : std::runtime_error(msg) {}

    ~FileFormatException() noexcept override = default;
};

/**
 * \brief Exception that is thrown when the formula is recognized to be
 * unsatisfiable.
 */
class UNSATException : public std::runtime_error
{
   public:
    explicit UNSATException() noexcept : std::runtime_error("") {}

    /**
     * \brief Creates a new UNSATException with a description of its reason
     * \param msg description of the exception's reaon
     */
    explicit UNSATException(const char* msg) noexcept : std::runtime_error(msg) {}

    ~UNSATException() noexcept override = default;
};

/**
 * \brief Exception that is thrown when the formula is recognized to be
 * satisfiable.
 */
class SATException : public std::runtime_error
{
   public:
    explicit SATException() noexcept : std::runtime_error("") {}

    /**
     * \brief Creates a new SATException with a description of its reason
     * \param msg description of the exception's reaon
     */
    explicit SATException(const char* msg) noexcept : std::runtime_error(msg) {}

    ~SATException() noexcept override = default;
};

/**
 * \brief Exception that is thrown when an elimination set cannot be computed
 * due to its excessive costs.
 */
class ElimSetException : public std::runtime_error
{
   public:
    explicit ElimSetException() noexcept : std::runtime_error("") {}

    /**
     * \brief Creates a new ElimSetException with a description of its reason
     * \param msg description of the exception's reaon
     */
    explicit ElimSetException(const char* msg) noexcept : std::runtime_error(msg) {}

    ~ElimSetException() noexcept override = default;
};

}  // end namespace hqspre

#endif
