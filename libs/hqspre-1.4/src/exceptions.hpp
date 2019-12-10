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

/**
 * \brief Exception that is thrown when the formula is recognized to be unsatisfiable.
 */
class UNSATException: public std::exception
{
public:
    explicit UNSATException():
        std::exception(),
        message()
    { }

    /**
     * \brief Creates a new UNSATException with a description of its reason
     * \param msg description of the exception's reaon
     */
    explicit UNSATException(const char* msg):
        std::exception(), message(msg)
    { }

    /**
     * \brief Creates a new UNSATException with a description of its reason
     * \param msg description of the exception's reaon
     */
    explicit UNSATException(const std::string& msg):
        std::exception(),
        message(msg)
    { }

    explicit UNSATException(std::string&& msg):
        std::exception(),
        message(std::move(msg))
    { }

    virtual ~UNSATException() noexcept = default;

    /**
     * \brief Returns a description of the reasons why the exception was thrown.
     */
    const char* what() const noexcept override { return message.c_str(); }

private:
    std::string message; ///< Reason for the exception
};


/**
 * \brief Exception that is thrown when the formula is recognized to be satisfiable.
 */
class SATException: public std::exception
{
public:
    explicit SATException():
        std::exception(),
        message()
    { }

    /**
     * \brief Creates a new SATException with a description of its reason
     * \param msg description of the exception's reaon
     */
    explicit SATException(const char * msg):
        std::exception(),
        message(msg)
    { }

    /**
     * \brief Creates a new SATException with a description of its reason
     * \param msg description of the exception's reaon
     */
    explicit SATException(const std::string& msg):
        std::exception(),
        message(msg)
    { }

    /**
     * \brief Creates a new SATException with a description of its reason
     * \param msg description of the exception's reaon
     */
    explicit SATException(std::string&& msg):
        std::exception(),
        message(std::move(msg))
    { }

    virtual ~SATException() noexcept = default;

    /**
     * \brief Returns a description of the reasons why the exception was thrown.
     */
    const char* what() const noexcept override { return message.c_str(); }

private:
    std::string message; ///< Reason for the exception
};


/**
 * \brief Exception that is thrown when an elimination set cannot be computed due to its excessive costs.
 */
class ElimSetException: public std::exception
{
public:
    explicit ElimSetException():
        std::exception(),
        message()
    { }

    /**
     * \brief Creates a new ElimSetException with a description of its reason
     * \param msg description of the exception's reaon
     */
    explicit ElimSetException(const char* msg):
        std::exception(),
        message(msg)
    { }

    /**
     * \brief Creates a new ElimSetException with a description of its reason
     * \param msg description of the exception's reaon
     */
    explicit ElimSetException(const std::string& msg):
        std::exception(),
        message(msg)
    { }

    /**
     * \brief Creates a new ElimSetException with a description of its reason
     * \param msg description of the exception's reaon
     */
    explicit ElimSetException(std::string&& msg):
        std::exception(),
        message(std::move(msg))
    { }

    virtual ~ElimSetException() noexcept = default;

    /**
     * \brief Returns a description of the reasons why the exception was thrown.
     */
    const char* what() const noexcept override { return message.c_str(); }

private:
    std::string message; ///< Reason for the exception
};

} // end namespace hqspre

#endif
