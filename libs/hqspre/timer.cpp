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

#include "timer.hpp"

#include <sys/time.h>
#include <csignal>
#include <cstdlib>
#include <ostream>

namespace hqspre {

/**
 * \file timer.cpp
 * \brief Implementation of the Timer class and related functions
 * \author Ralf Wimmer, Albert-Ludwigs-University Freiburg, Germany
 */

/**
   \brief Returns the time used by the current process.
   \note If _user_plus_system_time is set to true, the sum of
     user and system time are returned, otherwise only the user time.
   \return The time used by the current process, measure in CLOCKS_PER_SEC.
*/
static double
gettime() noexcept
{
    timeval start{};
    gettimeofday(&start, nullptr);
    return static_cast<double>(start.tv_sec) + static_cast<double>(start.tv_usec) / 1000000.0;
}

/**
   \brief Starts or continues time measurement.
*/
void
Timer::start() noexcept
{
    if (!_running) {
        _running    = true;
        _start_time = gettime();
    }
}

/**
   \brief Stops time measurement if running.
   Otherwise calling this function has no effect.
*/
void
Timer::stop() noexcept
{
    if (_running) {
        _running = false;
        _current_time += (gettime() - _start_time);
    }
}

/**
   \brief Resets the time to 0. The timer is stopped if it is currently running.
*/
void
Timer::reset() noexcept
{
    _running      = false;
    _current_time = 0.0;
    _start_time   = 0.0;
}

/**
  \brief Returns the measured time in seconds.
*/
double
Timer::read() const noexcept
{
    if (_running) {
        return (_current_time + (gettime() - _start_time));
    } else {
        return _current_time;
    }
}

/**
   \brief Prints the measured time to the given output stream.
   \param stream the stream the measured time is to be printed to
   \param timer the time measurement object whose values is to be printed.
   \note The value of the object is treated as a double value.
*/
std::ostream&
operator<<(std::ostream& stream, const Timer& timer)
{
    stream << timer.read();
    return stream;
}

}  // end namespace hqspre
