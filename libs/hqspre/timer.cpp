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
   \brief Returns the time used by the current process.
   \note If _user_plus_system_time is set to true, the sum of
     user and system time are returned, otherwise only the user time.
   \return The time used by the current process, measure in CLOCKS_PER_SEC.
*/
double
Timer::gettime() const noexcept
{
    static timeval start;
    gettimeofday(&start, nullptr);
    return static_cast<double>(start.tv_sec) + static_cast<double>(start.tv_usec) / 1000000.0;
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

#ifdef linux

/**
 * \brief Sets up a timeout after which a signal handler is called
 *
 * @param[in] timeout amount of time in seconds after which the handler is to be
 * called
 * @param[in] timeoutHandler function pointer pointing to the function called
 * after the timeout
 * @return true iff estabilishing the timeout was successful
 */
bool
createTimeout(double timeout, void (*timeoutHandler)(int))
{
    struct sigaction sigact;
    struct sigevent  sigev;
    timer_t          timerid;

    sigact.sa_handler = timeoutHandler;
    sigemptyset(&sigact.sa_mask);
    sigact.sa_flags = 0;
    if (sigaction(SIGUSR1, &sigact, nullptr) != 0) {
        // 'sigaction' failed
        return false;
    }

    sigev.sigev_notify = SIGEV_SIGNAL;
    sigev.sigev_signo  = SIGUSR1;
    if (timer_create(CLOCK_PROCESS_CPUTIME_ID, &sigev, &timerid) != 0) {
        // 'timer_create' failed
        return false;
    }

    struct itimerspec interval;
    interval.it_value.tv_sec = (int)timeout;
    interval.it_value.tv_nsec
        = static_cast<long int>((timeout - static_cast<double>(interval.it_value.tv_sec)) * 1000000000.0);
    interval.it_interval.tv_sec  = 0;
    interval.it_interval.tv_nsec = 0;
    if (timer_settime(timerid, 0, &interval, nullptr) != 0) {
        // 'timer_settime' failed
        return false;
    }

    // everything seems successful
    return true;
}

#endif

}  // end namespace hqspre
