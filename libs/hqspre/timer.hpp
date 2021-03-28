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
 *
 * HQSpre is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with HQSpre. If not, see <http://www.gnu.org/licenses/>.
 */

// $Id: timer.hpp 2352 2019-01-27 19:34:01Z wimmer $

#ifndef HQSPRE_TIMER_HPP_
#define HQSPRE_TIMER_HPP_

#include <iosfwd>

namespace hqspre {

/**
 * \class Timer
 * \brief Class to measure elapsed time.
 * \author Ralf Wimmer, Albert-Ludwigs-University Freiburg, Germany
 */
class Timer
{
   public:
    void   start() noexcept;
    void   stop() noexcept;
    void   reset() noexcept;
    double read() const noexcept;

   private:
    double gettime() const noexcept;

    double _current_time{0.0};  ///< The amount of time that has already passed
    double _start_time{0.0};    ///< The point in time when the timer has been started the last time
    bool   _running{false};     ///< Is the Timer running?
};

/**
 * \brief Starts a given timer upon creation and stops it upon destruction
 */
class ScopeTimer
{
   public:
    explicit ScopeTimer(Timer& timer) noexcept : _timer(timer) { _timer.start(); }

    ~ScopeTimer() noexcept { _timer.stop(); }

    ScopeTimer(const ScopeTimer& other) = delete;
    ScopeTimer(ScopeTimer&& other)      = delete;
    ScopeTimer& operator=(const ScopeTimer& other) = delete;
    ScopeTimer& operator=(ScopeTimer&& other) = delete;

   private:
    Timer& _timer;  ///< The timer to be started and stopped automatically
};

std::ostream& operator<<(std::ostream& stream, const Timer& timer);

#ifdef linux
bool createTimeout(double timeout, void (*timeoutHandler)(int));
#endif

}  // end namespace hqspre

#endif
