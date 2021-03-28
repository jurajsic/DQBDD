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

#ifndef HQSPRE_PROCESS_LIMITS_HPP_
#define HQSPRE_PROCESS_LIMITS_HPP_

#include <cstdint>
#include <ostream>

/**
 * \file process_limits.hpp
 * \author Sven Reimer
 * \date 2016
 * \brief Limiting the number of steps spent in a certain preprocessing
 * operation
 */

namespace hqspre {

enum class PreproMethod : unsigned int
{
    UNIV_REDUCT      = 0,
    UNIV_EXPANSION   = 1,
    EQUIV            = 2,
    BCE              = 3,
    IMPL_CHAINS      = 4,
    CONTRA           = 5,
    HEC              = 6,
    SUBSTIUTION      = 7,
    SUBSUMPTION      = 8,
    SELF_SUBSUMPTION = 9,
    RESOLUTION       = 10,
    UNIT_BY_SAT      = 11,
    IMPL_BY_SAT      = 12,
    UPLA             = 13,
    VIVIFICATION     = 14
};

/**
 * \brief Manages the limits of the differents processes.
 * Works like a "deterministic" timeout
 */
class ProcessLimit
{
   public:
    ProcessLimit() = default;
    ~ProcessLimit() = default;

    void useLimit(bool value) noexcept { _check_limits = value; }

    void increaseLimitBy(std::int64_t value) noexcept { _process_limit += value; }

    void decreaseLimitBy(std::int64_t value) noexcept { _process_limit -= value; }

    void decreaseLimitBy(std::int64_t factor, std::int64_t value) noexcept { _process_limit -= factor * value; }

    void decrementLimit() noexcept { --_process_limit; }

    bool reachedLimit() const noexcept { return (_check_limits && _process_limit < 0); }

    void setLimit(PreproMethod method) { _process_limit = _limits[static_cast<unsigned int>(method)]; }

	ProcessLimit& operator=(const ProcessLimit& other) {
		_process_limit = other._process_limit;
		_check_limits = other._check_limits;
		return *this;
	}

   private:
    const std::int64_t _limits[15] = {
        1ul << 31u,  // universal reduction limit - currently not used
        1ul << 25u,  // universal expansion limit
        1ul << 31u,  // equivalence checking limit - currently not used
        1ul << 29u,  // blocked clause elimination limit
        1ul << 31u,  // implication chain limit - currently not used
        1ul << 28u,  // contradiction limit
        1ul << 31u,  // hidden equivalences and contradiction limit
        1ul << 31u,  // substitution limit - currently not used
        1ul << 31u,  // subsumption limit - currently not used
        1ul << 26u,  // self subsumption limit
        1ul << 25u,  // resolution limit
        1ul << 26u,  // find units by SAT limit
        1ul << 31u,  // find implications by SAT limit
        1ul << 20u,  // upla
        1ul << 15u   // vivification
    };

    /// The given limit
    std::int64_t _process_limit = 0;

    /// Do we check the limits? This may be obsolete (always active) in future
    bool _check_limits = true;

    friend std::ostream& operator<<(std::ostream& stream, const ProcessLimit& pl);
};

inline std::ostream&
operator<<(std::ostream& stream, const ProcessLimit& pl)
{
    stream << pl._process_limit;
    return stream;
}

}  // end namespace hqspre

#endif
