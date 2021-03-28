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

#include "settings.hpp"
#include <ostream>

namespace hqspre {

static inline const char*
format(bool v)
{
    return (v ? "yes" : "no");
}

std::ostream&
operator<<(std::ostream& stream, const Settings& s)
{
    stream << "[CONFIG] Maximal number of preprocessor loops: " << s.max_loops << '\n'
           << "[CONFIG] Apply universal reduction: " << format(s.univ_reduction) << '\n'
           << "[CONFIG] Apply blocked clause elimination: " << s.bce << '\n'
           << "[CONFIG] Add hidden literals: " << s.hidden << '\n'
           << "[CONFIG] Add covered literals: " << format(s.covered) << '\n'
           << "[CONFIG] Apply blocked literal elimination: " << format(s.ble) << '\n'
           << "[CONFIG] Apply blocked literal addition: " << format(s.bla) << '\n'
           << "[CONFIG] Apply blocked implication addition: " << format(s.bia) << '\n'
           << "[CONFIG] Maximimal hidden/covered clause size: " << s.max_clause_size << '\n'
           << "[CONFIG] Apply hidden subsumption elimination: " << format(s.hse) << '\n'
           << "[CONFIG] Detect hidden equivalences and constants: " << s.hec << '\n'
           << "[CONFIG] Detect implication chains: " << s.impl_chains << '\n'
           << "[CONFIG] Detect contradictions: " << format(s.contradictions) << '\n'
           << "[CONFIG] Use SAT-based semantic gate detection: " << format(s.semantic_gates) << '\n'
           << "[CONFIG] Apply gate substitutions: " << format(s.substitution) << " with maximal costs "
           << s.max_substitution_cost << '\n'
           << "[CONFIG] Apply gate rewriting: " << format(s.rewrite) << '\n'
           << "[CONFIG] Apply self subsumption: " << format(s.self_subsumption) << '\n'
           << "[CONFIG] Apply subsumption: " << format(s.subsumption) << '\n'
           << "[CONFIG] Apply resolution: " << format(s.resolution) << " with maximal costs " << s.max_resolution_cost
           << '\n'
           << "[CONFIG] Apply SAT-based constant checks: " << format(s.sat_const) << '\n'
           << "[CONFIG] Apply SAT-based implication checks: " << format(s.sat_impl) << '\n'
           << "[CONFIG] Apply SAT-based incomplete decision methods: " << format(s.sat_incomplete) << '\n'
           << "[CONFIG] Apply vivification: " << format(s.vivify) << '\n'
           << "[CONFIG] Delete clauses using vivification: " << format(s.vivify_delete) << '\n'
           << "[CONFIG] Number of random patterns for UNSAT checks: " << s.num_random_patterns << '\n'
           << "[CONFIG] Apply universal expansion: " << s.univ_expand << '\n'
           << "[CONFIG] Try to preserve gates: " << s.preserve_gates << '\n'
           << "[CONFIG] Apply consistency checks: " << s.consistency_check << '\n'
#ifdef SKOLEM
           << "[CONFIG] Compute Skolem functions: " << format(s.skolem) << '\n'
#endif
           << std::flush;
    return stream;
}

}  // end namespace hqspre
