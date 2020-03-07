/*
 * This file is part of DQBDD.
 *
 * Copyright 2020 Juraj Síč
 *
 * DQBDD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * DQBDD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with DQBDD. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DQBF_BDD_SOLVER_HPP
#define DQBF_BDD_SOLVER_HPP

#include <fstream>

#include "cuddObj.hh"

class Solver {
protected:
    Cudd mgr;
public:
    Solver() = delete;
    Solver(const Cudd &mgr);
    virtual ~Solver() = default;
    virtual void readFile(std::ifstream& file) = 0;
    virtual bool solve() = 0;
    virtual void runTests() = 0;
};

#endif