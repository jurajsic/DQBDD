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

#ifndef DQBF_BDD_TREESOLVER_HPP
#define DQBF_BDD_TREESOLVER_HPP

#include "solver.hpp"
#include "quantifiertree.hpp"

class TreeSolver : public Solver {
private:
    QuantifiedVariablesManager qvMgr;

    QuantifierTreeNode *root = nullptr;
public:
    TreeSolver(const Cudd &mgr);
    ~TreeSolver();
    void readFile(std::ifstream& file);
    bool solve();
    void setTest1Formula();
    void setTest2Formula();
    void setTest3Formula();
    void runTests();
};

#endif