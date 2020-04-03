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

#ifndef DQBDD_OPTIONS_HPP
#define DQBDD_OPTIONS_HPP

// heuristic by which the next universal variable for elimination is chosen
enum UnivVarElimChoice {
    NumOfDependenciesOnce, // by number of existential variables dependent on u. vars decided at beginning
    NumOfDependenciesContinuous, // by number of existential variables dependent on u. vars decided each time u. var to eliminate is chosen 
    NumOfLeftoverVarsInConjuncts, // by doing  phi[x=0] and phi[x=1] (as BDDs) and checking how much variables are leftover in them
};


class Options {
public:
    // heuristic by which the next universal variable for elimination is chosen
    UnivVarElimChoice uVarElimChoice = UnivVarElimChoice::NumOfDependenciesOnce;
    // remove also universal variables during transforming quantifier tree to formula, 
    bool removeUnivVarsInTree = true;
};

#endif