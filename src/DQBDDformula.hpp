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

#ifndef DQBDD_FORMULA_HPP
#define DQBDD_FORMULA_HPP

#include <unordered_map>
#include <unordered_set>

#include <cuddObj.hh>

#include "DQBDDvariable.hpp"
#include "quantifiedvariablesmanipulator.hpp"

// heuristic by which next universal variable for elimination is chosen
enum UnivVarElimHeuristic {
    NumOfDependenciesOnce, // by number of existential variables dependent on u. vars decided at beginning
    NumOfDependenciesContinuous, // by number of existential variables dependent on u. vars decided each time u. var to eliminate is chosen 
    NumOfLeftoverVarsInConjuncts, // by doing  phi[x=0] and phi[x=1] (as BDDs) and checking how much variables are leftover in them
};

class Formula : public virtual QuantifiedVariablesManipulator {
private:

    Cudd mgr;
    // propositional predicate inside DQBF formula in BDD form
    BDD matrix = mgr.bddZero();

    std::ostream& print(std::ostream& out) const override;

    bool needToRecomputeSupportSet = true;

    UnivVarElimHeuristic uVarElimHeur;
    std::vector<Variable> univVarsOrderToRemove;
    void initializeUnivVarEliminationOrder();
    Variable getUnivVarToEliminate();

public:
    Formula() = delete;
    //Formula(const Cudd &mgr);
    Formula(const Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    Formula(const Cudd &mgr, QuantifiedVariablesManipulator &qvManipulator);
    Formula(const Formula &f) = delete;

    virtual ~Formula() = default;

    VariableSet const &getSupportSet() override;

    BDD getMatrix() const;
    void setMatrix(const BDD &matrix);

    void eliminateUnivVar(Variable uVarToEliminate);

    void eliminateExistVar(Variable existVarToEliminate);
    void eliminateExistVars(VariableSet existVarsToEliminate);
    VariableSet getPossibleExistVarsToEliminate();
    // eliminates all existential variables that are possible to eliminate based on Theorem 5 from DQBF localization paper
    // returns number of eliminated existential variables
    //int eliminatePossibleExistVars();

    // TODO implemenet -> should eliminate all universal variables and all possible exist (can add new exist)
    void eliminatePossibleVars();

    void printStats();
};

#endif