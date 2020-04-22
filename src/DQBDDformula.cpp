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

#include <iostream>
#include <algorithm>

#include "DQBDDexceptions.hpp"
#include "DQBDDformula.hpp"

Formula::Formula(const Cudd &mgr, QuantifiedVariablesManager &qvmgr) : QuantifiedVariablesManipulator(qvmgr), mgr(mgr) {}

Formula::Formula(const Cudd &mgr, QuantifiedVariablesManipulator &qvManipulator) : QuantifiedVariablesManipulator(qvManipulator), mgr(mgr) {}

BDD Formula::getMatrix() const {
    return matrix;
}

void Formula::setMatrix(const BDD &matrix) {
    if (mgr.getManager() != matrix.manager()) {
        throw DQBDDexception("Not possible to set matrix of formula with BDD from different CUDD manager.");
    }
    this->matrix = matrix;
    needToRecomputeSupportSet = true;
}

bool Formula::needsToRecomputeSupportSet() const {
    return needToRecomputeSupportSet;
}

VariableSet const &Formula::getSupportSet() {
    if (needToRecomputeSupportSet) {
        supportSet.clear();
        for (unsigned int index : matrix.SupportIndices()) {
            supportSet.insert(Variable(index, mgr));
        }
    }
    needToRecomputeSupportSet = false;
    return supportSet;
}

void Formula::eliminateUnivVar(Variable uVarToEliminate) {
    eliminateUnivVar(uVarToEliminate, false);
}

void Formula::eliminateUnivVar(Variable uVarToEliminate, bool useAlreadyComputedf1f2) {
    if (!getSupportSet().contains(uVarToEliminate)) {
        return;
    }
    
    std::cout << "Eliminating universal variable " << uVarToEliminate.getId() << std::endl;

    // find existential variables that will be duplicated
    VariableSet eVarsToDuplicate;
    VariableSet dependentVars = getUnivVarDependencies(uVarToEliminate);
    for (const Variable &dependentVar : dependentVars) {
        if (getSupportSet().contains(dependentVar)) {
            eVarsToDuplicate.insert(dependentVar);
            removeDependency(dependentVar, uVarToEliminate);
        }
    }

    removeUnivVar(uVarToEliminate);
    
    // pair used for replacing existential variables that depend on uVarToEliminate with new ones
    std::vector<BDD> varsToBeReplaced;
    std::vector<BDD> varsToReplaceWith;

    //std::cout << "Duplicating vars ";
    for (Variable eVarToDuplicate : eVarsToDuplicate) {
        //std::cout << eVarToDuplicate.getId() << " ";
        Variable newExistVar = eVarToDuplicate.newVarAtSameLevel();
        addExistVar(newExistVar, getExistVarDependencies(eVarToDuplicate));
        varsToBeReplaced.push_back(eVarToDuplicate);
        varsToReplaceWith.push_back(newExistVar);
    }
    //std::cout << std::endl;

    //std::cout << "Creating BDDs" << std::endl;
    // uVarToEliminate=false where we have old existential variables
    BDD f1 = (useAlreadyComputedf1f2) // if I already have the restriciton saved
                    ? minf1 : matrix.Restrict(!uVarToEliminate.getBDD());
    //std::cout << "Restriction 1 finished" << std::endl;
    // uVarToEliminate=true where we have new existential variables
    BDD f2 = (useAlreadyComputedf1f2)
                    ? minf2 : matrix.Restrict(uVarToEliminate);
    //std::cout << "Restriction 2 finished" << std::endl;
    f2 = f2.SwapVariables(varsToBeReplaced, varsToReplaceWith);
    //std::cout << "Replacing finished" << std::endl;
    // get their conjuction and thus remove univ_id from the formula
    setMatrix(f1 & f2);
    //std::cout << "BDD created" << std::endl;
}

void Formula::eliminateExistVar(Variable existVarToEliminate) {
    setMatrix(matrix.ExistAbstract(existVarToEliminate));
    removeExistVar(existVarToEliminate);
}

void Formula::eliminateExistVars(VariableSet existVarsToEliminate) {
    if (existVarsToEliminate.empty())
        return;
    
    /* for eliminating multiple existential variables it was just a tiny bit
     * better on benchmarks from QBFEVAL'19 to eliminate them at the same time 
     * using ExistAbstract with the cube of exist vars to eliminate
     */
    //std::cout << "Eliminating exist variables ";
    BDD CubeToRemove = mgr.bddOne();
    for (const Variable &eVarToEliminate : existVarsToEliminate) {
        CubeToRemove = CubeToRemove & eVarToEliminate;
        removeExistVar(eVarToEliminate);
        //std::cout << eVarToEliminate.getId() << " ";
    }
    //std::cout << std::endl;
    setMatrix(matrix.ExistAbstract(CubeToRemove));

    // this is if we do it one by one, it was sometimes a tiny bit slower
    //for (const Variable &eVarToEliminate : existVarsToEliminate) {
    //    eliminateExistVar(eVarToEliminate);
    //}
}

VariableSet Formula::getPossibleExistVarsToEliminate() {
    // the set of univ vars on which exist vars that are possible to eliminate need to depend on
    VariableSet univVarsNeededToDependOn;

    // only exist vars in support set need to be eliminated (assume that other will be removed by removeUnusedVars())
    VariableSet existVarsInSupportSet;

    for (const Variable &var : getSupportSet()) {
        if (isVarUniv(var)) {
            // exist vars possible to eliminate have to depend on all univ vars in support set
            univVarsNeededToDependOn.insert(var);
        } else if (isVarExist(var)) {
            existVarsInSupportSet.insert(var);
        }
    }

    VariableSet possibleExistVarsToEliminate;

    for (const Variable &eVar : existVarsInSupportSet) {
        // check if support set does not contain all universal variables
        if (univVarsNeededToDependOn.size() != qvMgr->getNumberOfUnivVars()) {
            // if it does not there can be some univ vars on which some exist var in support set depends on
            for (const Variable &uVar : getExistVarDependencies(eVar)) {
                // if uVar was not in univVarsNeededToDependOn
                if (univVarsNeededToDependOn.insert(uVar).second) {
                    // all exist vars we went trough already do not depend on this uVar and so they are not possible to eliminate
                    possibleExistVarsToEliminate.clear();
                }
            }
        }
        // only those exist vars that are quantified in this subformula and depend on all univ vars we already went trough can possibly be eliminated
        if (isVarHereQuantified(eVar) && getExistVarDependencies(eVar).size() == univVarsNeededToDependOn.size()) {
            possibleExistVarsToEliminate.insert(eVar);
        }
        
    }

    return possibleExistVarsToEliminate;
}

void Formula::initializeUnivVarEliminationOrder() {
    switch (qvMgr->options.uVarElimChoice)
    {
    case UnivVarElimChoice::NumOfDependenciesOnce:
    {
        univVarsOrderToRemove.assign(getUnivVars().begin(), getUnivVars().end());
        std::sort(univVarsOrderToRemove.begin(), univVarsOrderToRemove.end(),
                    [&](const Variable &a, const Variable &b) {
                        auto aDependenciesSize = getUnivVarDependencies(a).size();
                        auto bDependenciesSize = getUnivVarDependencies(b).size();
                        // variable with more dependencies should be earlier in the vector (so the one with less is chosen sooner)
                        return (aDependenciesSize > bDependenciesSize
                        // and in the case the dependencies have same size, the one with lower id should be eliminated sooner
                                || (aDependenciesSize == bDependenciesSize && a.getId() > b.getId()));
                    }
                );
        break;
    }
    case UnivVarElimChoice::NumOfLeftoverVarsInConjuncts:
    case UnivVarElimChoice::NumOfDependenciesContinuous:
        return;
    default:
        throw DQBDDexception("Chosen heuristic to choose next universal variable to eliminate is not implemented.");
        break;
    }
}

Variable Formula::getUnivVarToEliminate() {
    switch (qvMgr->options.uVarElimChoice)
    {
    case UnivVarElimChoice::NumOfDependenciesOnce:
    {
        Variable v = univVarsOrderToRemove.back();
        univVarsOrderToRemove.pop_back();
        while (!getUnivVars().contains(v)) {
            v = univVarsOrderToRemove.back();
            univVarsOrderToRemove.pop_back();
        }
        return v;
    }
    case UnivVarElimChoice::NumOfDependenciesContinuous:
    {
        return *std::min_element(getUnivVars().begin(), getUnivVars().end(),
                    [&](const Variable &a, const Variable &b) {
                        auto aDependenciesSize = getUnivVarDependencies(a).size();
                        auto bDependenciesSize = getUnivVarDependencies(b).size();
                        return (aDependenciesSize < bDependenciesSize
                                || (aDependenciesSize == bDependenciesSize && a.getId() < b.getId()));
                    }
                );
    }
    case UnivVarElimChoice::NumOfLeftoverVarsInConjuncts:
    {
        auto univVarsIter = getUnivVars().begin();
        Variable minUnivVar = *univVarsIter;
        // minUnivVar=false where we have old existential variables
        minf1 = matrix.Restrict(!minUnivVar);
        // minUnivVar=true where we have new existential variables
        minf2 = matrix.Restrict(minUnivVar);
        VariableSet minf1SupportSet, minf2SupportSet;
        for (unsigned int index : minf1.SupportIndices()) {
            minf1SupportSet.insert(Variable(index, mgr));
        }
        for (unsigned int index : minf2.SupportIndices()) {
            minf2SupportSet.insert(Variable(index, mgr));
        }
        // we take the size of union of variables in both conjucts and add to that..
        auto minNumberOfVariablesInConjuncts = minf1SupportSet.unite(minf2SupportSet.minus(getUnivVarDependencies(minUnivVar))).size() 
                        // ..what would be replaced existential variables in second conjuct
                        + getUnivVarDependencies(minUnivVar).intersect(minf2SupportSet).size();

        ++univVarsIter;
        for (; univVarsIter != getUnivVars().end(); ++univVarsIter) {
            Variable univVar = *univVarsIter;
            // univVar=false where we have old existential variables
            BDD f1 = matrix.Restrict(!univVar);
            // univVar=true where we have new existential variables
            BDD f2 = matrix.Restrict(univVar);

            VariableSet f1SupportSet, f2SupportSet;
            for (unsigned int index : f1.SupportIndices()) {
                f1SupportSet.insert(Variable(index, mgr));
            }
            for (unsigned int index : f2.SupportIndices()) {
                f2SupportSet.insert(Variable(index, mgr));
            }

            // we take the size of union of variables in both conjucts and add to that..
            auto numberOfVariablesInConjuncts = f1SupportSet.unite(f2SupportSet.minus(getUnivVarDependencies(univVar))).size() 
                            // ..what would be replaced existential variables in second conjuct
                            + getUnivVarDependencies(univVar).intersect(f2SupportSet).size();

            if (numberOfVariablesInConjuncts < minNumberOfVariablesInConjuncts
                    || (numberOfVariablesInConjuncts == minNumberOfVariablesInConjuncts
                            && univVar.getId() < minUnivVar.getId())) {
                minUnivVar = univVar;
                minf1 = f1;
                minf2 = f2;
                minNumberOfVariablesInConjuncts = numberOfVariablesInConjuncts;
            }
        }
        return minUnivVar;
    }
    default:
        throw DQBDDexception("Chosen heuristic to choose next universal variable to eliminate is not implemented.");
        break;
    }
}

bool Formula::eliminatePossibleExistVars() {
    bool somethingWasEliminated = false;
    VariableSet existVarsToEliminate = getPossibleExistVarsToEliminate();
    while (existVarsToEliminate.size() !=0) {
        eliminateExistVars(existVarsToEliminate);
        somethingWasEliminated = true;
        existVarsToEliminate = getPossibleExistVarsToEliminate();
    }
    return somethingWasEliminated;
}

bool Formula::eliminateSimpleUnivVars() {
    bool somethingWasEliminated = false;
    // function that returns the universal variables whose existential dependencies are not occuring in this formula
    auto getUnivVarsToEliminate = [&] {
        VariableSet uVarsToEliminate;
        for (const auto &uVar : getUnivVars()) {
            if (getUnivVarDependencies(uVar).intersect(getSupportSet()).empty()) {
                uVarsToEliminate.insert(uVar);
            }
        }
        return uVarsToEliminate;
    };
    auto univVarsToEliminate = getUnivVarsToEliminate();
    while (!univVarsToEliminate.empty()) {
        /*std::cout << "Eliminating universal variables ";
        BDD CubeToRemove = mgr.bddOne();
        for (const Variable &uVarToEliminate : univVarsToEliminate) {
            CubeToRemove = CubeToRemove & uVarToEliminate;
            removeUnivVar(uVarToEliminate);
            std::cout << uVarToEliminate.getId() << " ";
        }
        std::cout << "without dependencies" << std::endl;
        setMatrix(matrix.UnivAbstract(CubeToRemove));
        */
        std::vector<Variable> orderedVarsToEliminate(univVarsToEliminate.begin(), univVarsToEliminate.end()); 
        std::sort(orderedVarsToEliminate.begin(), orderedVarsToEliminate.end(),
                        [](const Variable &a, const Variable &b) 
                                { return a.getId() < b.getId(); });
        for (const Variable &uVarToEliminate : orderedVarsToEliminate) {
            eliminateUnivVar(uVarToEliminate);
        }
        somethingWasEliminated = true;
        removeUnusedVars();
        univVarsToEliminate = getUnivVarsToEliminate();
    }
    return somethingWasEliminated;
}

/**
 * @brief Iteratively remove universal variables and all possible existential variables in this (sub)formula
 * 
 */
void Formula::eliminatePossibleVars() {
    removeUnusedVars();

    initializeUnivVarEliminationOrder();

    while (!getUnivVars().empty()) {
        //printFormulaStats();
        
        VariableSet existVarsToEliminate = getPossibleExistVarsToEliminate();
        while (existVarsToEliminate.size() !=0) {
            eliminateExistVars(existVarsToEliminate);
            removeUnusedVars();
            if (getUnivVars().empty()) {
                break;
            }
            existVarsToEliminate = getPossibleExistVarsToEliminate();
        }
        
        if (getUnivVars().empty()) {
            break;
        }

        //printFormulaStats();
        //reorder();

        Variable uVarToEliminate = getUnivVarToEliminate();
        while (getUnivVarDependencies(uVarToEliminate).size() == 0) {
            // if we have universal variable without dependencies, 
            // we can eliminate multiple universal variables (without dependencies) at once
            eliminateSimpleUnivVars();
            if (getUnivVars().empty()) {
                break;
            }
            uVarToEliminate = getUnivVarToEliminate();
        }

        if (getUnivVars().empty()) {
            break;
        }

        if (qvMgr->options.uVarElimChoice == UnivVarElimChoice::NumOfLeftoverVarsInConjuncts) {
            eliminateUnivVar(uVarToEliminate, true);
        } else {
            eliminateUnivVar(uVarToEliminate, false);
        }
        
        removeUnusedVars();
    }

    // If all exist vars are quantified here -> it means we just need to check if BDD != 0
    // because otherwise there is a path to 1 resulting in assignment which satisfies this formula.
    // But we also have to check if there are no universal vars in support set.
    // If size of support set is larger than the number of quantified exist vars
    // then some univ var that is earlier in formula is in support set.
    if (getExistVars().size() == qvMgr->getNumberOfExistVars() // all existential variables are in this formula --> it is not a subformula but whole formula
         && getSupportSet().size() == getExistVars().size()) { // only existential variables are here
        if (!getMatrix().IsZero()) {
            setMatrix(mgr.bddOne());
            clear();
        }
    } else { // this is just a subformula, we cannot just remove existential variables
        // but we can check if there are some leftover exist vars possible to eliminate
        VariableSet existVarsToEliminate = getPossibleExistVarsToEliminate();
        while (existVarsToEliminate.size() !=0) {
            eliminateExistVars(existVarsToEliminate);
            removeUnusedVars();
            existVarsToEliminate = getPossibleExistVarsToEliminate();
        }
    }

    //printFormulaStats();
}

std::ostream& Formula::print(std::ostream& out) const {
    return out << getMatrix();
}

void Formula::printStats() {
    std::cout << "The BDD has " << matrix.nodeCount() 
                << " nodes and there are " << getUnivVars().size() << " universal and "
                << getExistVars().size() << " existential variables quantified in the formula." << std::endl;
}