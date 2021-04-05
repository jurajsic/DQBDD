/*
 * This file is part of DQBDD.
 *
 * Copyright 2020, 2021 Juraj Síč
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

#include "quantifiertree.hpp"

QuantifierTreeFormula::QuantifierTreeFormula(const Cudd &mgr, QuantifiedVariablesManager &qvmgr)
                            : QuantifiedVariablesManipulator(qvmgr), QuantifierTreeNode(qvmgr), Formula(mgr, qvmgr) 
                            {
                            }

QuantifierTreeFormula::QuantifierTreeFormula(const Cudd &mgr, QuantifiedVariablesManipulator &qvManipulator)
                            : QuantifiedVariablesManipulator(qvManipulator), QuantifierTreeNode(*qvManipulator.getManager()), Formula(mgr, qvManipulator)
                            {
                            }

void QuantifierTreeFormula::localise(const VariableSet &uVarsOutsideThisSubtree) {
    removeUnusedVars();
}

QuantifierTreeFormula* QuantifierTreeFormula::changeToFormula(Cudd &) {
    return this;
}

void QuantifierTreeFormula::negate() {
    setMatrix(!getMatrix());
}

// computes supportSet and uVarsSupportSet
void QuantifierTreeFormula::computeSupportSets() {
    if (needsToRecomputeSupportSet()) { // only if we need to compute it though
        // this computes supportSet
        Formula::getSupportSet();
        // this computes uVarsSupportSet
        uVarsSupportSet.clear();
        for (const Variable &var : supportSet) {
            if (isVarUniv(var)) {
                uVarsSupportSet.insert(var);
            } else if (isVarExist(var)) {
                uVarsSupportSet.insert(getExistVarDependencies(var).begin(),
                                       getExistVarDependencies(var).end());
            }
        }
    }
}

VariableSet const& QuantifierTreeFormula::getSupportSet() {
    computeSupportSets();
    return supportSet;
}

VariableSet const& QuantifierTreeFormula::getUVarsSupportSet() {
    computeSupportSets();
    return uVarsSupportSet;
}