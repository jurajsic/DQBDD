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

#include <algorithm>

#include "quantifiertree.hpp"

QuantifierTreeNode::QuantifierTreeNode(QuantifiedVariablesManager &qvmgr) : QuantifiedVariablesManipulator(qvmgr) {}

void QuantifierTreeNode::pushExistVar(Variable var) {
    // only put var in this node if it is used in this node
    if (getSupportSet().contains(var)) {
        addExistVar(var);
    }
}

void QuantifierTreeNode::pushUnivVar(Variable var) {
    /* We do not need to rename universal variables,
     * because we assume all existential variables
     * depending on it in this subtree are not in 
     * different subtrees (which is always true if we
     * created this tree by localising, because only 
     * for disjunction we create two subtrees with
     * same existential variables, but for universal
     * variable we cannot push two copies inside the subtrees
     * for disjunction.
     */ 
    if (getSupportSet().contains(var)) { // only put var in this node if it is used in this node
        addUnivVar(var);
    } else { // otherwise we can basically eliminate it by removing dependencies 
        VariableSet dependentExistVars = getUnivVarDependencies(var);
        for (const Variable &dependentExistVar : dependentExistVars) {
            //std::cout << dependentExistVar << " is in support set:" << getSupportSet().contains(dependentExistVar) << std::endl;
            if (getSupportSet().contains(dependentExistVar)) { // we assume this is the only tree that contains dependentExistVar
                removeDependency(dependentExistVar,var);
            }
        }
    }
}

VariableSet const& QuantifierTreeNode::getUVarsSupportSet() {
    return uVarsSupportSet;
}