/*
 * This file is part of DQBDD.
 *
 * Copyright 2021 Juraj Síč
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
#include "DQBDDexceptions.hpp"


QuantifierTreeConnection::QuantifierTreeConnection(std::shared_ptr<QuantifierTreeNode> node, bool isNegated) {
    // TODO update iterToChildToParent or something
}

QuantifierTreeConnection::QuantifierTreeConnection(Variable var, bool isNegated) {
    // TODO update iterToChildToParent or something
}

QuantifierTreeConnection QuantifierTreeConnection::andQT(bool isNegated, const std::list<QuantifierTreeConnection> &operands, QuantifiedVariablesManager &qvMgr) {
    std::shared_ptr<QuantifierTree> andQuantTree = std::make_shared<QuantifierTree>(true, operands, qvMgr);
    return QuantifierTreeConnection(andQuantTree, isNegated);
}

QuantifierTreeConnection QuantifierTreeConnection::orQT(bool isNegated, const std::list<QuantifierTreeConnection> &operands, QuantifiedVariablesManager &qvMgr) {
    std::shared_ptr<QuantifierTree> orQuantTree = std::make_shared<QuantifierTree>(false, operands, qvMgr);
    return QuantifierTreeConnection(orQuantTree, isNegated);
}
QuantifierTreeConnection QuantifierTreeConnection::varQT(bool isNegated, const Variable &var, const Cudd &mgr, QuantifiedVariablesManager &qvMgr) {
    std::shared_ptr<QuantifierTreeFormula> varFormula = std::make_shared<QuantifierTreeFormula>(mgr, qvMgr);
    return QuantifierTreeConnection(varFormula, isNegated);
}
QuantifierTreeConnection QuantifierTreeConnection::negatedQT(const QuantifierTreeConnection &operand) {
    QuantifierTreeConnection result(operand);
    result.isNegated = !operand.isNegated;
    return result;
}

    //bool isNegated() { return isNegated; }
    //std::shared_ptr<QuantifierTreeNode> getNode() { return node; }

void QuantifierTreeConnection::changeNodeToFormula();
void QuantifierTreeConnection::turnToNNF();
void QuantifierTreeConnection::localise(const VariableSet &uVarsOutsideThisSubtree);

    //TODO add operations here

    // HMMMM, if we only keep this as children connections and only update the child in the connection, there is no problem, right??? no need to delete from children and we do not keep the set of parents
