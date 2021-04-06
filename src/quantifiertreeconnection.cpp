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


QuantifierTreeConnection::QuantifierTreeConnection(bool isNegated, std::shared_ptr<QuantifierTreeNode> node) : isNegated(isNegated), node(node) {
    // TODO update iterToChildToParent or something
}

QuantifierTreeConnection QuantifierTreeConnection::andQT(bool isNegated, const std::list<QuantifierTreeConnection> &operands, QuantifiedVariablesManager &qvMgr) {
    std::shared_ptr<QuantifierTree> andQuantTree = std::make_shared<QuantifierTree>(QuantifierTreeOperator::andOperator, operands, qvMgr);
    return QuantifierTreeConnection(isNegated, andQuantTree);
}

QuantifierTreeConnection QuantifierTreeConnection::orQT(bool isNegated, const std::list<QuantifierTreeConnection> &operands, QuantifiedVariablesManager &qvMgr) {
    std::shared_ptr<QuantifierTree> orQuantTree = std::make_shared<QuantifierTree>(QuantifierTreeOperator::orOperator, operands, qvMgr);
    return QuantifierTreeConnection(isNegated, orQuantTree);
}
QuantifierTreeConnection QuantifierTreeConnection::varQT(bool isNegated, const Variable &var, const Cudd &mgr, QuantifiedVariablesManager &qvMgr) {
    std::shared_ptr<QuantifierTreeFormula> varFormula = std::make_shared<QuantifierTreeFormula>(mgr, qvMgr);
    varFormula->setMatrix(var);
    return QuantifierTreeConnection(isNegated, varFormula);
}
QuantifierTreeConnection QuantifierTreeConnection::negatedQT(const QuantifierTreeConnection &operand) {
    QuantifierTreeConnection result(operand);
    result.isNegated = !operand.isNegated;
    return result;
}

bool QuantifierTreeConnection::isNodeNegated() const {
    return isNegated;
}

std::shared_ptr<QuantifierTreeNode> QuantifierTreeConnection::getNode() const {
    return node;
}

// void QuantifierTreeConnection::changeNodeToFormula();
// void QuantifierTreeConnection::turnToNNF();
// void QuantifierTreeConnection::localise(const VariableSet &uVarsOutsideThisSubtree);

    //TODO add operations here

    // HMMMM, if we only keep this as children connections and only update the child in the connection, there is no problem, right??? no need to delete from children and we do not keep the set of parents
