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

#include <functional>
#include <vector>

#include "gateparser.hpp"
#include "dqbddexceptions.hpp"

namespace dqbdd {

GateParser::GateParser(Cudd &mgr, QuantifiedVariablesManager &qvmgr) : mgr(mgr), DQBFPrefix(qvmgr) {}


void GateParser::addGate(unsigned long gateID, GateType type) {
    std::vector<GateLiteral> emptyVector;
    addGate(gateID, type, emptyVector);
}

void GateParser::addGate(unsigned long gateID, GateType type, const std::vector<GateLiteral> &operands) {
    if (isFormulaParsed) { // we cannot add gates after parsing, we need to first clear parser
        throw dqbddException("Cannot add new gates after parsing is finished");
    }

    if (gateIDtoGate.count(gateID) > 0) { // we already added a gate with same ID
        throw dqbddException("You cannot have two gates with the same ID");
    } else {
        for (const GateLiteral &operand : operands) {
            if (gateIDtoGate.count(operand.second) == 0) { // we assume that gates which do not yet exist are variable gates
                addGate(operand.second, GateType::VAR);
                //throw dqbddException("Operands of a newly added gate have to already exist");
            }
        }
        switch (type)
        {
            case GateType::AND:
            case GateType::OR:
            {
                // we do not need to check anything here
                break;
            }

            case GateType::MUX:
            {
                if (operands.size() != 3) {
                    throw dqbddException("Wrong number of operands for MUX gate (there should be 3)");
                }
                break;
            }

            case GateType::XOR:
            {
                if (operands.size() != 2) {
                    throw dqbddException("Wrong number of operands for XOR gate (there should be 2)");
                }
                break;
            }

            case GateType::VAR:
            {
                if (operands.size() != 0) {
                    throw dqbddException("Wrong number of operands for VAR gate (there should be none)");
                }
                
                // variables which are not already in prefix are set as existential variables without dependencies
                Variable addedVariable = Variable(gateID, mgr);
                if (!DQBFPrefix.isVarHereQuantified(addedVariable)) {
                    DQBFPrefix.addExistVar(addedVariable);
                }

                break;
            }
            
            default:
            {
                throw dqbddException("Unsupported gate type during adding a new gate");
                break;
            }
        }

        gateIDtoGate[gateID] = {gateID, type, operands};
        gateInputOrder.push_back(gateID);
    }
}

void GateParser::finishedParsing(bool outputGateNegation, unsigned long outputGateID) {
    if (gateIDtoGate.count(outputGateID) == 0) {
        throw dqbddException("Trying to finish parsing with an output gate ID that does not denote any gate");
    }
    outputGateLiteral = std::make_pair(outputGateNegation, outputGateID);
    isFormulaParsed = true;
}

void GateParser::clearParser() {
    DQBFPrefix.clear();
    gateIDtoGate.clear();
    gateInputOrder.clear();
    isFormulaParsed = false;
}

Formula* GateParser::getFormula() {
    if (!isFormulaParsed) {
        throw dqbddException("A file must be parsed first before it is possible to get formula");
    }

    std::unordered_map<unsigned long, BDD> gateIDtoBDD;

    std::function<BDD(const GateLiteral&)> getBDDfromGateLiteral;
    getBDDfromGateLiteral = [&](const GateLiteral &processedGateLit)->BDD { // transforms gate to BDD
        auto &processedGateID = processedGateLit.second;
        BDD result;
        if (gateIDtoBDD.count(processedGateID) > 0) { // we already transformed this gate to BDD
            result = gateIDtoBDD[processedGateID];
        } else { // we need to transform this gate to BDD
            auto &processedGate = gateIDtoGate[processedGateID];
            switch (processedGate.type)
            {
                case GateType::AND:
                {
                    result = mgr.bddOne();
                    for (const auto &operandLiteral : processedGate.operands) {
                        result &= getBDDfromGateLiteral(operandLiteral);
                    }
                    break;
                }

                case GateType::OR:
                {
                    result = mgr.bddZero();
                    for (const auto &operandLiteral : processedGate.operands) {
                        result |= getBDDfromGateLiteral(operandLiteral);
                    }
                    break;
                }

                case GateType::MUX:
                {
                    // TODO check if Ite works properly and if using MUX(A,B,C) = (A AND B) OR (!A AND C) is not better
                    result = getBDDfromGateLiteral(processedGate.operands[0]).Ite(
                                                            getBDDfromGateLiteral(processedGate.operands[1]),
                                                            getBDDfromGateLiteral(processedGate.operands[2]));
                    break;
                }

                case GateType::XOR:
                {
                    result = getBDDfromGateLiteral(processedGate.operands[0]).Xor(getBDDfromGateLiteral(processedGate.operands[1]));
                    break;
                }

                case GateType::VAR:
                {
                    // TODO check if i'm adding variables to DQBFprefix during addGate
                    result = Variable(processedGateID, mgr);
                    break;
                }

                default:
                {
                    // this should not happen
                    throw dqbddException("Unsupported gate type during transformation to formula");
                    break;
                }
            }
            gateIDtoBDD[processedGateID] = result;
        }

        return (processedGateLit.first ? result : !result);
    };

    Formula *DQBFformula = new Formula(mgr, DQBFPrefix);
    DQBFformula->setMatrix(getBDDfromGateLiteral(outputGateLiteral));

    clearParser();
    return DQBFformula;
}

QuantifierTreeNode* GateParser::getQuantifierTree() {
    if (!isFormulaParsed) {
        throw dqbddException("A file must be parsed first before it is possible to get quantifier tree");
    }

    transformToNNF();

    std::function<QuantifierTreeNode*(const GateLiteral&, bool)> gateLitToTree;
    // isRoot is ugly, but needed to copy variables from DQBFPrefix to the root of tree
    gateLitToTree = [&](const GateLiteral &processedGateLit, bool isRoot)->QuantifierTreeNode* {
        QuantifierTreeNode* result;
        auto &processedGateID = processedGateLit.second;
        auto &processedGate = gateIDtoGate[processedGateID];
        switch (processedGate.type)
        {
            case GateType::AND:
            {
                if (!processedGateLit.first) {
                    throw dqbddException("During transformation of AND gate to quantifier tree, there should not be negation before it");
                }
                std::list<QuantifierTreeNode*> treeOperands;
                for (const auto &operandGateLiteral : processedGate.operands) {
                    treeOperands.push_back(gateLitToTree(operandGateLiteral, false));
                }

                if (treeOperands.size() == 0) { // if AND gate does not have operands, it represents the constant true
                    QuantifierTreeFormula *DQBFtrue;
                    if (isRoot) {
                        DQBFtrue = new QuantifierTreeFormula(mgr, DQBFPrefix);
                    } else {
                        DQBFtrue = new QuantifierTreeFormula(mgr, *DQBFPrefix.getManager());
                    }
                    DQBFtrue->setMatrix(mgr.bddOne());
                    result = DQBFtrue;
                } else if (treeOperands.size() == 1) { // for one operand we just return the tree of this operand
                    result = *treeOperands.begin();
                } else {
                    if (isRoot) {
                        result = new QuantifierTree(true, treeOperands, DQBFPrefix);
                    } else {
                        result = new QuantifierTree(true, treeOperands, *DQBFPrefix.getManager());
                    }
                }
                break;
            }

            case GateType::OR:
            {
                if (!processedGateLit.first) {
                    throw dqbddException("During transformation of OR gate to quantifier tree, there should not be negation before it");
                }
                std::list<QuantifierTreeNode*> treeOperands;
                for (const auto &operandGateLiteral : processedGate.operands) {
                    treeOperands.push_back(gateLitToTree(operandGateLiteral, false));
                }

                if (treeOperands.size() == 0) { // if OR gate does not have operands, it represents the constant false
                    QuantifierTreeFormula *DQBFfalse;
                    if (isRoot) {
                        DQBFfalse = new QuantifierTreeFormula(mgr, DQBFPrefix);
                    } else {
                        DQBFfalse = new QuantifierTreeFormula(mgr, *DQBFPrefix.getManager());
                    }
                    DQBFfalse->setMatrix(mgr.bddZero());
                    result = DQBFfalse;
                } else if (treeOperands.size() == 1) { // for one operand we just return the tree of this operand
                    result = *treeOperands.begin();
                } else {
                    if (isRoot) {
                        result = new QuantifierTree(false, treeOperands, DQBFPrefix);
                    } else {
                        result = new QuantifierTree(false, treeOperands, *DQBFPrefix.getManager());
                    }
                }
                break;
            }

            case GateType::VAR:
            {
                QuantifierTreeFormula *varFormula;
                if (isRoot) {
                    varFormula = new QuantifierTreeFormula(mgr, DQBFPrefix);
                } else {
                    varFormula = new QuantifierTreeFormula(mgr, *DQBFPrefix.getManager());
                }
                Variable resultingVar = Variable(processedGateID, mgr);
                varFormula->setMatrix((processedGateLit.first ? resultingVar : !resultingVar));
                result = varFormula;
                break;
            }

            default:
            {
                // this should not happen
                throw dqbddException("An unsupported gate was encountered during transformation to quantifier tree");
                break;
            }
        }
        return result;
    };

    auto res = gateLitToTree(outputGateLiteral, true);
    clearParser();
    return res;
}

void GateParser::printPrenexDQCIR(std::ostream &output) {
    // HEADER
    output << "#QCIR-G14" << std::endl;
    
    // QUANTIFIER PREFIX
    // first universal variables
    if (!DQBFPrefix.getUnivVars().empty()) {
        output << "forall(";
        for (auto uVarIter = DQBFPrefix.getUnivVars().begin(); uVarIter != DQBFPrefix.getUnivVars().end(); ++uVarIter) {
            if (uVarIter == DQBFPrefix.getUnivVars().begin()) {
                output << uVarIter->getId();
            } else {
                output << std::string(", ") << uVarIter->getId();
            }
        }
        output << std::string(")") << std::endl;
    }
    //then existential variables
    for (const Variable &eVar : DQBFPrefix.getExistVars()) {
        output << std::string("depend(") << eVar.getId();
        for (const Variable &depVar : DQBFPrefix.getExistVarDependencies(eVar)) {
            output << std::string(", ") << depVar.getId();
        }
        output << std::string(")") << std::endl;
    }

    // GATES
    // helping function for printing GateLiteral
    auto getGateLiteralString = [&](const GateLiteral &literal)->std::string {
        std::string result = std::to_string(literal.second);
        return (literal.first ? result : (std::string("-") + result));
    };
    // define output gate
    output << "output(" << getGateLiteralString(outputGateLiteral) << ")" << std::endl;
    // we print gates in the order they were added (later added gates depend on earlier ones, we need to therefore print earlier sooner)
    for (unsigned long processedGateID : gateInputOrder) {
        const Gate &processedGate = gateIDtoGate[processedGateID];
        switch (processedGate.type)
        {
            case GateType::AND:
            {
                output << processedGateID << " = and(";
                for (auto operandIter = processedGate.operands.begin(); operandIter != processedGate.operands.end(); ++operandIter) {
                    if (operandIter == processedGate.operands.begin()) {
                        output << getGateLiteralString(*operandIter);
                    } else {
                        output << ", " << getGateLiteralString(*operandIter);
                    }
                }
                output << ")" << std::endl;
                break;
            }

            case GateType::OR:
            {
                output << processedGateID << " = or(";
                break;
            }

            case GateType::MUX:
            {
                output << processedGateID << " = ite(";
                break;
            }

            case GateType::XOR:
            {
                output << processedGateID << " = xor(";
                break;
            }

            case GateType::VAR:
            {
                // we do not print var gates
                break;
            }

            default:
            {
                // this should not happen
                throw dqbddException("An unsupported gate was encountered during printing DQCIR");
                break;
            }
        }

        // print the operands
        if (processedGate.type != GateType::VAR) {
            for (auto operandIter = processedGate.operands.begin(); operandIter != processedGate.operands.end(); ++operandIter) {
                if (operandIter == processedGate.operands.begin()) {
                    output << getGateLiteralString(*operandIter);
                } else {
                    output << ", " << getGateLiteralString(*operandIter);
                }
            }
            output << ")" << std::endl;
        }
    }
}

void GateParser::printPrenexCleansedDQCIR(std::ostream &output) {
    // TODO
}


void GateParser::removeMUXAndXORGates() {
    // TODO
    for (auto gateIter = gateInputOrder.begin(); gateIter != gateInputOrder.end(); ++gateIter) {
        const Gate &processedGate = gateIDtoGate[*gateIter];
        if (processedGate.type == GateType::MUX) {
            // MUX(A,B,C) = (A AND B) OR (!A AND C)
            
        } else if (processedGate.type == GateType::XOR) {
            // A XOR B = (A AND !B) OR (!A AND B)

        }
    }
}

void transformToNNF() {
    // TODO
}

} // namespace dqbdd