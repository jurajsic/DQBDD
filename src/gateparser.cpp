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

void GateParser::addExistVar(unsigned long existVarID, bool dependsOnAllDefinedUnivVars) {
    if (dependsOnAllDefinedUnivVars) {
        std::vector<unsigned long> univVarsGateIDs;
        for (const Variable &univVar : DQBFPrefix.getUnivVars()) {
            univVarsGateIDs.push_back(univVar.getId());
        }
        addExistVar(existVarID, univVarsGateIDs);
    } else {
        addExistVar(existVarID, std::vector<unsigned long>{});
    }
}

void GateParser::addExistVar(unsigned long existVarID, const std::vector<unsigned long> &dependencySetVarIDs) {
    Variable varToAdd(existVarID, mgr);

    if (DQBFPrefix.isVarHereQuantified(varToAdd)) {
        throw dqbddException(std::string("Trying to add an existential variable with ID ") + std::to_string(existVarID) + std::string(" in GateParser while it was already added before"));
    } else {
        DQBFPrefix.addExistVar(varToAdd);
        addGate(existVarID, GateType::VAR);
        for (unsigned long univVarID : dependencySetVarIDs) {
            Variable dependentVar(univVarID, mgr);
            if (!DQBFPrefix.isVarUniv(dependentVar)) {
                throw dqbddException(std::string("Trying to add existential variable ") + std::to_string(existVarID) + std::string(" with variable with ID ") + std::to_string(univVarID) + std::string(" which is not defined as universal"));
            } else {
                DQBFPrefix.addDependency(varToAdd, dependentVar);
            }
        }
    }
}

void GateParser::addUnivVar(unsigned long univVarID) {
    Variable varToAdd(univVarID, mgr);

    if (DQBFPrefix.isVarHereQuantified(varToAdd)) {
        throw dqbddException(std::string("Trying to add a universal variable with ID ") + std::to_string(univVarID) + std::string(" in GateParser while it was already added before"));
    } else {
        DQBFPrefix.addUnivVar(varToAdd);
        addGate(univVarID, GateType::VAR);
    }
}

unsigned long GateParser::addGate(GateType type) {
    ++maxGateID;
    addGate(maxGateID, type, std::vector<GateLiteral>{});
    return maxGateID;
}

unsigned long GateParser::addGate(GateType type, const std::vector<GateLiteral> &operands) {
    ++maxGateID;
    addGate(maxGateID, type, operands);
    return maxGateID;
}

void GateParser::addGate(unsigned long gateID, GateType type) {
    addGate(gateID, type, std::vector<GateLiteral>{});
}

void GateParser::addGate(unsigned long gateID, GateType type, const std::vector<GateLiteral> &operands) {
    if (isFormulaParsed) { // we cannot add gates after parsing, we need to first clear parser
        throw dqbddException("Cannot add new gates after parsing is finished");
    }

    if (gateIDToGate.count(gateID) > 0) { // we already added a gate with same ID
        throw dqbddException("You cannot have two gates with the same ID");
    } else {
        for (const GateLiteral &operand : operands) {
            if (gateIDToGate.count(operand.second) == 0) { // we assume that gates which do not yet exist are variable gates
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
                Variable varToAdd(gateID, mgr);
                if (!DQBFPrefix.isVarHereQuantified(varToAdd)) {
                    DQBFPrefix.addExistVar(varToAdd);
                }

                break;
            }
            
            default:
            {
                throw dqbddException("Unsupported gate type during adding a new gate");
                break;
            }
        }

        gateIDToGate[gateID] = {gateID, type, operands};
        gateInputOrder.push_back(gateID);
        if (gateID > maxGateID) {
            maxGateID = gateID;
        }
    }
}

void GateParser::finishedParsing(bool isOutputGatePositive, unsigned long outputGateID) {
    if (gateIDToGate.count(outputGateID) == 0) {
        throw dqbddException("Trying to finish parsing with an output gate ID that does not denote any gate");
    }
    outputGateLiteral = std::make_pair(isOutputGatePositive, outputGateID);
    isFormulaParsed = true;
}

void GateParser::clearParser() {
    DQBFPrefix.clear();
    gateIDToGate.clear();
    gateInputOrder.clear();
    isFormulaParsed = false;
    maxGateID = 0;
    gateIDToNegatedGateID.clear();
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
            auto &processedGate = gateIDToGate[processedGateID];
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

    std::function<QuantifierTreeNode*(const GateLiteral&)> gateLitToTree;
    // isRoot is ugly, but needed to copy variables from DQBFPrefix to the root of tree
    gateLitToTree = [&](const GateLiteral &processedGateLit)->QuantifierTreeNode* {
        QuantifierTreeNode* result;
        auto &processedGateID = processedGateLit.second;
        auto &processedGate = gateIDToGate[processedGateID];
        switch (processedGate.type)
        {
            case GateType::AND:
            {
                if (!processedGateLit.first) {
                    throw dqbddException(std::string("During transformation of AND gate ") + std::to_string(processedGateID) + std::string(" to quantifier tree, there was a negation before it"));
                }
                std::list<QuantifierTreeNode*> treeOperands;
                for (const auto &operandGateLiteral : processedGate.operands) {
                    treeOperands.push_back(gateLitToTree(operandGateLiteral));
                }

                if (treeOperands.size() == 0) { // if AND gate does not have operands, it represents the constant true
                    QuantifierTreeFormula *DQBFtrue = new QuantifierTreeFormula(mgr, *DQBFPrefix.getManager());
                    DQBFtrue->setMatrix(mgr.bddOne());
                    result = DQBFtrue;
                } else if (treeOperands.size() == 1) { // for one operand we just return the tree of this operand
                    result = *treeOperands.begin();
                } else {
                    result = new QuantifierTree(true, treeOperands, *DQBFPrefix.getManager());
                }
                break;
            }

            case GateType::OR:
            {
                if (!processedGateLit.first) {
                    throw dqbddException(std::string("During transformation of OR gate ") + std::to_string(processedGateID) + std::string(" to quantifier tree, there was a negation before it"));
                }
                std::list<QuantifierTreeNode*> treeOperands;
                for (const auto &operandGateLiteral : processedGate.operands) {
                    treeOperands.push_back(gateLitToTree(operandGateLiteral));
                }

                if (treeOperands.size() == 0) { // if OR gate does not have operands, it represents the constant false
                    QuantifierTreeFormula *DQBFfalse = new QuantifierTreeFormula(mgr, *DQBFPrefix.getManager());
                    DQBFfalse->setMatrix(mgr.bddZero());
                    result = DQBFfalse;
                } else if (treeOperands.size() == 1) { // for one operand we just return the tree of this operand
                    result = *treeOperands.begin();
                } else {
                    result = new QuantifierTree(false, treeOperands, *DQBFPrefix.getManager());
                }
                break;
            }

            case GateType::VAR:
            {
                QuantifierTreeFormula *varFormula = new QuantifierTreeFormula(mgr, *DQBFPrefix.getManager());
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

    auto res = gateLitToTree(outputGateLiteral);

    // we move variables from DQBFPrefix to root of quantifier tree (as is done in constructor that takes quantifier manipulator)
    for (Variable uVar : DQBFPrefix.getUnivVars()) {
        res->addUnivVar(uVar);
    }
    for (Variable eVar : DQBFPrefix.getExistVars()) { // dependencies are not needed to move, they are saved in qvMgr of DQBFPrefix
        res->addExistVar(eVar);
    }

    clearParser();
    return res;
}

void GateParser::printPrenexDQCIR(std::ostream &output) {
    // HEADER
    output << "#QCIR-G14" << std::endl;

    printPrefixAndGates(output);
}

void GateParser::printPrenexCleansedDQCIR(std::ostream &output) {
    removeMUXAndXORGates();
    
    // HEADER
    output << "#QCIR-G14 " << maxGateID << std::endl;

    printPrefixAndGates(output);
}

void GateParser::printPrefixAndGates(std::ostream &output) { 
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
        const Gate &processedGate = gateIDToGate[processedGateID];
        switch (processedGate.type)
        {
            case GateType::AND:
            {
                output << processedGateID << " = and(";
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


unsigned long GateParser::addNewGateAtPositionWithoutChecks(GateType type, const std::vector<GateLiteral> &operands, const std::list<unsigned long>::iterator &position) {
    ++maxGateID;
    gateIDToGate[maxGateID] = {maxGateID, type, operands};
    gateInputOrder.insert(position, maxGateID);
    return maxGateID;
}

void GateParser::removeMUXAndXORGates() {
    for (auto gateIter = gateInputOrder.begin(); gateIter != gateInputOrder.end(); ++gateIter) {
        Gate &processedGate = gateIDToGate[*gateIter];
        if (processedGate.type == GateType::MUX) {
            // MUX(A,B,C) = (A AND B) OR (!A AND C)
            const GateLiteral& gateLitA = processedGate.operands[0];
            GateLiteral gateLitNegA = std::make_pair(!gateLitA.first, gateLitA.second);
            const GateLiteral& gateLitB = processedGate.operands[1];
            const GateLiteral& gateLitC = processedGate.operands[2];
            // A AND B - gate is added just before the current gate in gateInputOrder
            unsigned long firstConjGateID = addNewGateAtPositionWithoutChecks(GateType::AND, std::vector<GateLiteral> {gateLitA, gateLitB}, gateIter);
            // !A AND C - gate is added just before the current gate in gateInputOrder
            unsigned long secConjGateID = addNewGateAtPositionWithoutChecks(GateType::AND, std::vector<GateLiteral> {gateLitNegA, gateLitC}, gateIter);
            // change this gate to OR gate with the two new operands
            processedGate.type = GateType::OR;
            processedGate.operands = std::vector<GateLiteral> {std::make_pair(true, firstConjGateID), std::make_pair(true, secConjGateID)};
        } else if (processedGate.type == GateType::XOR) {
            // A XOR B = (A AND !B) OR (!A AND B)
            const GateLiteral& gateLitA = processedGate.operands[0];
            GateLiteral gateLitNegA = std::make_pair(!gateLitA.first, gateLitA.second);
            const GateLiteral& gateLitB = processedGate.operands[1];
            GateLiteral gateLitNegB = std::make_pair(!gateLitB.first, gateLitB.second);
            // A AND !B - gate is added just before the current gate in gateInputOrder
            unsigned long firstConjGateID = addNewGateAtPositionWithoutChecks(GateType::AND, std::vector<GateLiteral> {gateLitA, gateLitNegB}, gateIter);
            // !A AND B - gate is added just before the current gate in gateInputOrder
            unsigned long secConjGateID = addNewGateAtPositionWithoutChecks(GateType::AND, std::vector<GateLiteral> {gateLitNegA, gateLitB}, gateIter);
            // change this gate to OR gate with the two new operands
            processedGate.type = GateType::OR;
            processedGate.operands = std::vector<GateLiteral> {std::make_pair(true, firstConjGateID), std::make_pair(true, secConjGateID)};
        }
    }
}

void GateParser::transformToNNF() {
    removeMUXAndXORGates();
    pushNegation(outputGateLiteral);
}

void GateParser::pushNegation(GateLiteral &gateLiteralToPushNegation) {
    if ((!gateLiteralToPushNegation.first) && gateIDToGate[gateLiteralToPushNegation.second].type != GateType::VAR) {
        gateLiteralToPushNegation = std::make_pair(true, getNegatedGateID(gateLiteralToPushNegation.second));
    }

    for (GateLiteral &operand : gateIDToGate[gateLiteralToPushNegation.second].operands) {
        pushNegation(operand);
    }
}

unsigned long GateParser::getNegatedGateID(unsigned long gateIDToNegate) {
    if (gateIDToNegatedGateID.count(gateIDToNegate) > 0) { // if we already negated given gate
        return gateIDToNegatedGateID[gateIDToNegate];
    } else { // we need to create a negated gate
        const Gate &gateToNegate = gateIDToGate[gateIDToNegate];
        std::vector<GateLiteral> negatedOperands;
        for (const GateLiteral &originalOperand : gateToNegate.operands) {
            negatedOperands.push_back(std::make_pair(!originalOperand.first, originalOperand.second));
        }
        switch(gateToNegate.type) {
            case GateType::AND:
            {
                gateIDToNegatedGateID[gateIDToNegate] = addNewGateAtPositionWithoutChecks(GateType::OR, negatedOperands, gateInputOrder.end());
                break;
            }
            case GateType::OR:
            {
                gateIDToNegatedGateID[gateIDToNegate] = addNewGateAtPositionWithoutChecks(GateType::AND, negatedOperands, gateInputOrder.end());
                break;
            }
            case GateType::VAR:
            {
                // we do not need to negate var gate
                gateIDToNegatedGateID[gateIDToNegate] = gateIDToNegate;
                break;
            }
            default:
            {
                throw dqbddException("We can negate only AND and OR gates");
            }
        }
        return gateIDToNegatedGateID[gateIDToNegate];
    }
}

} // namespace dqbdd