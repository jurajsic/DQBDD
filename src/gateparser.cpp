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
#include <unordered_set>

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
    collapseGates();

    std::unordered_map<unsigned long, QuantifierTreeConnection*> gateIDtoQuantifierTreeConnection;
    // for VAR gates, we want to save both non-negated and negated versions, other gates should not have negation in front
    std::unordered_map<unsigned long, QuantifierTreeConnection*> gateIDtoNegatedQuantifierTreeConnection;

    std::function<QuantifierTreeConnection*(const GateLiteral&)> gateLitToTree;
    // isRoot is ugly, but needed to copy variables from DQBFPrefix to the root of tree
    gateLitToTree = [&](const GateLiteral &processedGateLit)->QuantifierTreeConnection* {
        bool isProcessedGatePositive = processedGateLit.first;
        auto &processedGateID = processedGateLit.second;
        auto &processedGate = gateIDToGate[processedGateID];

        if (isProcessedGatePositive && (gateIDtoQuantifierTreeConnection.count(processedGateID) > 0)) { // if gate is not negated and we already processed it before...
            return gateIDtoQuantifierTreeConnection[processedGateID];
        } else if (!isProcessedGatePositive && processedGate.type == GateType::VAR 
                    && gateIDtoNegatedQuantifierTreeConnection.count(processedGateID) > 0) { // ...or if gate is negated var we already processed before....
            return gateIDtoNegatedQuantifierTreeConnection[processedGateID];
        } else { //... otherwise, we have to compute it
            QuantifierTreeConnection* result;
            switch (processedGate.type)
            {
                case GateType::AND:
                {
                    if (!processedGateLit.first) {
                        throw dqbddException(std::string("During transformation of AND gate ") + std::to_string(processedGateID) + std::string(" to quantifier tree, there was a negation before it"));
                    }
                    std::list<QuantifierTreeConnection*> treeOperands;
                    for (const auto &operandGateLiteral : processedGate.operands) {
                        treeOperands.push_back(gateLitToTree(operandGateLiteral));
                    }

                    if (treeOperands.size() == 0) { // if AND gate does not have operands, it represents the constant true
                        QuantifierTreeFormula *DQBFtrue = new QuantifierTreeFormula(mgr, *DQBFPrefix.getManager());
                        DQBFtrue->setMatrix(mgr.bddOne());
                        result = new QuantifierTreeConnection(DQBFtrue);
                    } else if (treeOperands.size() == 1) { // for one operand we just return the tree of this operand
                        result = *treeOperands.begin();
                    } else {
                        result = new QuantifierTreeConnection(new QuantifierTree(true, treeOperands, *DQBFPrefix.getManager()));
                    }
                    break;
                }

                case GateType::OR:
                {
                    if (!processedGateLit.first) {
                        throw dqbddException(std::string("During transformation of OR gate ") + std::to_string(processedGateID) + std::string(" to quantifier tree, there was a negation before it"));
                    }
                    std::list<QuantifierTreeConnection*> treeOperands;
                    for (const auto &operandGateLiteral : processedGate.operands) {
                        treeOperands.push_back(gateLitToTree(operandGateLiteral));
                    }

                    if (treeOperands.size() == 0) { // if OR gate does not have operands, it represents the constant false
                        QuantifierTreeFormula *DQBFfalse = new QuantifierTreeFormula(mgr, *DQBFPrefix.getManager());
                        DQBFfalse->setMatrix(mgr.bddZero());
                        result = new QuantifierTreeConnection(DQBFfalse);
                    } else if (treeOperands.size() == 1) { // for one operand we just return the tree of this operand
                        result = *treeOperands.begin();
                    } else {
                        result = new QuantifierTreeConnection(new QuantifierTree(false, treeOperands, *DQBFPrefix.getManager()));
                    }
                    break;
                }

                case GateType::VAR:
                {
                    QuantifierTreeFormula *varFormula = new QuantifierTreeFormula(mgr, *DQBFPrefix.getManager());
                    Variable resultingVar = Variable(processedGateID, mgr);
                    varFormula->setMatrix((processedGateLit.first ? resultingVar : !resultingVar));
                    result = new QuantifierTreeConnection(varFormula);
                    break;
                }

                default:
                {
                    // this should not happen
                    throw dqbddException("An unsupported gate was encountered during transformation to quantifier tree");
                    break;
                }
            }

            if (isProcessedGatePositive) {
                gateIDtoQuantifierTreeConnection[processedGateID] = result;
            } else {
                gateIDtoNegatedQuantifierTreeConnection[processedGateID] = result;
            }

            return result;
        }
    };

    QuantifierTreeConnection *rootConnection = gateLitToTree(outputGateLiteral);
    QuantifierTreeNode *root = rootConnection->child;

    // we move variables from DQBFPrefix to root of quantifier tree (as is done in constructor that takes quantifier manipulator)
    for (Variable uVar : DQBFPrefix.getUnivVars()) {
        root->addUnivVar(uVar);
    }
    for (Variable eVar : DQBFPrefix.getExistVars()) { // dependencies are not needed to move, they are saved in qvMgr of DQBFPrefix
        root->addExistVar(eVar);
    }

    clearParser();
    
    // we need to delete root connection, as we do not need it and it cannot be used as a child
    delete rootConnection;

    return root;
}

void GateParser::printPrenexDQCIR(std::ostream &output) {
    // HEADER
    output << "#QCIR-G14" << std::endl;

    printPrefix(output);
    printGates(output);
}

void GateParser::printPrenexCleansedDQCIR(std::ostream &output) {
    removeMUXAndXORGates();

    // HEADER
    output << "#QCIR-G14 " << maxGateID << std::endl;

    printPrefix(output);
    printGates(output);
}

void GateParser::printPrefix(std::ostream &output) {
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
}

void GateParser::printGates(std::ostream &output) {
    // helper function for printing GateLiteral
    auto getGateLiteralString = [&](const GateLiteral &literal)->std::string {
        std::string result = std::to_string(literal.second);
        return (literal.first ? result : (std::string("-") + result));
    };
    // define output gate
    output << "output(" << getGateLiteralString(outputGateLiteral) << ")" << std::endl;

    // wasGateIDAlreadyPrinted[gateID] will be true, if we already printed gate to dqcir
    std::unordered_set<unsigned long> alreadyPrintedGateIDs;

    // printGatesRecursivelly(gateID, output) will prints the gate with ID gateID and all gates on which it depends (bottom up)
    std::function<void(unsigned long)> printGatesRecursivelly;
    printGatesRecursivelly = [&](unsigned long gateIDToPrint) {
        if (alreadyPrintedGateIDs.count(gateIDToPrint) > 0) { // if we already printed this gate, we do not do it again
            return;
        } else {
            const Gate &gateToPrint = gateIDToGate[gateIDToPrint];
            // print gates of operands
            for (const GateLiteral &gateToPrintOperand : gateToPrint.operands) {
                printGatesRecursivelly(gateToPrintOperand.second);
            }

            // print the gate
            switch (gateToPrint.type)
            {
                case GateType::AND:
                {
                    output << gateIDToPrint << " = and(";
                    break;
                }

                case GateType::OR:
                {
                    output << gateIDToPrint << " = or(";
                    break;
                }

                case GateType::MUX:
                {
                    output << gateIDToPrint << " = ite(";
                    break;
                }

                case GateType::XOR:
                {
                    output << gateIDToPrint << " = xor(";
                    break;
                }

                case GateType::VAR:
                {
                    // we do not print VAR gates
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
            if (gateToPrint.type != GateType::VAR) {
                for (auto operandIter = gateToPrint.operands.begin(); operandIter != gateToPrint.operands.end(); ++operandIter) {
                    if (operandIter == gateToPrint.operands.begin()) {
                        output << getGateLiteralString(*operandIter);
                    } else {
                        output << ", " << getGateLiteralString(*operandIter);
                    }
                }
                output << ")" << std::endl;
            }

            alreadyPrintedGateIDs.insert(gateIDToPrint);
        }
    };

    printGatesRecursivelly(outputGateLiteral.second);
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

    // the set of gate IDs that are in NNF
    std::unordered_set<unsigned long> gateIDsThatAreInNNF;
    // this function will push the negation (if there is) of the gate literal into its gate and transforms the gate into NNF
    std::function<void(GateLiteral&)> tranformGateLiteralToNNF;
    tranformGateLiteralToNNF = [&](GateLiteral &GateLiteralToTranformToNNF) {
        // the gate which will need to be changed to NNF
        unsigned long gateIDToTranformToNNF = GateLiteralToTranformToNNF.second;

        // if the literal is negated, we push the negation into the gate of the literal (except for VAR gates)
        if ((!GateLiteralToTranformToNNF.first) && gateIDToGate[gateIDToTranformToNNF].type != GateType::VAR) {
            gateIDToTranformToNNF = getNegatedGateID(gateIDToTranformToNNF);
            GateLiteralToTranformToNNF = std::make_pair(true, gateIDToTranformToNNF);
        }

        // if gateIDToTranformToNNF is not in NNF, we transform to NNF its operands (thus making the gate in NNF too)
        if (gateIDsThatAreInNNF.count(gateIDToTranformToNNF) == 0) {
            for (GateLiteral &operand : gateIDToGate[gateIDToTranformToNNF].operands) {
                tranformGateLiteralToNNF(operand);
            }
            gateIDsThatAreInNNF.insert(gateIDToTranformToNNF);
        }
    };

    tranformGateLiteralToNNF(outputGateLiteral);
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
        unsigned long negatedGateID;
        switch(gateToNegate.type) {
            case GateType::AND:
            {
                negatedGateID = addNewGateAtPositionWithoutChecks(GateType::OR, negatedOperands, gateInputOrder.end());
                break;
            }
            case GateType::OR:
            {
                negatedGateID = addNewGateAtPositionWithoutChecks(GateType::AND, negatedOperands, gateInputOrder.end());
                break;
            }
            case GateType::VAR:
            {
                // we do not need to negate var gate
                negatedGateID = gateIDToNegate;
                break;
            }
            default:
            {
                throw dqbddException("We can negate only AND and OR gates");
            }
        }
        gateIDToNegatedGateID[gateIDToNegate] = negatedGateID;
        gateIDToNegatedGateID[negatedGateID] = gateIDToNegate;
        return negatedGateID;
    }
}

void GateParser::collapseGates() {
    // we assume that gates are in NNF, i.e. only AND, OR, and VAR gates + negation only before VAR gates

    std::unordered_set<unsigned long> collapsedGateIDs;
    std::function<void(unsigned long)> collapseGateRecursively;
    collapseGateRecursively = [&](unsigned long gateIDToCollapse) {
        if (collapsedGateIDs.count(gateIDToCollapse) > 0) {
            return;
        } else {
            Gate &gateToCollapse = gateIDToGate[gateIDToCollapse];
            std::vector<GateLiteral> newOperands;
            for (GateLiteral &gateToCollapseOperand : gateToCollapse.operands) {
                collapseGateRecursively(gateToCollapseOperand.second);
                const Gate &operandGate = gateIDToGate[gateToCollapseOperand.second];
                if (gateToCollapse.type == operandGate.type // both are AND or both are OR
                    //&& !gateToCollapse.operands.empty() // and operand does not represent constant value true/false
                ) {  
                    // we collapse here, we add to new operands the operands of operandGate
                    newOperands.insert(newOperands.end(), operandGate.operands.begin(), operandGate.operands.end());
                } else { // we cannot collapse here
                    newOperands.push_back(gateToCollapseOperand);
                }
            }
            gateToCollapse.operands = newOperands;
            collapsedGateIDs.insert(gateIDToCollapse);
        }
    };

    collapseGateRecursively(outputGateLiteral.second);
}

} // namespace dqbdd