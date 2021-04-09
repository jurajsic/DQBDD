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

#include <sstream>
#include <fstream>
#include <algorithm>

#include "prenexcleansedqcirparser.hpp"
#include "dqbddexceptions.hpp"

namespace dqbdd {

PrenexCleansedQCIRParser::PrenexCleansedQCIRParser(Cudd &mgr, QuantifiedVariablesManager &qvmgr) : mgr(mgr), DQBFPrefix(qvmgr) {}

PrenexCleansedQCIRParser::Literal PrenexCleansedQCIRParser::getLiteralFromString(std::string LiteralStr) {
    long litNum = std::stol(LiteralStr);
    if (litNum < 0) {
        return Literal(false, -litNum);
    } else {
        return Literal(true, litNum);
    }
}

BDD PrenexCleansedQCIRParser::getBDDFromGate(unsigned long gate) {
    if (gates.count(gate) > 0) {
        OperationAndOperands &ops = gates[gate];
        if (ops.first) { // for conjunction...
            BDD ret = mgr.bddOne();
            for (Literal &operand : ops.second) {
                if (operand.first) {
                    ret &= getBDDFromGate(operand.second);
                } else {
                    ret &= !getBDDFromGate(operand.second);
                }
            }
            return ret;
        } else { //...for disjunction
            BDD ret = mgr.bddZero();
            for (Literal &operand : ops.second) {
                if (operand.first) {
                    ret |= getBDDFromGate(operand.second);
                } else {
                    ret |= !getBDDFromGate(operand.second);
                }
            }
            return ret;
        }
    } else {
        Variable ret = Variable(gate, mgr);
        if (!DQBFPrefix.isVarHereQuantified(ret)) {
            DQBFPrefix.addExistVar(ret);
        }
        return ret;
    }
}

QuantifierTreeNode* PrenexCleansedQCIRParser::getQTFromGate(unsigned long gate) {
    if (gates.count(gate) > 0) {
        OperationAndOperands &ops = gates[gate];
        if (ops.first) { // for conjunction...
            if (ops.second.size() == 0) { // for no operands 'and' gate represents constant true
                auto trueTree = new QuantifierTreeFormula(mgr, DQBFPrefix);
                trueTree->setMatrix(mgr.bddOne());
                return trueTree;
            } else if (ops.second.size() == 1) { 
                // for one operand, we just return the (negation of the) tree of the operand
                auto operandTree = getQTFromGate(ops.second[0].second);
                if (!ops.second[0].first) {
                    operandTree->negate();
                }
                return operandTree;
            } else { // for more operands, we create new tree
                std::list<QuantifierTreeNode*> operands;
                for (Literal &operand : ops.second) {
                    auto operandTree = getQTFromGate(operand.second);
                    if (!operand.first) {
                        operandTree->negate();
                    }
                    operands.push_back(operandTree);
                }
                return (new QuantifierTree(true, operands, *DQBFPrefix.getManager()));
            }
        } else { //...for disjunction
            if (ops.second.size() == 0) { // for no operands 'or' gate represents constant true
                auto trueTree = new QuantifierTreeFormula(mgr, DQBFPrefix);
                trueTree->setMatrix(mgr.bddZero());
                return trueTree;
            } else if (ops.second.size() == 1) { 
                // for one operand, we just return the (negation of the) tree of the operand
                auto operandTree = getQTFromGate(ops.second[0].second);
                if (!ops.second[0].first) {
                    operandTree->negate();
                }
                return operandTree;
            } else { // for more operands, we create new tree
                std::list<QuantifierTreeNode*> operands;
                for (Literal &operand : ops.second) {
                    auto operandTree = getQTFromGate(operand.second);
                    if (!operand.first) {
                        operandTree->negate();
                    }
                    operands.push_back(operandTree);
                }
                return (new QuantifierTree(false, operands, *DQBFPrefix.getManager()));
            }
        }
    } else {
        Variable normalVar = Variable(gate, mgr);
        if (!DQBFPrefix.isVarHereQuantified(normalVar)) {
            DQBFPrefix.addExistVar(normalVar);
        }
        QuantifierTreeFormula *varFormula = new QuantifierTreeFormula(mgr, *DQBFPrefix.getManager());
        varFormula->setMatrix(normalVar);
        return varFormula;
    }
}

bool PrenexCleansedQCIRParser::parse(std::string fileName) {
    std::string line;
    std::ifstream inputFile(fileName);

    if (!inputFile.is_open()) {
        std::string errorMes = "Could not open file '";
        errorMes += fileName + "'.";
        throw dqbddException(errorMes);
    }

    bool prefixFinished = false;

    while(std::getline(inputFile, line)) {
        if (line == "" || line[0] == '#') {
            continue;
        }

        std::replace_if(line.begin(), line.end(), [](char c) {return (c == '(' || c == ')' || c == ',');}, ' ');
        std::istringstream streamLine(line);
        std::string token;
        streamLine >> token;

        if (!prefixFinished) {
            // processing quantifier prefix
            if (token == "exists") {
                while (streamLine >> token) {
                    Variable existVar(std::stoul(token), mgr);
                    if (DQBFPrefix.isVarUniv(existVar)) {
                        throw dqbddException("Cannot have the same variable as both universal and existential.");
                    }
                    DQBFPrefix.addExistVar(existVar, DQBFPrefix.getUnivVars());
                }
            } else if (token == "forall") {
                while (streamLine >> token) {
                    Variable univVar(std::stoul(token), mgr);
                    if (DQBFPrefix.isVarExist(univVar)) {
                        throw dqbddException("Cannot have the same variable as both universal and existential.");
                    }
                    DQBFPrefix.addUnivVar(univVar);
                }
            } else if (token == "depend") {
                streamLine >> token;
                Variable existVar(std::stoul(token), mgr);
                if (DQBFPrefix.isVarUniv(existVar)) {
                    throw dqbddException("Cannot have the same variable as both universal and existential.");
                }
                DQBFPrefix.addExistVar(existVar);
                while (streamLine >> token) {
                    Variable univVar(std::stoul(token), mgr);
                    if (!DQBFPrefix.isVarUniv(univVar)) {
                        throw dqbddException("Not able to add existential variable which has non universal variable in dependency list.");
                    }
                    DQBFPrefix.addDependency(existVar, univVar);
                }
            } else if (token == "output") {
                streamLine >> token;
                outputGate = getLiteralFromString(token);
                prefixFinished = true;
            } else {
                throw dqbddException("Unexpected token found in the quantifier prefix of the input file (maybe forgotten output gate?).");
            }
        } else {
            // processing gates
            unsigned long inputGate = std::stoul(token);
            
            streamLine >> token;
            if (token != "=") {
                throw dqbddException("Unexpected token in input file");
            }
            
            OperationAndOperands ops;
            
            streamLine >> token;
            if (token == "and") {
                ops.first = true;
            } else if (token == "or") {
                ops.first = false;
            } else {
                throw dqbddException("Only operations 'and' and 'or' are allowed in cleansed QCIR");
            }
            
            while (streamLine >> token) {
                ops.second.push_back(getLiteralFromString(token));
            }

            if (gates.count(inputGate) > 0) {
                throw dqbddException("There cannot be two definitions of the same gate");
            }
            gates[inputGate] = ops;
        }
    }

    return false;
}

Formula* PrenexCleansedQCIRParser::getFormula() {
    BDD matrix;
    if (outputGate.first) {
        matrix = getBDDFromGate(outputGate.second); 
    } else {
        matrix = !getBDDFromGate(outputGate.second);
    }

    Formula *DQBFformula = new Formula(mgr, DQBFPrefix);
    DQBFPrefix.clear();
    DQBFformula->setMatrix(matrix);
    return DQBFformula;
}

QuantifierTreeNode* PrenexCleansedQCIRParser::getQuantifierTree() {
    auto outputGateTree = getQTFromGate(outputGate.second);
    if (!outputGate.first) {
        outputGateTree->negate();
    }
    for (Variable uVar : DQBFPrefix.getUnivVars()) {
        outputGateTree->addUnivVar(uVar);
    }

    for (Variable eVar : DQBFPrefix.getExistVars()) {
        outputGateTree->addExistVar(eVar);
    }
    DQBFPrefix.clear();
    return outputGateTree;
}

} // namespace dqbdd