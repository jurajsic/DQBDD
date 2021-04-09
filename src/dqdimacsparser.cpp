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

#include <sstream>
#include <fstream>

// for warnings
#include <iostream>

#include "dqdimacsparser.hpp"
#include "dqbddexceptions.hpp"

namespace dqbdd {

DQDIMACSParser::DQDIMACSParser(Cudd &mgr, QuantifiedVariablesManager &qvmgr) : mgr(mgr), DQBFPrefix(qvmgr) {}

bool DQDIMACSParser::parse(std::string fileName) {
    std::string line;
    std::ifstream inputFile(fileName);

    if (!inputFile.is_open()) {
        std::string errorMes = "Could not open file '";
        errorMes += fileName + "'.";
        throw dqbddException(errorMes);
    }

    bool pLineProcessed = false;
    bool prefixFinished = false;
    unsigned expectedNumOfClauses = 0;
    //int maximumVariable = 0;
    std::string lastToken = "";

    while(std::getline(inputFile, line)) {
        if (line == "") { // TODO whitespaces maybe also can be on empty line???
            // TODO can we have empty lines? if not warning
            continue;
        }
        std::istringstream streamline(line);
        std::string token;
        streamline >> token;
        std::string currentFirstToken = token;
        if (token == "c") {
            continue;
        } else if (token == "p") {
            if (pLineProcessed) {
                throw dqbddException("There are multiple problem lines (lines starting with 'p') in input (DQ)DIMACS.");
            }
            streamline >> token;
            if (token != "cnf") {
                std::cerr << "WARNING: The problem line (i.e. the first line after comments) in input (DQ)DIMACS should have the form 'p cnf <num> <num>'" << std::endl;
            }
            streamline >> token;
            //maximumVariable = std::stoi(token);
            streamline >> token;
            expectedNumOfClauses = std::stoul(token);
            pLineProcessed = true;
        } else if (!pLineProcessed && (token == "a" || token == "e" || token == "d")) {
            throw dqbddException("Input DQDIMACS file is missing the problem line (line starting with 'p').");
        } else if (prefixFinished && (token == "a" || token == "e" || token == "d")) {
            throw dqbddException("Prefix in input DQDIMACS file cannot be between the definition of clauses.");
        } else if (token == "a") {
            if (lastToken == "a") {
                std::cerr << "WARNING: Multiple 'a' lines in input after each other." << std::endl;
            }
            while (streamline >> token) {
                if (token != "0") {
                    Variable univVar(std::stoi(token), mgr);
                    if (DQBFPrefix.isVarExist(univVar)) {
                        throw dqbddException("Cannot have the same variable as both universal and existential.");
                    }
                    DQBFPrefix.addUnivVar(univVar);
                }
            }
        } else if (token == "e") {
            if (lastToken == "e") {
                std::cerr << "WARNING: Multiple 'e' lines in input after each other." << std::endl;
            }
            while (streamline >> token) {
                if (token != "0") {
                    Variable existVar(std::stoi(token), mgr);
                    if (DQBFPrefix.isVarUniv(existVar)) {
                        throw dqbddException("Cannot have the same variable as both universal and existential.");
                    }
                    DQBFPrefix.addExistVar(existVar, DQBFPrefix.getUnivVars());
                }
            }    
        } else if (token == "d") {
            streamline >> token;
            Variable existVar(std::stoi(token), mgr);
            if (DQBFPrefix.isVarUniv(existVar)) {
                throw dqbddException("Cannot have the same variable as both universal and existential.");
            }
            DQBFPrefix.addExistVar(existVar);
            while (streamline >> token) {
                if (token != "0") {
                    Variable univVar(std::stoi(token), mgr);
                    if (!DQBFPrefix.isVarUniv(univVar)) {
                        throw dqbddException("Not able to add existential variable which has non universal variable in dependency list.");
                    }
                    DQBFPrefix.addDependency(existVar, univVar);
                }
            }
        } else { // parse clause (disjunction of literals)
            prefixFinished = true;
            auto getLiteralFromStr = [&](std::string tok) {
                int i = std::stoi(tok);
                if (i < 0) {
                    Variable var = Variable(-i,mgr);
                    if (!DQBFPrefix.isVarHereQuantified(var)) {
                        DQBFPrefix.addExistVar(var);
                    }
                    return Literal(false, var);
                } else {
                    Variable var = Variable(i,mgr);
                    if (!DQBFPrefix.isVarHereQuantified(var)) {
                        DQBFPrefix.addExistVar(var);
                    }
                    return Literal(true, var);
                }
            };
            std::vector<Literal> disj;
            disj.push_back(getLiteralFromStr(token));
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                disj.push_back(getLiteralFromStr(token));
            }
            clauses.push_back(disj);
        }
        lastToken = currentFirstToken;
    }

    if (expectedNumOfClauses != clauses.size()) {
        std::cout << "WARNING: Expected number of clauses is different from the real number of clauses in input DQDIMACS file." << std::endl;
    }
    
    // TODO check if maximumVariable is larger than the number of maximal variable in DQBFPrefix - if not, put warning

    return false;
}

Formula* DQDIMACSParser::getFormula() {
    BDD matrix = mgr.bddOne();
    for (auto &clause : clauses) {
        BDD clauseBDD = mgr.bddZero();
        for (auto &lit : clause) {
            // literal is...
            if (lit.first) { // ...non-negated variable
                clauseBDD |= lit.second;
            } else { // ...negated variable
                clauseBDD |= !lit.second;
            }
        }
        matrix &= clauseBDD;
    }

    Formula *DQBFformula = new Formula(mgr, DQBFPrefix);
    DQBFPrefix.clear();
    DQBFformula->setMatrix(matrix);
    return DQBFformula;
}

QuantifierTreeNode* DQDIMACSParser::getQuantifierTree() {
    std::list<QuantifierTreeNode*> qtClauses;

    if (clauses.size() == 0) {
        auto trueTree = new QuantifierTreeFormula(mgr, DQBFPrefix);
        trueTree->setMatrix(mgr.bddOne());
        DQBFPrefix.clear();
        return trueTree;
    }

    if (clauses.size() == 1) {
        // if we have only one clause, it will be a root
        std::list<QuantifierTreeNode*> literals;
        for (auto &lit : clauses[0]) {
            QuantifierTreeFormula *varFormula;
            if (clauses[0].size() == 1) {
                varFormula = new QuantifierTreeFormula(mgr, DQBFPrefix);
            } else {
                varFormula = new QuantifierTreeFormula(mgr, *DQBFPrefix.getManager());
            }
            // literal is...
            if (lit.first) { // ...non-negated variable
                varFormula->setMatrix(lit.second);
            } else { // ...negated variable
                varFormula->setMatrix(!lit.second);
            }
            literals.push_back(varFormula);
        }
        QuantifierTreeNode *clauseTree;
        if (literals.size() == 0) {
            auto qtzero = new QuantifierTreeFormula(mgr, DQBFPrefix);
            qtzero->setMatrix(mgr.bddZero());
            clauseTree = qtzero;
        } else if (literals.size() == 1) {
            clauseTree = *literals.begin();
        } else {
            clauseTree = new QuantifierTree(false, literals, DQBFPrefix);
        }
        DQBFPrefix.clear();
        return clauseTree;
    }

    for (auto &clause : clauses) {
        std::list<QuantifierTreeNode*> literals;
        for (auto &lit : clause) {
            QuantifierTreeFormula *varFormula = new QuantifierTreeFormula(mgr, *DQBFPrefix.getManager());
            // literal is...
            if (lit.first) { // ...non-negated variable
                varFormula->setMatrix(lit.second);
            } else { // ...negated variable
                varFormula->setMatrix(!lit.second);
            }
            literals.push_back(varFormula);
        }
        QuantifierTreeNode *clauseTree;
        if (literals.size() == 0) {
            auto qtzero = new QuantifierTreeFormula(mgr, *DQBFPrefix.getManager());
            qtzero->setMatrix(mgr.bddZero());
            clauseTree = qtzero;
        } else if (literals.size() == 1) {
            clauseTree = *literals.begin();
        } else {
            clauseTree = new QuantifierTree(false, literals, *DQBFPrefix.getManager());
        }
        qtClauses.push_back(clauseTree);
    }

    auto qt = new QuantifierTree(true, qtClauses, DQBFPrefix);
    DQBFPrefix.clear();
    return qt;
}

} // namespace dqbdd