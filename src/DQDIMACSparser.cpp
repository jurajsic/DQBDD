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

#include <sstream>
#include <fstream>

#include "DQDIMACSparser.hpp"
#include "DQBDDexceptions.hpp"

DQDIMACSParser::DQDIMACSParser(Cudd &mgr, QuantifiedVariablesManager &qvmgr) : mgr(mgr), DQBFPrefix(qvmgr) {}

bool DQDIMACSParser::parse(std::string fileName) {
    std::string line;
    std::ifstream inputFile(fileName);

    while(std::getline(inputFile, line)) {
        if (line == "") {
            continue;
        }
        std::istringstream streamline(line);
        std::string token;
        streamline >> token;
        if (token == "p" || token == "c") {
            continue;
        } else if (token == "a") {
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                Variable univVar(std::stoi(token), mgr);
                DQBFPrefix.addUnivVar(univVar);
            }
        } else if (token == "e") {
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                Variable existVar(std::stoi(token), mgr);
                DQBFPrefix.addExistVar(existVar, DQBFPrefix.getUnivVars());
            }    
        } else if (token == "d") {
            streamline >> token;
            Variable existVar(std::stoi(token), mgr);
            DQBFPrefix.addExistVar(existVar);
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                Variable univVar(std::stoi(token), mgr);
                if (!DQBFPrefix.isVarUniv(univVar)) {
                    throw DQBDDexception("Not able to add existential variable which has non universal variable in dependency list.");
                }
                DQBFPrefix.addDependency(existVar, univVar);
            }
        } else { // parse clause (disjunction of literals)
            auto getLiteralFromStr = [&](std::string tok) {
                int i = std::stoi(tok);
                if (i < 0) {
                    return Literal(false, Variable(i,mgr));
                } else {
                    return Literal(true, Variable(i,mgr));
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
    }
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
        }
        QuantifierTree *clauseTree = new QuantifierTree(false, literals, *DQBFPrefix.getManager());
        qtClauses.push_back(clauseTree);
    }

    auto qt = new QuantifierTree(true, qtClauses, DQBFPrefix);
    DQBFPrefix.clear();
    return qt;
}