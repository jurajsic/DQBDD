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

#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iostream>

#include "simplesolver.hpp"

/*
SimpleSolver::SimpleSolver(Formula formula) : {
    this->formula = formula;
}

SimpleSolver::SimpleSolver(std::ifstream& file) : {
    readFile(file);
}
*/

SimpleSolver::SimpleSolver(const Cudd &mgr) : Solver(mgr), qvMgr(), formula(mgr, qvMgr) {}

/*
void SimpleSolver::setFormula(Formula formula) {
    this->formula = Formula(formula);
}
*/

void SimpleSolver::setTest1Formula() {
    formula.clear();
    // TODO constructor variable(BDD) was deleted, can be fixed?
    Variable y1(0,mgr);
    Variable y2(1,mgr);
    Variable x1(2,mgr);
    Variable x2(3,mgr);
    formula.addExistVar(y1);
    formula.addExistVar(y2);
    formula.addUnivVar(x1);
    formula.addDependency(y1, x1);
    formula.addUnivVar(x2);
    formula.addDependency(y2, x2);
    BDD m = x1 & x2;
    m = m.Xnor(y1.getBDD().Xnor(y2));
    formula.setMatrix(m);
}

void SimpleSolver::setTest2Formula() {
    formula.clear();
    // TODO constructor variable(BDD) was deleted, can be fixed?
    Variable y1(0,mgr);
    Variable y2(1,mgr);
    Variable x1(2,mgr);
    Variable x2(3,mgr);
    formula.addExistVar(y1);
    formula.addExistVar(y2);
    formula.addUnivVar(x1);
    formula.addDependency(y1, x1);
    formula.addUnivVar(x2);
    formula.addDependency(y2, x1);
    formula.addDependency(y2, x2);
    BDD m = x1 & x2;
    m = m.Xnor(y1.getBDD().Xnor(y2));
    formula.setMatrix(m);
}

void SimpleSolver::runTests() {
    int numOfCorrect = 0;
    int numOfAll = 0;
    ++numOfAll;
    std::cout << "Test " << numOfAll << std::endl;
    setTest1Formula();
    if (solve() == false) {
        ++numOfCorrect;
    }
    ++numOfAll;
    std::cout << "Test " << numOfAll << std::endl;
    setTest2Formula();
    if (solve() == true) {
        ++numOfCorrect;
    }
    std::cout << numOfCorrect << "/" << numOfAll << std::endl;
}

void SimpleSolver::readFile(std::ifstream& file) {
    std::string line;
    BDD matrix = mgr.bddOne();

    while(std::getline(file, line)) {
        if (line == "") {
            continue;
        }
        std::istringstream streamline(line);
        std::string token;
        streamline >> token;
        if (token == "p" || token == "c") {
            continue;
            // TODO maybe initialize manager here based on the size?? - dont forget to ignore comment then
            /*
            streamline >> token; // ignore "cnf"
            streamline >> token; // number of variables
            // TODO decide initial number of variables
            int numOfVariables = std::stoi(token) + 1;
            streamline >> token; // number of CNF conjuncts
            // TODO decide iniatiliazon of BDD based on the size of formula
            bddProcessor.initialize(100000,10000);
            bddProcessor.setNumOfVars(numOfVariables);
            */
        } else if (token == "a") {
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                Variable univVar(std::stoi(token), mgr);
                formula.addUnivVar(univVar);
            }
        } else if (token == "e") {
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                Variable existVar(std::stoi(token), mgr);
                formula.addExistVar(existVar, formula.getUnivVars());
            }    
        } else if (token == "d") {
            streamline >> token;
            Variable existVar(std::stoi(token), mgr);
            formula.addExistVar(existVar);
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                Variable univVar(std::stoi(token), mgr);
                formula.addDependency(existVar, univVar);
            }
        } else { // parse clause (disjunction of literals)
            // n -> returns BDD variable with index n
            // -n -> returns negated BDD variable with index n
            auto getVarFromStr = [&](std::string tok) {
                int i = std::stoi(tok);
                if (i < 0) {
                    return !mgr.bddVar(-i);
                } else {
                    return mgr.bddVar(i);
                }
            };
            BDD disj = getVarFromStr(token);
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                disj = disj | getVarFromStr(token);
            }
            matrix = matrix & disj;
        }
    }
    formula.setMatrix(matrix);
}

Variable SimpleSolver::getSomeUnivVar(int choice) {
    if (choice == 0) {
        Variable v = univVarsOrderToRemove.back();
        univVarsOrderToRemove.pop_back();
        while (formula.getUnivVars().count(v) == 0) {
            v = univVarsOrderToRemove.back();
            univVarsOrderToRemove.pop_back();
        }
        return v;
    } else if (choice == 1) {
        Variable minVar = *(formula.getUnivVars().begin());
        auto min = formula.getUnivVarDependencies(minVar).size();
        for (Variable uVar : formula.getUnivVars()) {
            if (formula.getUnivVarDependencies(uVar).size() < min) {
                minVar = uVar;
            }
        }
        return minVar;
    } else {
        return *formula.getUnivVars().begin();
    }
}

void SimpleSolver::reorder() {
    std::cout << "Reordering" << std::endl;
    mgr.ReduceHeap();
    formula.printStats();
}

void SimpleSolver::setUnivVarsOrder() {
    univVarsOrderToRemove.assign(formula.getUnivVars().begin(), formula.getUnivVars().end());
    std::sort(univVarsOrderToRemove.begin(), univVarsOrderToRemove.end(),
                [&](Variable a, Variable b) {
                    return (formula.getUnivVarDependencies(a).size() > formula.getUnivVarDependencies(b).size());
                }
            );
}

bool SimpleSolver::solve() {
    //std::cout << formula << std::endl;
    formula.eliminatePossibleVars();
    return (formula.getMatrix().IsOne());

    //TODO delete this:

    setUnivVarsOrder();
    while (!formula.getUnivVars().empty()) {
        //formula.printStats();
        
        VariableSet existVarsToEliminate = formula.getPossibleExistVarsToEliminate();
        while (existVarsToEliminate.size() !=0) {
            formula.eliminateExistVars(existVarsToEliminate);
            formula.removeUnusedVars();
            if (formula.getUnivVars().empty()) {
                break;
            }
            existVarsToEliminate = formula.getPossibleExistVarsToEliminate();
        }

        if (formula.getUnivVars().empty()) {
            break;
        }
        
        //formula.printStats();
        //reorder();

        // find the universal variable to remove next
        Variable uVarToEliminate = getSomeUnivVar();
        std::cout << "Processing univ variable " << uVarToEliminate.getId() << std::endl;
        formula.eliminateUnivVar(uVarToEliminate);
        
        formula.removeUnusedVars();
    }

    //formula.printStats();

    // check if matrix of formula is 0
    // if it is -> UNSAT
    // it it is not -> there exists a way in BDD to get to 1 
    //                   -> because we only have existential variables left
    //                   -> SAT
    return !formula.getMatrix().IsZero();
}