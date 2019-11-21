#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iostream>
#include "solver.hpp"

/*
Solver::Solver(Formula formula) : {
    this->formula = formula;
}

Solver::Solver(std::ifstream& file) : {
    readFile(file);
}
*/

Solver::Solver(const Cudd &mgr) : mgr(mgr), qvMgr(), formula(mgr, qvMgr) {}

void Solver::setFormula(Formula formula) {
    this->formula = Formula(formula);
}

void Solver::setTest1Formula() {
    Formula f(mgr, qvMgr);
    // TODO constructor variable(BDD) was deleted, can be fixed?
    Variable y1(0,mgr);
    Variable y2(1,mgr);
    Variable x1(2,mgr);
    Variable x2(3,mgr);
    f.addExistVar(y1);
    f.addExistVar(y2);
    f.addUnivVar(x1);
    f.addDependency(y1, x1);
    f.addUnivVar(x2);
    f.addDependency(y2, x2);
    BDD m = x1 & x2;
    m = m.Xnor(y1.getBDD().Xnor(y2));
    f.setMatrix(m);
    setFormula(f);
}

void Solver::setTest2Formula() {
    Formula f(mgr, qvMgr);
    // TODO constructor variable(BDD) was deleted, can be fixed?
    Variable y1(0,mgr);
    Variable y2(1,mgr);
    Variable x1(2,mgr);
    Variable x2(3,mgr);
    f.addExistVar(y1);
    f.addExistVar(y2);
    f.addUnivVar(x1);
    f.addDependency(y1, x1);
    f.addUnivVar(x2);
    f.addDependency(y2, x1);
    f.addDependency(y2, x2);
    BDD m = x1 & x2;
    m = m.Xnor(y1.getBDD().Xnor(y2));
    f.setMatrix(m);
    setFormula(f);
}

void Solver::readFile(std::ifstream& file) {
    std::string line;
    BDD matrix = mgr.bddOne();

    while(std::getline(file, line)) {
        std::istringstream streamline(line);
        std::string token;
        streamline >> token;
        if (token == "p") {
            continue;
            // TODO maybe initialize manager here based on the size??
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

Variable Solver::getSomeUnivVar(int choice) {
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

void Solver::printFormulaStats() {
    std::cout << "Formula BDD have " << formula.getMatrix().nodeCount() 
                << " nodes with " << formula.getUnivVars().size() << " universal variables and "
                << formula.getExistVars().size() << " existential variables." << std::endl;
}

void Solver::reorder() {
    std::cout << "Reordering" << std::endl;
    mgr.ReduceHeap();
    printFormulaStats();
}

void Solver::setUnivVarsOrder() {
    univVarsOrderToRemove.assign(formula.getUnivVars().begin(), formula.getUnivVars().end());
    std::sort(univVarsOrderToRemove.begin(), univVarsOrderToRemove.end(),
                [&](Variable a, Variable b) {
                    return (formula.getUnivVarDependencies(a).size() > formula.getUnivVarDependencies(b).size());
                }
            );
}

bool Solver::solve() {
    setUnivVarsOrder();
    while (!formula.getUnivVars().empty()) {
        printFormulaStats();
        
        while (formula.eliminatePossibleExistVars() !=0) { // while something was eliminated
            formula.removeUnusedVars();
        }
        
        if (formula.getUnivVars().empty()) {
            break;
        }
        printFormulaStats();
        reorder();

        // find the universal variable to remove next
        Variable uVarToEliminate = getSomeUnivVar();
        std::cout << "Processing univ variable " << uVarToEliminate.getId() << std::endl;
        formula.eliminateUnivVar(uVarToEliminate);
        
        formula.removeUnusedVars();
    }

    printFormulaStats();

    // check if matrix of formula is 0
    // if it is -> UNSAT
    // it it is not -> there exists a way in BDD to get to 1 
    //                   -> because we only have existential variables left
    //                   -> SAT
    return !formula.getMatrix().IsZero();
}