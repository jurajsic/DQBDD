#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <limits>
#include "bdd.h"
#include "Solver.hpp"
#include "BDDPair.hpp"

Solver::Solver(Formula formula) : Solver() {
    this->formula = formula;
}

Solver::Solver(std::ifstream& file) : Solver() {
    readFile(file);
}

void Solver::setFormula(Formula formula) {
    this->formula = formula;
}

bdd Solver::getVarBDDFromStr(std::string strVar) {
    int var_id = std::stoi(strVar);
    if (var_id < 0) {
        return bddProcessor.getBDDReprNeg(Variable(-var_id));
    } else {
        return bddProcessor.getBDDRepr(Variable(var_id));
    }
}

void Solver::readFile(std::ifstream& file) {
    std::string line;
    bdd matrix = bddtrue;

    while(std::getline(file, line)) {
        std::istringstream streamline(line);
        std::string token;
        streamline >> token;
        if (token == "p") {
            streamline >> token; // ignore "cnf"
            streamline >> token; // number of variables
            // TODO decide initial number of variables
            int numOfVariables = std::stoi(token) + 1;
            streamline >> token; // number of CNF conjuncts
            // TODO decide iniatiliazon of BDD based on the size of formula
            bddProcessor.initialize(100000,10000);
            bddProcessor.setNumOfVars(numOfVariables);
        } else if (token == "a") {
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                Variable univVar(std::stoi(token));
                formula.addUnivVar(univVar);
            }
        } else if (token == "e") {
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                Variable existVar(std::stoi(token));
                formula.addExistVar(existVar);
                formula.addDependency(existVar, formula.getUnivVars());
            }    
        } else if (token == "d") {
            streamline >> token;
            Variable existVar(std::stoi(token));
            formula.addExistVar(existVar);
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                Variable univVar(std::stoi(token));
                formula.addDependency(existVar, univVar);
            }
        } else {
            bdd disj = getVarBDDFromStr(token);
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                disj = disj | getVarBDDFromStr(token);
            }
            matrix = matrix & disj;
        }
    }
    formula.setMatrix(matrix);
}

Variable Solver::getSomeUnivVar() {
    int min = std::numeric_limits<int>::max();
    Variable minVar(0);
    for (Variable uVar : formula.getUnivVars()) {
        if (formula.getUnivVarDependencies(uVar).size() < min) {
            minVar = uVar;
        }
    }
    return minVar;
    //return *formula.getUnivVars().begin();
}

bool Solver::solve() {
    while (!formula.getUnivVars().empty()) {
        // remove existential variables that depend on every universal variable
        auto eVars = formula.getExistVars();
        for (Variable eVar : eVars) {
            if (formula.dependsOnEverything(eVar)) {
                formula.removeExistVar(eVar);
                std::cout << "Removing exist variable " << eVar.getId() << std::endl;
                // eliminate from bdd
                bdd newMatrix = bdd_exist(formula.getMatrix(), bddProcessor.getBDDRepr(eVar));
                formula.setMatrix(newMatrix);
                if (formula.isTrue() || formula.isFalse())
                    return formula.isTrue();
            }
        }
        
        // find the universal variable to remove next
        Variable uVarToEliminate = getSomeUnivVar();
        std::cout << "Processing univ variable " << uVarToEliminate.getId() << std::endl;
        auto eVarsToDuplicate = formula.getUnivVarDependencies(uVarToEliminate);
        formula.removeUnivVar(uVarToEliminate);
        
        // pair used for replacing existential variables that depend on uVarToEliminate with new ones
        BDDPair pairToRepl;

        for (Variable eVarToDuplicate : eVarsToDuplicate) {
            //std::cout << "Duplicating var " << eVarToDuplicate.getId() << std::endl;
            Variable newExistVar = bddProcessor.getFreeVariable();
            formula.addExistVar(newExistVar);
            formula.addDependency(newExistVar, formula.getExistVarDependencies(eVarToDuplicate));
            pairToRepl.addToPair(eVarToDuplicate, newExistVar);
            //bddProcessor.addToPairAndSetOrder(pairToRepl, eVarToDuplicate, newExistVar);
        }
        std::cout << "Setting better order with " << formula.getExistVars().size() + formula.getUnivVars().size() << " variables in use" << std::endl;
        //bddProcessor.setNewOrder(pairToRepl);


        std::cout << "Creating BDDs" << std::endl;
        // univ_id=0 where we have old existential variables
        bdd f1 = bdd_restrict(formula.getMatrix(), bddProcessor.getBDDReprNeg(uVarToEliminate));
        std::cout << "Restriction 1 finished" << std::endl;
        // univ_id=1 where we have new existential variables
        bdd f2  = bdd_restrict(formula.getMatrix(), bddProcessor.getBDDRepr(uVarToEliminate));
        std::cout << "Restriction 2 finished" << std::endl;
        f2 = bdd_replace(f2, pairToRepl.getPair());
        std::cout << "Replacing finished" << std::endl;
        // get their conjuction and thus remove univ_id from the formula
        formula.setMatrix(f1 & f2);
        std::cout << "BDD created" << std::endl;

        if (formula.isTrue() || formula.isFalse())
            return formula.isTrue();
    }

    // check if matrix of formula is false
    // if it is -> UNSAT
    // it it is not -> there exists a way in BDD to get to 1 
    //      -> because we only have existential variables left
    //      -> SAT
    return !formula.isFalse();
    
    /*
    for (Variable eVar : formula.getExistVars()) {
        std::cout << "Removing exist variable " << eVar.getId() << std::endl;
        // eliminate from bdd
        bdd newMatrix = bdd_exist(formula.getMatrix(), bddProcessor.getBDDRepr(eVar));
        formula.setMatrix(newMatrix);
        if (formula.isTrue() || formula.isFalse())
            return formula.isTrue();
    }

    return formula.isTrue();
    */
}