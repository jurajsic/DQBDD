#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
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
            streamline >> token; // ignore cnf
            streamline >> token; // number of variables
            // TODO decide initial number of variables
            int numOfVariables = std::stoi(token);
            streamline >> token; // number of CNF conjuncts
            // TODO decide iniatiliazon of BDD based on the size of formula
            bddProcessor.initialize(1000000,1000000);
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
    return *formula.getUnivVars().begin();
}

bool Solver::solve() {
    while (!formula.getUnivVars().empty()) {
        // find the universal variable to remove next
        Variable uVarToEliminate = getSomeUnivVar();
        auto eVarsToDuplicate = formula.getUnivVarDependencies(uVarToEliminate);
        formula.removeUnivVar(uVarToEliminate);
        
        // pair used for replacing existential variables that depend on uVarToEliminate with new ones
        BDDPair pairToRepl;

        for (Variable eVarToDuplicate : eVarsToDuplicate) {
            Variable newExistVar = bddProcessor.getFreeVariable();
            formula.addExistVar(newExistVar);
            formula.addDependency(newExistVar, formula.getExistVarDependencies(eVarToDuplicate));
            pairToRepl.addToPair(eVarToDuplicate, newExistVar);
        }
        
        // univ_id=0 where we have old existential variables
        bdd f1 = bdd_restrict(formula.getMatrix(), bddProcessor.getBDDReprNeg(uVarToEliminate));
        // univ_id=1 where we have new existential variables
        bdd f2  = bdd_replace(formula.getMatrix(), pairToRepl.getPair());
        f2 = bdd_restrict(f2, bddProcessor.getBDDRepr(uVarToEliminate));
        // get their conjuction and thus remove univ_id from the formula
        formula.setMatrix(f1 & f2);

        auto eVars = formula.getExistVars();
        for (Variable eVar : eVars) {
            if (formula.dependsOnEverything(eVar)) {
                formula.removeExistVar(eVar);
                // eliminate from bdd
                bdd newMatrix = bdd_exist(formula.getMatrix(), bddProcessor.getBDDRepr(eVar));
                formula.setMatrix(newMatrix);
            }
        }
    }
}