#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <limits>
#include <iostream>
//#include "bdd.h"
#include "Solver.hpp"
//#include "BDDPair.hpp"

Solver::Solver(Formula formula) : {
    this->formula = formula;
}

Solver::Solver(std::ifstream& file) : {
    readFile(file);
}

void Solver::setFormula(Formula formula) {
    this->formula = formula;
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
                BDD newMatrix = formula.getMatrix().ExistAbstract(eVar.getRepr());
                if (newMatrix.IsOne() || newMatrix.IsZero()){
                    return newMatrix.IsOne();
                }
                formula.setMatrix(newMatrix);
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
            int eVarToDuplicateLevel = mgr.ReadPerm(eVarToDuplicate.getId());
            Variable newExistVar = mgr.bddNewVarAtLevel(eVarToDuplicateLevel);
            formula.addExistVar(newExistVar);
            formula.addDependency(newExistVar, formula.getExistVarDependencies(eVarToDuplicate));
            pairToRepl.addToPair(eVarToDuplicate, newExistVar);
        }


        // TODO what is the FUCKING difference between constrain and restrict
        std::cout << "Creating BDDs" << std::endl;
        BDD matrix = formula.getMatrix();
        // univ_id=0 where we have old existential variables
        BDD f1 = matrix.Restrict(!uVarToEliminate.getRepr());
        std::cout << "Restriction 1 finished" << std::endl;
        // univ_id=1 where we have new existential variables
        BDD f2  = matrix.Restrict(uVarToEliminate);
        std::cout << "Restriction 2 finished" << std::endl;
        f2 = bdd_replace(f2, pairToRepl.getPair());
        std::cout << "Replacing finished" << std::endl;
        // get their conjuction and thus remove univ_id from the formula
        BDD res = f1 & f2;
        if (res.IsOne() || res.IsOne())
            return res.IsOne();
        formula.setMatrix(res);
        std::cout << "BDD created" << std::endl;
    }

    // check if matrix of formula is 0
    // if it is -> UNSAT
    // it it is not -> there exists a way in BDD to get to 1 
    //                   -> because we only have existential variables left
    //                   -> SAT
    return !formula.getMatrix().IsZero();
    
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