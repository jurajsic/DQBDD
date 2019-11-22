#include <iostream>
#include "formula.hpp"

//Formula::Formula(const Cudd &mgr) : mgr(mgr) {}

Formula::Formula(const Cudd &mgr, QuantifiedVariablesManager &qvmgr) : QuantifiedVariablesManipulator(qvmgr), mgr(mgr) {}

BDD Formula::getMatrix() const {
    return matrix;
}

void Formula::setMatrix(const BDD &matrix) {
    if (mgr.getManager() != matrix.manager()) {
        throw "Managers are fucking different mate";
    }
    this->matrix = matrix;
}

VariableSet Formula::getSupportSet() {
    VariableSet supportSet;
    for (unsigned int index : matrix.SupportIndices()) {
        supportSet.insert(Variable(index, mgr));
    }
    return supportSet;
}

void Formula::eliminateUnivVar(Variable uVarToEliminate) {
    // TODO duplicate only those that are in the bdd????? -> pozri ako som zrobil pushUnivVar ci to fakt treba
    VariableSet eVarsToDuplicate;
    VariableSet supportSet = getSupportSet();
    VariableSet dependentVars = getUnivVarDependencies(uVarToEliminate);
    for (Variable dependentVar : dependentVars) {
        if (supportSet.contains(dependentVar)) {
            eVarsToDuplicate.insert(dependentVar);
            removeDependency(dependentVar, uVarToEliminate);
        }
    }

    removeUnivVar(uVarToEliminate);
    
    // pair used for replacing existential variables that depend on uVarToEliminate with new ones
    std::vector<BDD> varsToBeReplaced;
    std::vector<BDD> varsToReplaceWith;

    std::cout << "Duplicating vars ";
    for (Variable eVarToDuplicate : eVarsToDuplicate) {
        std::cout << eVarToDuplicate.getId() << " ";
        Variable newExistVar = eVarToDuplicate.newVarAtSameLevel();
        addExistVar(newExistVar, getExistVarDependencies(eVarToDuplicate));
        varsToBeReplaced.push_back(eVarToDuplicate);
        varsToReplaceWith.push_back(newExistVar);
    }
    std::cout << std::endl;


    // TODO what is the FUCKING difference between constrain and restrict
    std::cout << "Creating BDDs" << std::endl;
    // uVarToEliminate=false where we have old existential variables
    BDD f1 = matrix.Restrict(!uVarToEliminate.getBDD());
    std::cout << "Restriction 1 finished" << std::endl;
    // uVarToEliminate=true where we have new existential variables
    BDD f2 = matrix.Restrict(uVarToEliminate);
    std::cout << "Restriction 2 finished" << std::endl;
    f2 = f2.SwapVariables(varsToBeReplaced, varsToReplaceWith);
    std::cout << "Replacing finished" << std::endl;
    // get their conjuction and thus remove univ_id from the formula
    setMatrix(f1 & f2);
    std::cout << "BDD created" << std::endl;
}

void Formula::eliminateExistVar(Variable existVarToEliminate) {
    matrix = matrix.ExistAbstract(existVarToEliminate);
    removeExistVar(existVarToEliminate);
}

void Formula::eliminateExistVars(VariableSet existVarsToEliminate) {
    if (existVarsToEliminate.empty())
        return;
    
    std::cout << "Eliminating exist variables ";
    BDD CubeToRemove = mgr.bddOne();
    for (const Variable &eVarToEliminate : existVarsToEliminate) {
        CubeToRemove = CubeToRemove & eVarToEliminate;
        removeExistVar(eVarToEliminate);
        std::cout << eVarToEliminate.getId() << " ";
    }
    std::cout << std::endl;
    matrix = matrix.ExistAbstract(CubeToRemove);
}

// TODO recreate it to getter of possible exist vars to eliminat

VariableSet Formula::getPossibleExistVarsToEliminate() {
    VariableSet supportSet = getSupportSet();
    VariableSet univVarsNeededToDependOn;
    VariableSet possibleExistVarsToEliminate;

    for (const Variable &var : supportSet) {
        if (isVarUniv(var)) {
            if (univVarsNeededToDependOn.insert(var).second) { // if var was not in univVarsNeededToDependOn
                possibleExistVarsToEliminate.clear();
            }
        } else if (isVarExist(var)) {
            for (const Variable &uVar : getExistVarDependencies(var)) {
                if (univVarsNeededToDependOn.insert(uVar).second) { // if uVar was not in univVarsNeededToDependOn
                    possibleExistVarsToEliminate.clear();
                }
            }
            if (isVarHere(var) && getExistVarDependencies(var).size() == univVarsNeededToDependOn.size()) {
                possibleExistVarsToEliminate.insert(var);
            }
        }
    }

    return possibleExistVarsToEliminate;
}

void Formula::eliminatePossibleVars() {
    //setUnivVarsOrder();
    removeUnusedVars();
    while (!getUnivVars().empty()) {
        //printFormulaStats();
        
        VariableSet existVarsToEliminate = getPossibleExistVarsToEliminate();
        while (existVarsToEliminate.size() !=0) {
            eliminateExistVars(existVarsToEliminate);
            removeUnusedVars();
            existVarsToEliminate = getPossibleExistVarsToEliminate();
        }
        
        if (getUnivVars().empty()) {
            break;
        }
        //printFormulaStats();
        //reorder();

        // find the universal variable to remove next
        Variable uVarToEliminate = *getUnivVars().begin();//getSomeUnivVar();
        std::cout << "Processing univ variable " << uVarToEliminate.getId() << std::endl;
        eliminateUnivVar(uVarToEliminate);
        
        removeUnusedVars();
    }

    // TODO add logic for deleting all exist vars (just check if matrix == 0 in that case)

    VariableSet existVarsToEliminate = getPossibleExistVarsToEliminate();
    while (existVarsToEliminate.size() !=0) {
        eliminateExistVars(existVarsToEliminate);
        removeUnusedVars();
        existVarsToEliminate = getPossibleExistVarsToEliminate();
    }

    //printFormulaStats();
}

/*
int Formula::eliminatePossibleExistVars() {
    VariableSet supportSet = getSupportSet();
    VariableSet univVarsNeededToDependOn;
    VariableSet existVarsToEliminate;

    bool isUnivVarInSupportSet = false;
    bool isFreeVarInSupportSet = false;

    for (const Variable &var : supportSet) {
        if (isVarUniv(var)) {
            if (univVarsNeededToDependOn.insert(var).second) { // if var was not in univVarsNeededToDependOn
                existVarsToEliminate.clear();
            }
            isUnivVarInSupportSet = true;
        } else if (isVarExist(var)) {
            for (const Variable &uVar : getExistVarDependencies(var)) {
                if (univVarsNeededToDependOn.insert(uVar).second) { // if uVar was not in univVarsNeededToDependOn
                    existVarsToEliminate.clear();
                }
            }
            if (isVarHere(var)) {
                if (getExistVarDependencies(var).size() == univVarsNeededToDependOn.size()) {
                    existVarsToEliminate.insert(var);
                }
            } else {
                isFreeVarInSupportSet = true;
            }
        } else {
            isFreeVarInSupportSet = true;
        }
    }

    if (!isUnivVarInSupportSet || !isFreeVarInSupportSet) { // if matrix contains only existential variables that are in this formula
        // we can just look if bbd is Zero, if it is not, then eliminating all exist variables would just leave us with One
        if (!matrix.IsZero()) {
            matrix = mgr.bddOne();
        }
        int numberOfExistVars = getExistVars().size();
        clear();
        return numberOfExistVars;
    } else {
        eliminateExistVars(existVarsToEliminate);
        return existVarsToEliminate.size();
    }    
}
*/

/*
bool Formula::isMatrixOne() {
    return matrix == bddtrue;
}

bool Formula::isMatrixZero() {
    return matrix == bddfalse;
}
*/