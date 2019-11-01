#include "Formula.hpp"

Formula::Formula(const Cudd &mgr) : mgr(mgr) {

}

VariableSet Formula::getUnivVars() const {
    return univVars;
}

VariableSet Formula::getExistVars() const {
    return existVars;
}

BDD Formula::getMatrix() const {
    return matrix;
}

void Formula::setMatrix(const BDD &matrix) {
    if (mgr.getManager() != matrix.manager()) {
        throw "Managers are fucking different mate";
    }
    this->matrix = matrix;
}

void Formula::addUnivVar(Variable uVar) {
    auto result = univVars.insert(uVar);
    if (result.second) { // if uVar was not in univVars before
        univVarsDependencies[uVar] = {};
    }
}

void Formula::addExistVar(Variable eVar) {
    addExistVar(eVar, VariableSet());
}


void Formula::addExistVar(Variable eVar, VariableSet dependencies) {
    auto result = existVars.insert(eVar);
    if (result.second) { // if eVar was not in existVars before
        existVarsDependencies[eVar] = {};
    }
    addDependency(eVar, dependencies);
}

void Formula::addDependency(Variable eVar, VariableSet dependencies) {
    for (Variable uVar : dependencies) {
        existVarsDependencies[eVar].insert(uVar);
        univVarsDependencies[uVar].insert(eVar);
    }
}

void Formula::addDependency(Variable eVar, Variable uVar) {
    addDependency(eVar, VariableSet{uVar});
}

void Formula::removeDependency(Variable eVar, VariableSet dependencies) {
    //existVarsDependencies[eVar].erase(dependencies.begin(),dependencies.end());
    for (Variable remVar : dependencies) {
        existVarsDependencies[eVar].erase(remVar);
        univVarsDependencies[remVar].erase(eVar);
    }
}

void Formula::removeDependency(Variable eVar, Variable uVar) {
    removeDependency(eVar, VariableSet{uVar});
}

void Formula::removeUnivVar(Variable uVar) {
    auto existVarsToUpdate = univVarsDependencies[uVar];
    for (Variable existVarToUpdate : existVarsToUpdate) {
        removeDependency(existVarToUpdate, uVar);
    }
    univVars.erase(uVar);
    univVarsDependencies.erase(uVar);
}

void Formula::removeExistVar(Variable eVar) {
    removeDependency(eVar, existVarsDependencies[eVar]);
    existVars.erase(eVar);
    existVarsDependencies.erase(eVar);
}

void Formula::removeVar(Variable var) {
    if (univVars.count(var) != 0) { // var is universal var
        removeUnivVar(var);
    } else if (existVars.count(var) != 0) { // var is exist var
        removeExistVar(var);
    }
}

VariableSet Formula::getExistVarDependencies(Variable eVar) {
    return existVarsDependencies[eVar];
}

VariableSet Formula::getUnivVarDependencies(Variable uVar) {
    return univVarsDependencies[uVar];
}


bool Formula::dependsOnEverything(Variable eVar) {
    return (existVarsDependencies[eVar].size() == univVars.size());
}

void Formula::removeUnusedVars() {
    VariableSet usedVars;
    for (unsigned int index : matrix.SupportIndices()) {
        usedVars.insert(Variable(index, mgr));
    }

    VariableSet varsToRemove;
    for (Variable uVar : univVars) {
        if (usedVars.count(uVar) == 0) { // !usedVars.contains(uVar) in c++20
            varsToRemove.insert(uVar);
        }
    }

    for (Variable eVar : existVars) {
        if (usedVars.count(eVar) == 0) { // !usedVars.contains(eVar) in c++20
            varsToRemove.insert(eVar);
        }
    }

    for (Variable varToRemove : varsToRemove) {
        removeVar(varToRemove);
    }
}

/*
bool Formula::isMatrixOne() {
    return matrix == bddtrue;
}

bool Formula::isMatrixZero() {
    return matrix == bddfalse;
}
*/