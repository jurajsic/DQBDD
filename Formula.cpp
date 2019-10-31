#include "Formula.hpp"


VariableSet Formula::getUnivVars() const {
    return univVars;
}

VariableSet Formula::getExistVars() const {
    return existVars;
}

BDD Formula::getMatrix() const {
    return matrix;
}

void Formula::setMatrix(BDD matrix) {
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

VariableSet Formula::getExistVarDependencies(Variable eVar) {
    return existVarsDependencies[eVar];
}

VariableSet Formula::getUnivVarDependencies(Variable uVar) {
    return univVarsDependencies[uVar];
}


bool Formula::dependsOnEverything(Variable eVar) {
    return (existVarsDependencies[eVar].size() == univVars.size());
}

/*
bool Formula::isMatrixOne() {
    return matrix == bddtrue;
}

bool Formula::isMatrixZero() {
    return matrix == bddfalse;
}
*/