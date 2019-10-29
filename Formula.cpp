#include "Formula.hpp"

/*
Formula::Formula(const Formula &f) {
    univVars = f.univVars;
    existVars = f.existVars;
    univVarsDependencies = f.univVarsDependencies;
    existVarsDependencies = f.existVarsDependencies;
    matrix = f.matrix;
}*/

std::unordered_set<Variable> Formula::getUnivVars() const {
    return univVars;
}

std::unordered_set<Variable> Formula::getExistVars() const {
    return existVars;
}

bdd Formula::getMatrix() const {
    return matrix;
}

void Formula::setMatrix(bdd matrix) {
    this->matrix = matrix;
}

void Formula::addUnivVar(Variable uVar) {
    auto result = univVars.insert(uVar);
    if (result.second) { // if uVar was not in univVars before
        univVarsDependencies[uVar] = {};
    }
}

void Formula::addExistVar(Variable eVar) {
    auto result = existVars.insert(eVar);
    if (result.second) { // if eVar was not in existVars before
        existVarsDependencies[eVar] = {};
    }
}

void Formula::addDependency(Variable eVar, std::unordered_set<Variable> dependencies) {
    for (Variable uVar : dependencies) {
        existVarsDependencies[eVar].insert(uVar);
        univVarsDependencies[uVar].insert(eVar);
    }
}

void Formula::addDependency(Variable eVar, Variable uVar) {
    addDependency(eVar, std::unordered_set<Variable>{uVar});
}

void Formula::removeDependency(Variable eVar, std::unordered_set<Variable> dependencies) {
    //existVarsDependencies[eVar].erase(dependencies.begin(),dependencies.end());
    for (Variable remVar : dependencies) {
        existVarsDependencies[eVar].erase(remVar);
        univVarsDependencies[remVar].erase(eVar);
    }
}

void Formula::removeDependency(Variable eVar, Variable uVar) {
    removeDependency(eVar, std::unordered_set<Variable>{uVar});
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

std::unordered_set<Variable> Formula::getExistVarDependencies(Variable eVar) {
    return existVarsDependencies[eVar];
}

std::unordered_set<Variable> Formula::getUnivVarDependencies(Variable uVar) {
    return univVarsDependencies[uVar];
}


bool Formula::dependsOnEverything(Variable eVar) {
    return (existVarsDependencies[eVar].size() == univVars.size());
}

bool Formula::isTrue() {
    return matrix == bddtrue;
}

bool Formula::isFalse() {
    return matrix == bddfalse;
}