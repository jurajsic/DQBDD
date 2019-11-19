#include "quantifiedvariablesmanipulator.hpp"

VariableSet QuantifiedVariablesManipulator::getUnivVars() const {
    return univVars;
}

VariableSet QuantifiedVariablesManipulator::getExistVars() const {
    return existVars;
}

VariableSet QuantifiedVariablesManipulator::getVars() const{
    return vars;
}

void QuantifiedVariablesManipulator::addUnivVar(Variable uVar) {
    auto result = univVars.insert(uVar);
    if (result.second) { // if uVar was not in univVars before
        univVarsDependencies[uVar] = {};
        vars.insert(uVar);
    }
}

void QuantifiedVariablesManipulator::addExistVar(Variable eVar) {
    addExistVar(eVar, VariableSet());
}


void QuantifiedVariablesManipulator::addExistVar(Variable eVar, VariableSet dependencies) {
    auto result = existVars.insert(eVar);
    if (result.second) { // if eVar was not in existVars before
        existVarsDependencies[eVar] = {};
        vars.insert(eVar);
    }
    addDependency(eVar, dependencies);
}

void QuantifiedVariablesManipulator::addDependency(Variable eVar, VariableSet dependencies) {
    for (Variable uVar : dependencies) {
        existVarsDependencies[eVar].insert(uVar);
        univVarsDependencies[uVar].insert(eVar);
    }
}

void QuantifiedVariablesManipulator::addDependency(Variable eVar, Variable uVar) {
    addDependency(eVar, VariableSet{uVar});
}

void QuantifiedVariablesManipulator::removeDependency(Variable eVar, VariableSet dependencies) {
    //existVarsDependencies[eVar].erase(dependencies.begin(),dependencies.end());
    for (Variable remVar : dependencies) {
        existVarsDependencies[eVar].erase(remVar);
        univVarsDependencies[remVar].erase(eVar);
    }
}

void QuantifiedVariablesManipulator::removeDependency(Variable eVar, Variable uVar) {
    removeDependency(eVar, VariableSet{uVar});
}

void QuantifiedVariablesManipulator::removeUnivVar(Variable uVar) {
    auto existVarsToUpdate = univVarsDependencies[uVar];
    for (Variable existVarToUpdate : existVarsToUpdate) {
        removeDependency(existVarToUpdate, uVar);
    }
    univVars.erase(uVar);
    univVarsDependencies.erase(uVar);
    vars.erase(uVar);
}

void QuantifiedVariablesManipulator::removeExistVar(Variable eVar) {
    removeDependency(eVar, existVarsDependencies[eVar]);
    existVars.erase(eVar);
    existVarsDependencies.erase(eVar);
    vars.erase(eVar);
}

void QuantifiedVariablesManipulator::removeVar(Variable var) {
    auto result = vars.erase(var);
    if (result == 0) { // if var was not part of this manipulator
        return;
    }

    if (isUnivVar(var)) { // var is universal var
        removeUnivVar(var);
    } else if (isExistVar(var)) { // var is exist var
        removeExistVar(var);
    }
}

VariableSet QuantifiedVariablesManipulator::getExistVarDependencies(Variable eVar) {
    return existVarsDependencies[eVar];
}

VariableSet QuantifiedVariablesManipulator::getUnivVarDependencies(Variable uVar) {
    return univVarsDependencies[uVar];
}

bool QuantifiedVariablesManipulator::isUnivVar(Variable var) {
    return (univVars.count(var) != 0);
}

bool QuantifiedVariablesManipulator::isExistVar(Variable var) {
    return (existVars.count(var) != 0);
}

bool QuantifiedVariablesManipulator::dependsOnEverything(Variable eVar) {
    return (existVarsDependencies[eVar].size() == allUnivVars.size());
}
