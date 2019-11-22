#include "quantifiedvariablesmanipulator.hpp"

/********************************************/
/********************************************/
/***************** MANAGER ******************/
/********************************************/
/********************************************/

void QuantifiedVariablesManager::addExistVarInstance(Variable eVar) {
    ++numberOfUsedExistVars[eVar];
    if (numberOfUsedExistVars[eVar] == 1)
        existVarsDependencies[eVar] = {};
}

void QuantifiedVariablesManager::removeExistVarInstance(Variable eVar) {
    --numberOfUsedExistVars[eVar];
    if (numberOfUsedExistVars[eVar] == 0) {
        removeDependency(eVar, existVarsDependencies[eVar]);
        existVarsDependencies.erase(eVar);
    }
}

void QuantifiedVariablesManager::addUnivVarInstance(Variable uVar) {
    ++numberOfUsedUnivVars[uVar];
    if (numberOfUsedUnivVars[uVar] == 1)
        univVarsDependencies[uVar] = {};
    
}

void QuantifiedVariablesManager::removeUnivVarInstance(Variable uVar) {
    --numberOfUsedUnivVars[uVar];
    if (numberOfUsedUnivVars[uVar] == 0) {
        auto existVarsToUpdate = univVarsDependencies[uVar];
        for (Variable existVarToUpdate : existVarsToUpdate) {
            removeDependency(existVarToUpdate, VariableSet{uVar});
        }
        univVarsDependencies.erase(uVar);
    }
}

void QuantifiedVariablesManager::addDependency(Variable eVar, VariableSet dependencies) {
    for (Variable uVar : dependencies) {
        existVarsDependencies[eVar].insert(uVar);
        univVarsDependencies[uVar].insert(eVar);
    }
}

void QuantifiedVariablesManager::removeDependency(Variable eVar, VariableSet dependencies) {
    for (Variable remVar : dependencies) {
        existVarsDependencies[eVar].erase(remVar);
        univVarsDependencies[remVar].erase(eVar);
    }
}

VariableSet const &QuantifiedVariablesManager::getExistVarDependencies(Variable eVar) {
    return existVarsDependencies[eVar];
}

VariableSet const &QuantifiedVariablesManager::getUnivVarDependencies(Variable uVar) {
    return univVarsDependencies[uVar];
}

bool QuantifiedVariablesManager::isVarUniv(Variable var) {
    return (univVarsDependencies.find(var) != univVarsDependencies.end());
}

bool QuantifiedVariablesManager::isVarExist(Variable var) {
    return (existVarsDependencies.find(var) != existVarsDependencies.end());
}


/********************************************/
/********************************************/
/*************** MANIPULATOR ****************/
/********************************************/
/********************************************/

//QuantifiedVariablesManipulator::QuantifiedVariablesManipulator() : internalQVManager(), qvMgr(&internalQVManager) {}

QuantifiedVariablesManipulator::QuantifiedVariablesManipulator(QuantifiedVariablesManager &qvMgr) : qvMgr(&qvMgr) {}

//QuantifiedVariablesManipulator::QuantifiedVariablesManipulator(const QuantifiedVariablesManipulator &qvm) {}

QuantifiedVariablesManipulator::~QuantifiedVariablesManipulator() {
    for (const Variable &uVar : getUnivVars()) {
        qvMgr->removeUnivVarInstance(uVar);
    }

    for (const Variable &eVar : getExistVars()) {
        qvMgr->removeExistVarInstance(eVar);
    }
}

VariableSet const &QuantifiedVariablesManipulator::getUnivVars() const {
    return univVars;
}

VariableSet const &QuantifiedVariablesManipulator::getExistVars() const {
    return existVars;
}
/*
VariableSet QuantifiedVariablesManipulator::getVars() const{
    return vars;
}
*/
void QuantifiedVariablesManipulator::addUnivVar(Variable uVar) {
    if (univVars.insert(uVar).second == true)
        qvMgr->addUnivVarInstance(uVar);
}

void QuantifiedVariablesManipulator::addExistVar(Variable eVar) {
    addExistVar(eVar, VariableSet());
}

void QuantifiedVariablesManipulator::addExistVar(Variable eVar, VariableSet dependencies) {
    if (existVars.insert(eVar).second == true)
        qvMgr->addExistVarInstance(eVar);
    addDependency(eVar, dependencies);
}

void QuantifiedVariablesManipulator::addDependency(Variable eVar, VariableSet dependencies) {
    qvMgr->addDependency(eVar, dependencies);
}

void QuantifiedVariablesManipulator::addDependency(Variable eVar, Variable uVar) {
    addDependency(eVar, VariableSet{uVar});
}

void QuantifiedVariablesManipulator::removeDependency(Variable eVar, VariableSet dependencies) {
    qvMgr->removeDependency(eVar, dependencies);
}

void QuantifiedVariablesManipulator::removeDependency(Variable eVar, Variable uVar) {
    removeDependency(eVar, VariableSet{uVar});
}

void QuantifiedVariablesManipulator::removeUnivVar(Variable uVar) {
    if (univVars.erase(uVar) != 0) { // if uVar was in univVars before erasing
        qvMgr->removeUnivVarInstance(uVar);
    }
}

void QuantifiedVariablesManipulator::removeExistVar(Variable eVar) {
    if (existVars.erase(eVar) != 0) { // if eVar was in existVars before erasing
        qvMgr->removeExistVarInstance(eVar);
    }
}

void QuantifiedVariablesManipulator::removeVar(Variable var) {
    if (!isVarHere(var))
        return;

    if (isVarUniv(var)) { // var is universal var
        removeUnivVar(var);
    } else if (isVarExist(var)) { // var is exist var
        removeExistVar(var);
    }
}

VariableSet const &QuantifiedVariablesManipulator::getExistVarDependencies(Variable eVar) {
    if (!isVarExist(eVar))
        throw "Variable is not existential";
    return qvMgr->getExistVarDependencies(eVar);
}

VariableSet const &QuantifiedVariablesManipulator::getUnivVarDependencies(Variable uVar) {
    if (!isVarUniv(uVar))
        throw "Variable is not universal";
    return qvMgr->getUnivVarDependencies(uVar);
}

bool QuantifiedVariablesManipulator::isVarHere(Variable var) {
    return (univVars.contains(var) || existVars.contains(var));
}

bool QuantifiedVariablesManipulator::isVarUniv(Variable var) {
    return qvMgr->isVarUniv(var);
}

bool QuantifiedVariablesManipulator::isVarExist(Variable var) {
    return qvMgr->isVarExist(var);
}

void QuantifiedVariablesManipulator::clear() {
    for (const Variable &uVar : getUnivVars()) {
        removeUnivVar(uVar);
    }

    for (const Variable &eVar : getExistVars()) {
        removeExistVar(eVar);
    }
}

// removes variables that are not in the support set of matrix (we can do this because 3a) rule )
void QuantifiedVariablesManipulator::removeUnusedVars() {
    VariableSet usedVars = getSupportSet();

    VariableSet varsToRemove;
    for (Variable uVar : getUnivVars()) {
        if (!usedVars.contains(uVar)) {
            varsToRemove.insert(uVar);
        }
    }

    for (Variable eVar : getExistVars()) {
        if (usedVars.contains(eVar)) {
            varsToRemove.insert(eVar);
        }
    }

    for (Variable varToRemove : varsToRemove) {
        removeVar(varToRemove);
    }
}