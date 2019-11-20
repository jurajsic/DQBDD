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

VariableSet QuantifiedVariablesManager::getExistVarDependencies(Variable eVar) {
    return existVarsDependencies[eVar];
}

VariableSet QuantifiedVariablesManager::getUnivVarDependencies(Variable uVar) {
    return univVarsDependencies[uVar];
}



/********************************************/
/********************************************/
/*************** MANIPULATOR ****************/
/********************************************/
/********************************************/

QuantifiedVariablesManipulator::QuantifiedVariablesManipulator() : internalQVManager(), qvMgr(&internalQVManager) {}

QuantifiedVariablesManipulator::QuantifiedVariablesManipulator(QuantifiedVariablesManager &qvMgr) : qvMgr(&qvMgr) {}

VariableSet QuantifiedVariablesManipulator::getUnivVars() const {
    return univVars;
}

VariableSet QuantifiedVariablesManipulator::getExistVars() const {
    return existVars;
}
/*
VariableSet QuantifiedVariablesManipulator::getVars() const{
    return vars;
}
*/
void QuantifiedVariablesManipulator::addUnivVar(Variable uVar) {
    univVars.insert(uVar);
    qvMgr->addUnivVarInstance(uVar);
}

void QuantifiedVariablesManipulator::addExistVar(Variable eVar) {
    addExistVar(eVar, VariableSet());
}


void QuantifiedVariablesManipulator::addExistVar(Variable eVar, VariableSet dependencies) {
    existVars.insert(eVar);
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
    //existVarsDependencies[eVar].erase(dependencies.begin(),dependencies.end());
    qvMgr->removeDependency(eVar, dependencies);
}

void QuantifiedVariablesManipulator::removeDependency(Variable eVar, Variable uVar) {
    removeDependency(eVar, VariableSet{uVar});
}

void QuantifiedVariablesManipulator::removeUnivVar(Variable uVar) {
    if (!isUnivVar(uVar))
        return;
    
    univVars.erase(uVar);
    qvMgr->removeUnivVarInstance(uVar);
}

void QuantifiedVariablesManipulator::removeExistVar(Variable eVar) {
    if (!isExistVar(eVar))
        return;

    existVars.erase(eVar);
    qvMgr->removeExistVarInstance(eVar);
}

void QuantifiedVariablesManipulator::removeVar(Variable var) {
    /*auto result = vars.erase(var);
    if (result == 0) { // if var was not part of this manipulator
        return;
    }*/

    if (isUnivVar(var)) { // var is universal var
        removeUnivVar(var);
    } else if (isExistVar(var)) { // var is exist var
        removeExistVar(var);
    }
}

VariableSet QuantifiedVariablesManipulator::getExistVarDependencies(Variable eVar) {
    if (!isExistVar(eVar))
        throw "This manipualtor does not have this variable as existential";
    return qvMgr->getExistVarDependencies(eVar);
}

VariableSet QuantifiedVariablesManipulator::getUnivVarDependencies(Variable uVar) {
    if (!isUnivVar(uVar))
        throw "This manipualtor does not have this variable as universal";
    return qvMgr->getUnivVarDependencies(uVar);
}

bool QuantifiedVariablesManipulator::isUnivVar(Variable var) {
    return (univVars.contains(var));
}

bool QuantifiedVariablesManipulator::isExistVar(Variable var) {
    return (existVars.contains(var));
}

/*
// TODO delete this because it does not do what you think it does
bool QuantifiedVariablesManipulator::dependsOnEverything(Variable eVar) {
    return (existVarsDependencies[eVar].size() == univVars.size());
}
*/