/*
 * This file is part of DQBDD.
 *
 * Copyright 2020 Juraj Síč
 *
 * DQBDD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * DQBDD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with DQBDD. If not, see <http://www.gnu.org/licenses/>.
 */

#include "quantifiedvariablesmanipulator.hpp"
#include "DQBDDexceptions.hpp"

bool VariableSet::contains(Variable const &var) const {
    return (this->count(var) != 0);
}

VariableSet VariableSet::intersect(const VariableSet &vs) const {
    VariableSet intersection = { };
    for (auto iter = this->begin(); iter != this->end(); ++iter) {
        if (vs.contains(*iter)) {
            intersection.insert(*iter);
        }
    }
    return intersection;
}

VariableSet VariableSet::unite(const VariableSet &vs) const {
    VariableSet unionn = vs;
    unionn.insert(this->begin(),this->end());
    return unionn;
}

VariableSet VariableSet::minus(const VariableSet &vs) const {
    VariableSet minus = {};
    for (auto iter = this->begin(); iter != this->end(); ++iter) {
        if (!vs.contains(*iter)) {
            minus.insert(*iter);
        }
    }
    return minus;
}

VariableSet VariableSet::getSupportSetOfBDD(const BDD &b, Cudd &mgr) {
    VariableSet bSupportSet;
    for (unsigned int index : b.SupportIndices()) {
        bSupportSet.insert(Variable(index, mgr));
    }
    return bSupportSet;
}

std::ostream& operator<<(std::ostream& os, const VariableSet& variableSet) {
    os << std::string("{");
    for (auto iter = variableSet.begin(); iter != variableSet.end(); ++iter) {
        if (iter == variableSet.begin()) {
            os << *iter;
        } else {
            os << std::string(", ") << *iter;
        }
    }
    os << std::string("}");
    return os;
}

/********************************************/
/********************************************/
/***************** MANAGER ******************/
/********************************************/
/********************************************/

QuantifiedVariablesManager::QuantifiedVariablesManager(Options options) : options(options) {}

void QuantifiedVariablesManager::addExistVarInstance(Variable eVar) {
    ++numberOfUsedExistVars[eVar];
    if (numberOfUsedExistVars[eVar] == 1) { // if eVar is newly added
        ++numberOfExistVars;
        existVarsDependencies[eVar] = {};
    }
}

void QuantifiedVariablesManager::removeExistVarInstance(Variable eVar) {
    --numberOfUsedExistVars[eVar];
    if (numberOfUsedExistVars[eVar] == 0) { // if eVar is not in manager anymore
        --numberOfExistVars;
        removeDependency(eVar, existVarsDependencies[eVar]);
        existVarsDependencies.erase(eVar);
    }
}

void QuantifiedVariablesManager::addUnivVarInstance(Variable uVar) {
    ++numberOfUsedUnivVars[uVar];
    if (numberOfUsedUnivVars[uVar] == 1) { // if uVar is newly added
        ++numberOfUnivVars;
        univVarsDependencies[uVar] = {};
    }
    
}

void QuantifiedVariablesManager::removeUnivVarInstance(Variable uVar) {
    --numberOfUsedUnivVars[uVar];
    if (numberOfUsedUnivVars[uVar] == 0) { // if uVar is not in manager anymore
        --numberOfUnivVars;
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

unsigned QuantifiedVariablesManager::getNumberOfUnivVars() {
    return numberOfUnivVars;
}

unsigned QuantifiedVariablesManager::getNumberOfExistVars() {
    return numberOfExistVars;
}

/********************************************/
/********************************************/
/*************** MANIPULATOR ****************/
/********************************************/
/********************************************/

//QuantifiedVariablesManipulator::QuantifiedVariablesManipulator() : internalQVManager(), qvMgr(&internalQVManager) {}

QuantifiedVariablesManipulator::QuantifiedVariablesManipulator(QuantifiedVariablesManager &qvMgr) : qvMgr(&qvMgr) {}

QuantifiedVariablesManipulator::QuantifiedVariablesManipulator(const QuantifiedVariablesManipulator &qvm) {
    qvMgr = qvm.qvMgr;
    for (Variable uVar : qvm.getUnivVars()) {
        addUnivVar(uVar);
    }

    for (Variable eVar : qvm.getExistVars()) {
        addExistVar(eVar);
    }

    //supportSet = qvm.supportSet;
}


QuantifiedVariablesManipulator::~QuantifiedVariablesManipulator() {
    for (const Variable &uVar : getUnivVars()) {
        qvMgr->removeUnivVarInstance(uVar);
    }

    for (const Variable &eVar : getExistVars()) {
        qvMgr->removeExistVarInstance(eVar);
    }
}

QuantifiedVariablesManager* QuantifiedVariablesManipulator::getManager() {
    return qvMgr;
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
    if (!isVarHereQuantified(var))
        return;

    if (isVarUniv(var)) { // var is universal var
        removeUnivVar(var);
    } else if (isVarExist(var)) { // var is exist var
        removeExistVar(var);
    }
}

VariableSet const &QuantifiedVariablesManipulator::getExistVarDependencies(Variable eVar) const {
    if (!isVarExist(eVar))
        throw DQBDDexception("Trying to get the dependency set of variable which is not existential.");
    return qvMgr->getExistVarDependencies(eVar);
}

VariableSet const &QuantifiedVariablesManipulator::getUnivVarDependencies(Variable uVar) const {
    if (!isVarUniv(uVar))
        throw DQBDDexception("Trying to get the set of dependent variable for variable which is not universal.");
    return qvMgr->getUnivVarDependencies(uVar);
}

bool QuantifiedVariablesManipulator::isVarHereQuantified(Variable var) const {
    return (univVars.contains(var) || existVars.contains(var));
}

bool QuantifiedVariablesManipulator::isVarUniv(Variable var) const {
    return qvMgr->isVarUniv(var);
}

bool QuantifiedVariablesManipulator::isVarExist(Variable var) const {
    return qvMgr->isVarExist(var);
}

void QuantifiedVariablesManipulator::clear() {
    VariableSet varsCopy = getUnivVars();
    for (const Variable &uVar : varsCopy) {
        removeUnivVar(uVar);
    }

    varsCopy = getExistVars();
    for (const Variable &eVar : varsCopy) {
        removeExistVar(eVar);
    }
}

VariableSet const &QuantifiedVariablesManipulator::getSupportSet() {
    return supportSet;
}

// removes variables that are not in the support set of matrix
void QuantifiedVariablesManipulator::removeUnusedVars() {
    VariableSet usedVars = getSupportSet();

    VariableSet varsToRemove;
    for (const Variable &uVar : getUnivVars()) {
        if (!usedVars.contains(uVar)) {
            varsToRemove.insert(uVar);
        }
    }

    for (const Variable &eVar : getExistVars()) {
        if (!usedVars.contains(eVar)) {
            varsToRemove.insert(eVar);
        }
    }

    for (const Variable &varToRemove : varsToRemove) {
        removeVar(varToRemove);
    }
}

std::ostream& QuantifiedVariablesManipulator::print(std::ostream& out) const {
    return out;
}

std::ostream& operator<<(std::ostream& out, const QuantifiedVariablesManipulator& qvm)
{
    // print all univ variables
    for (const Variable &uVar : qvm.getUnivVars()) {
        out << std::string("∀") << uVar << std::string(" ");       
    }

    //print all exist variables
    for (const Variable &eVar : qvm.getExistVars()) {
        out << std::string("∃") << eVar << std::string("{");
        // print the set of dependencies of exist variable
        auto size = qvm.getExistVarDependencies(eVar).size();
        for (const Variable &uVar : qvm.getExistVarDependencies(eVar)) {
            out << uVar;
            // the last dependency does not have ',' after it
            if (size != 1) {
                out << std::string(",");
            }
            --size;
        }
        out << std::string("} ");
    }
    
    // print the inner formula
    return qvm.print(out);
}