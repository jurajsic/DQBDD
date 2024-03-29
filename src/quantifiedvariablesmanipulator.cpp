/*
 * This file is part of DQBDD.
 *
 * Copyright 2020, 2021 Juraj Síč
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

#include <cassert>

#include <easylogging++.hpp>

#include "quantifiedvariablesmanipulator.hpp"
#include "dqbddexceptions.hpp"

namespace dqbdd {

bool VariableSet::contains(Variable const &var) const {
    return (this->count(var) != 0);
}

bool VariableSet::isSubsetOf(const VariableSet &vs) const {
    for (auto iter = this->begin(); iter != this->end(); ++iter) {
        if (!vs.contains(*iter)) {
            return false;
        }
    }
    return true;
}

VariableSet VariableSet::intersect(const VariableSet &vs) const {
    VariableSet intersection = { };
    if (this->size() < vs.size()) {
        for (auto iter = this->begin(); iter != this->end(); ++iter) {
            if (vs.contains(*iter)) {
                intersection.insert(*iter);
            }
        }
    } else {
        for (const Variable &var : vs) {
            if (this->contains(var)) {
                intersection.insert(var);
            }
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
    VLOG(5) << "addExistVarInstance(" << eVar << ")";
    ++numberOfUsedExistVars[eVar];
    if (numberOfUsedExistVars[eVar] == 1) { // if eVar is newly added
        ++numberOfExistVars;
        existVarsDependencies[eVar] = {};
    }
}

void QuantifiedVariablesManager::removeExistVarInstance(Variable eVar) {
    --numberOfUsedExistVars[eVar];
    if (numberOfUsedExistVars[eVar] == 0) { // if eVar is not in manager anymore
        VLOG(5) << "removeExistVarInstance(" << eVar << "): removing last instance";
        --numberOfExistVars;
        removeDependency(eVar, existVarsDependencies[eVar]);
        existVarsDependencies.erase(eVar);
    } else {
        VLOG(5) << "removeExistVarInstance(" << eVar << ")";
    }
}

void QuantifiedVariablesManager::addUnivVarInstance(Variable uVar) {
    VLOG(5) << "addUnivVarInstance(" << uVar << ")";
    ++numberOfUsedUnivVars[uVar];
    if (numberOfUsedUnivVars[uVar] == 1) { // if uVar is newly added
        ++numberOfUnivVars;
        univVarsDependencies[uVar] = {};
    }
    
}

void QuantifiedVariablesManager::removeUnivVarInstance(Variable uVar) {
    --numberOfUsedUnivVars[uVar];
    if (numberOfUsedUnivVars[uVar] == 0) { // if uVar is not in manager anymore
        VLOG(5) << "removeUnivVarInstance(" << uVar << "): removing last instance";
        --numberOfUnivVars;
        // making a copy, because we will be changing univVarsDependencies[uVar] in the loop
        auto existVarsToUpdate = univVarsDependencies[uVar];
        for (Variable existVarToUpdate : existVarsToUpdate) {
            removeDependency(existVarToUpdate, VariableSet{uVar});
        }
        univVarsDependencies.erase(uVar);
    } else {
        VLOG(5) << "removeUnivVarInstance(" << uVar << ")";
    }
}

void QuantifiedVariablesManager::addDependency(Variable eVar, VariableSet dependencies) {
    VLOG(5) << "addDependency(" << eVar << ", " << dependencies << ")"; 
    for (Variable uVar : dependencies) {
        existVarsDependencies[eVar].insert(uVar);
        univVarsDependencies[uVar].insert(eVar);
    }
}

void QuantifiedVariablesManager::removeDependency(Variable eVar, VariableSet dependencies) {
    VLOG(5) << "removeDependency(" << eVar << ", " << dependencies << ")"; 
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

void QuantifiedVariablesManager::reorderVars(Cudd &mgr) {
    VariableSet univVars;
    std::vector<Variable> sortedUnivVars;
    for (auto &uVarWithDependencies : univVarsDependencies) {
        Variable univVar = uVarWithDependencies.first;
        univVars.insert(univVar);
        sortedUnivVars.push_back(univVar);
    }
    VariableSet existVars;
    for (auto &eVarWithDependencies : existVarsDependencies) {
        existVars.insert(eVarWithDependencies.first);
    }

    VLOG(5) << "reorderVars: Univ vars " << univVars;
    VLOG(5) << "reorderVars: Exist vars " << existVars;

    std::sort(sortedUnivVars.begin(), sortedUnivVars.end(),
                    [&](const Variable &a, const Variable &b) {
                        auto aDependenciesSize = getUnivVarDependencies(a).size();
                        auto bDependenciesSize = getUnivVarDependencies(b).size();
                        return (aDependenciesSize < bDependenciesSize
                                || (aDependenciesSize == bDependenciesSize && a.getId() < b.getId()));
                    }
                );
    
    VariableSet removedUnivVars;
    int i = mgr.ReadSize(); // number of BDD variables in mgr (it is possible there are more than those used in formulas)
    // the new ordering of variables
    std::vector<int> orderOfElimination(i);
    // on the lowest level are existential variables that depend on everything (they will be eliminated first)
    for (auto it = existVars.begin(); it != existVars.end(); ) {
        if (existVarsDependencies[*it].size() == univVars.size()) {
            --i; orderOfElimination[i] = it->getId();
            it = existVars.erase(it);
        } else {
            ++it;
        }
    }
    for (const Variable &univVar : sortedUnivVars) {
        removedUnivVars.insert(univVar);
        --i; orderOfElimination[i] = univVar.getId();
        univVars.erase(univVar);
        for (auto it = existVars.begin(); it != existVars.end(); ) {
            if (existVarsDependencies[*it].minus(removedUnivVars).size() == univVars.size()) {
                --i; orderOfElimination[i] = it->getId();
                it = existVars.erase(it);
            } else {
                ++it;
            }
        }
    }

    // the variables in CUDD are given from 0 to mgr.ReadSize(), it is possible that some of them are not used in formulas, but we still need to give them some order
    for (int bddVar = 0; bddVar < mgr.ReadSize(); ++bddVar) {
        if (!existVarsDependencies.count(Variable(bddVar, mgr)) && !univVarsDependencies.count(Variable(bddVar, mgr))) {
            --i; orderOfElimination[i] = bddVar;
        }
    }

    if (VLOG_IS_ON(5)) {
        std::stringstream orderLog;
        orderLog << "reorderVars: new order of BDD variables from highest to lowest: ";
        for (int j = i; j < orderOfElimination.size(); ++j) {
            orderLog << orderOfElimination[j] << ' ';
        }
        VLOG(5) << orderLog.str();
    }


    assert(i == 0);

    mgr.ShuffleHeap(orderOfElimination.data());
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
        throw dqbddException("Trying to get the dependency set of variable which is not existential.");
    return qvMgr->getExistVarDependencies(eVar);
}

VariableSet const &QuantifiedVariablesManipulator::getUnivVarDependencies(Variable uVar) const {
    if (!isVarUniv(uVar))
        throw dqbddException("Trying to get the set of dependent variable for variable which is not universal.");
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
    VLOG(5) << "removeUnusedVars for " << this;
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

void QuantifiedVariablesManipulator::printPrefix(std::ostream &out) const {
    // print all univ variables
    for (const Variable &uVar : getUnivVars()) {
        out << std::string("∀") << uVar << std::string(" ");       
    }

    //print all exist variables
    for (const Variable &eVar : getExistVars()) {
        out << std::string("∃") << eVar << std::string("{");
        // print the set of dependencies of exist variable
        auto size = getExistVarDependencies(eVar).size();
        for (const Variable &uVar : getExistVarDependencies(eVar)) {
            out << uVar;
            // the last dependency does not have ',' after it
            if (size != 1) {
                out << std::string(",");
            }
            --size;
        }
        out << std::string("} ");
    }
}

std::ostream& operator<<(std::ostream& out, const QuantifiedVariablesManipulator& qvm)
{
    // print the prefix
    qvm.printPrefix(out);
    // print the inner formula
    return qvm.print(out);
}

} // namespace dqbdd