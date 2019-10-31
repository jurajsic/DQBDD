#include <numeric>
#include <cassert>
#include <algorithm>
#include <list>
#include <unordered_map>
#include "BDDProcessor.hpp"

/*
BDDProcessor::BDDProcessor(int nodeNum, int cacheSize) : BDDProcessor() {
    initialize(nodeNum, cacheSize);
}

BDDProcessor::BDDProcessor(int nodeNum, int cacheSize, int numOfVars) : BDDProcessor() {
    initialize(nodeNum, cacheSize);
    setNumOfVars(numOfVars);
}
*/



BDDProcessor::~BDDProcessor() {
    if (isInitialized) {
        bdd_done();
    }
    isInitialized = false;
}

void BDDProcessor::initialize(int nodeNum, int cacheSize) {
    if (!isInitialized) {
        bdd_init(nodeNum, cacheSize);
    } else {
        throw "BDD processor is trying to be initialized again.";
    }

    bdd_setmaxincrease(500000);

    isInitialized = true;
}

void BDDProcessor::setNumOfVars(int numOfVars) {
    if (this->numOfVars > 0) {
        throw "Number of variables was already set.";
    } else {
        if (bdd_setvarnum(numOfVars) < 0) {
            throw "Could not set variables.";
        }
        this->numOfVars = numOfVars;
        nextFreeVarId = numOfVars;
        curVarOrder.resize(numOfVars,0);
        std::iota(curVarOrder.begin(), curVarOrder.end(), 0);
    }
}

void BDDProcessor::addVars(int numOfVars) {
    if (bdd_extvarnum(numOfVars) < 0) {
        throw "Could not add more variables.";
    }
    int oldNumOfVars = this->numOfVars;
    this->numOfVars += numOfVars;
    curVarOrder.resize(this->numOfVars);
    std::iota(curVarOrder.begin() + oldNumOfVars, curVarOrder.end(), oldNumOfVars);
}

Variable BDDProcessor::getFreeVariable() {
    if (numOfVars == 0) {
        setNumOfVars(10);
        return nextFreeVarId++;
    } else if (!(nextFreeVarId < numOfVars)) {
        addVars(numOfVars); // double the amount of variables
    }
    return Variable(nextFreeVarId++);
}

bool BDDProcessor::varExists(Variable var) {
    return (var.getId() < numOfVars);
}

bdd BDDProcessor::getBDDRepr(Variable var) {
    if (varExists(var)) {
        return bdd_ithvar(var.getId());
    } else {
        throw "Variable " + std::to_string(var.getId()) + " is outside of scope of defined variables.";
    }
}

bdd BDDProcessor::getBDDReprNeg(Variable var) {
    if (varExists(var)) {
        return bdd_nithvar(var.getId());
    } else {
        throw "Variable " + std::to_string(var.getId()) + " is outside of scope of defined variables.";
    }
}

void BDDProcessor::addToPairAndSetOrder(BDDPair& pair, Variable var1, Variable var2) {
    pair.addToPair(var1, var2);
    int levelOfVar1 = bdd_var2level(var1.getId());
    int levelOfVar2 = bdd_var2level(var2.getId());
    assert(curVarOrder[levelOfVar1] == var1.getId());
    assert(curVarOrder[levelOfVar2] == var2.getId());
    curVarOrder.erase(curVarOrder.begin() + levelOfVar2);
    curVarOrder.insert(curVarOrder.begin() + levelOfVar1, var2.getId());
    bdd_setvarorder(curVarOrder.data());
}

void BDDProcessor::setNewOrder(BDDPair& replacingPairs) {
    std::unordered_map<int,int> varToRepl;
    std::vector<int> oldVars;
    std::vector<int> newVars;
    for (auto pair : replacingPairs.getPairs()) {
        int oldVarId = pair.first.getId();
        int newVarId = pair.second.getId();
        varToRepl[oldVarId] = newVarId;
        oldVars.push_back(oldVarId);
        newVars.push_back(newVarId);
    }
    auto comp = [](const int& a, const int& b) -> bool {
        return (bdd_var2level(a) < bdd_var2level(b));
    };
    std::sort(oldVars.begin(), oldVars.end(), comp);
    std::sort(newVars.begin(), newVars.end(), comp);

    std::list<int> newOrder(curVarOrder.begin(),curVarOrder.end());

    auto oldVarsIt = oldVars.begin();
    auto newVarsIt = newVars.begin();
    for (auto it = newOrder.begin(); it != newOrder.end();) {
        int curVar = *it;
        if (newVarsIt != newVars.end() && *newVarsIt == curVar) {
            ++newVarsIt;
            it = newOrder.erase(it);
        } else {
            if (oldVarsIt != oldVars.end() && *oldVarsIt == curVar) {
                ++oldVarsIt;
                newOrder.insert(it, varToRepl[curVar]);
            }
            ++it;
        }
    }

    assert(oldVarsIt == oldVars.end() && newVarsIt == newVars.end());
    
    curVarOrder = std::vector<int>(newOrder.begin(),newOrder.end());
    std::cout << "setting it now for real with " << numOfVars << " variables" << std::endl;
    bdd_setvarorder(curVarOrder.data());
    // bbd_setvarordernochange(curVarOrder.data()); // does not work
}