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
}

void BDDProcessor::initialize(int nodeNum, int cacheSize) {
    if (!isInitialized) {
        bdd_init(nodeNum, cacheSize);
    } else {
        throw "BDD processor is trying to be initialized again.";
    }

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
        nextFreeVarId = 0;
    }
}

void BDDProcessor::addVars(int numOfVars) {
    if (bdd_extvarnum(numOfVars) != this->numOfVars) {
        throw "Could not add more variables.";
    }
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

bdd BDDProcessor::getBDDRepr(Variable var) {
    return bdd_ithvar(var.getId());
}

bdd BDDProcessor::getBDDReprNeg(Variable var) {
    return bdd_nithvar(var.getId());
}