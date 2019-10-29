#ifndef BDDPROCESSOR_HPP
#define BDDPROCESSOR_HPP

#include "bdd.h"
#include "Variable.hpp"

class BDDProcessor {
private:
    int nextFreeVarId = -1;
    int numOfVars = 0;
    bool isInitialized = false;
    void addVars(int numOfVars);

    

    BDDProcessor() = default;
    //BDDProcessor(int nodeNum, int cacheSize);
    //BDDProcessor(int nodeNum, int cacheSize, int numOfVars);

public:
    // making it a singleton
    static BDDProcessor& BDDProcessor::getInstance() {
        static BDDProcessor instance;
        return instance;
    }
    BDDProcessor(BDDProcessor const&) = delete;
    void operator=(BDDProcessor const&)  = delete;

    ~BDDProcessor();
    void initialize(int nodeNum, int cacheSize);
    void setNumOfVars(int numOfVars);
    Variable getFreeVariable();
    bdd getBDDRepr(Variable var);
    bdd getBDDReprNeg(Variable var);
};

#endif