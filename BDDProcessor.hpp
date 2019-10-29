#ifndef BDDPROCESSOR_HPP
#define BDDPROCESSOR_HPP

#include <vector>
#include "bdd.h"
#include "Variable.hpp"
#include "BDDPair.hpp"

class BDDProcessor {
private:
    int nextFreeVarId = -1;
    int numOfVars = 0;
    bool isInitialized = false;

    std::vector<int> curVarOrder;

    void addVars(int numOfVars);

    

    BDDProcessor() = default;
    //BDDProcessor(int nodeNum, int cacheSize);
    //BDDProcessor(int nodeNum, int cacheSize, int numOfVars);

public:
    // making it a singleton
    static BDDProcessor& getInstance() {
        static BDDProcessor instance;
        return instance;
    }
    BDDProcessor(BDDProcessor const&) = delete;
    void operator=(BDDProcessor const&)  = delete;

    ~BDDProcessor();
    void initialize(int nodeNum, int cacheSize);
    void setNumOfVars(int numOfVars);
    Variable getFreeVariable();
    bool varExists(Variable var);
    bdd getBDDRepr(Variable var);
    bdd getBDDReprNeg(Variable var);
    void addToPairAndSetOrder(BDDPair& pair, Variable var1, Variable var2);
    void setNewOrder(BDDPair& pair);
};

#endif