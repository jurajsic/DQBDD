#ifndef BDDPAIR_HPP
#define BDDPAIR_HPP

#include "bdd.h"
#include "Variable.hpp"

class BDDPair {
private:
    bddPair* pair = bdd_newpair();
public:
    BDDPair() = default;
    ~BDDPair();
    void addToPair(Variable oldVar, Variable newVar);
    bddPair* getPair();
};

#endif