#ifndef BDDPAIR_HPP
#define BDDPAIR_HPP

#include <vector>
#include "bdd.h"
#include "Variable.hpp"

class BDDPair {
private:
    bddPair* pair = bdd_newpair();
    std::vector<std::pair<Variable,Variable>> pairs;
public:
    BDDPair() = default;
    ~BDDPair();
    void addToPair(Variable oldVar, Variable newVar);
    bddPair* getPair() const;
    std::vector<std::pair<Variable,Variable>> getPairs() const;
};

#endif