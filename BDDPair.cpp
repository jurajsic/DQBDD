#include "BDDPair.hpp"

BDDPair::~BDDPair() {
    bdd_freepair(pair);
}

void BDDPair::addToPair(Variable oldVar, Variable newVar) {
    if (bdd_setpair(pair, oldVar.getId(), newVar.getId()) != 0) {
        throw "Not able to add pair of variables to bdd_pair.";
    }
}

bddPair* BDDPair::getPair() {
    return pair;
}