#include "BDDPair.hpp"

BDDPair::~BDDPair() {
    bdd_freepair(pair);
}

void BDDPair::addToPair(Variable oldVar, Variable newVar) {
    if (bdd_setpair(pair, oldVar.getId(), newVar.getId()) != 0) {
        throw "Not able to add pair of variables to bdd_pair.";
    }
    pairs.push_back(std::pair<Variable,Variable>(oldVar, newVar));
}

bddPair* BDDPair::getPair() const {
    return pair;
}

std::vector<std::pair<Variable,Variable>> BDDPair::getPairs() const {
    return pairs;
}