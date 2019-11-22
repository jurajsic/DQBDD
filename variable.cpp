#include <functional>
#include <iostream>
#include "variable.hpp"

Variable::Variable(int id, Cudd &mgr) : id(id), mgr(mgr) {
    representation = mgr.bddVar(id);
}

/*
Variable::Variable(BDD repr) : representation(repr) {
    id = repr.NodeReadIndex();
}
*/

int Variable::getLevel() const {
    return mgr.ReadPerm(id);
}

Variable Variable::newVarAtSameLevel() {
    return Variable(mgr.bddNewVarAtLevel(getLevel()).NodeReadIndex(), mgr);
}

Variable::operator BDD() const {
    return getBDD();
}

unsigned int Variable::getId() const {
    return id;
}


BDD Variable::getBDD() const {
    // TODO decide if we save the BDD - if not we have to tell mgr that id is taken!
    //return mgr.bddVar(id);
    return representation;
}


bool Variable::operator==(const Variable &anotherVariable) const {
    return (id == anotherVariable.id);
}

BDD Variable::operator&(const BDD& other) const {
    return getBDD() & other;
}

BDD Variable::operator|(const BDD& other) const {
    return getBDD() | other;
}

BDD Variable::operator!() const {
    return !getBDD();
}