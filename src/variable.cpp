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

std::ostream& operator<<(std::ostream& os, const Variable& obj)
{
    os << obj.getBDD();
    return os;
}