#include <functional>
#include "Variable.hpp"

Variable::Variable(int id, Cudd &mgr) : id(id) {
    representation = mgr.bddVar(id);
}

Variable::Variable(BDD repr) : representation(repr) {
    id = repr.NodeReadIndex();
}

Variable::operator BDD() {
    return representation;
}

unsigned int Variable::getId() const {
    return id;
}


BDD Variable::getRepr() const {
    return representation;
}


bool Variable::operator==(const Variable &anotherVariable) const {
    return (id == anotherVariable.id);
}

BDD Variable::operator&(const BDD& other) const {
    return representation & other;
}

BDD Variable::operator|(const BDD& other) const {
    return representation | other;
}