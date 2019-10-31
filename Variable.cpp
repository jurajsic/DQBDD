#include <functional>
#include "Variable.hpp"

Variable::Variable(int id, Cudd &mgr) : id(id) {
    mgr.bddVar(id);
}

Variable::Variable(BDD repr) : represantation(repr) {
    id = repr.NodeReadIndex();
}

Variable::operator BDD() {
    return represantation;
}

unsigned int Variable::getId() const {
    return id;
}


BDD Variable::getRepr() const {
    return represantation;
}


bool Variable::operator==(const Variable &anotherVariable) const {
    return (id == anotherVariable.id);
}
