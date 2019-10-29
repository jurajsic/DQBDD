#ifndef VARIABLE_HPP
#define VARIABLE_HPP

#include "bdd.h"

class Variable {
private:
    int id;

public:
    Variable(int id);
    int getId() const;
    bool operator==(const Variable &anotherVariable) const;
};
/*
class ExistVariable : Variable {
public:
    ExistVariable(int id);
};

class UnivVariable : Variable {
public:
    UnivVariable(int id);
};
*/
#endif