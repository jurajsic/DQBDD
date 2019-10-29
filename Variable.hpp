#ifndef VARIABLE_HPP
#define VARIABLE_HPP

#include "bdd.h"

class Variable {
private:
    int id;

public:
    Variable() = delete;
    Variable(int id);
    int getId() const;
    bool operator==(const Variable &anotherVariable) const;
};

namespace std
{
    template <>
    struct hash<Variable>
    {
        size_t operator()(const Variable& k) const
        {
            return hash<int>()(k.getId());
        }
    };
}

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