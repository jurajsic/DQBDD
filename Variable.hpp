#ifndef VARIABLE_HPP
#define VARIABLE_HPP

#include "cuddObj.hh"

class Variable {
private:
    unsigned int id;
    BDD represantation;

public:
    operator BDD();
    Variable() = delete;
    Variable(int id, Cudd &mgr);
    Variable(BDD repr);
    unsigned int getId() const;
    BDD getRepr() const;
    // TODO
    //int getLevel();
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

#endif