#ifndef VARIABLE_HPP
#define VARIABLE_HPP

#include "cuddObj.hh"

class Variable {
private:
    unsigned int id;
    BDD representation;

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
    //BDD operator&(const Variable& other) const;
    //BDD operator|(const Variable& other) const;
    BDD operator&(const BDD& other) const;
    BDD operator|(const BDD& other) const;
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