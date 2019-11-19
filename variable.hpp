#ifndef VARIABLE_HPP
#define VARIABLE_HPP

#include <unordered_set>
#include "cuddObj.hh"

class Variable {
private:
    unsigned int id;
    BDD representation;
    Cudd mgr;

public:
    operator BDD();
    Variable() = delete;
    Variable(int id, Cudd &mgr);
    //Variable(BDD repr);
    unsigned int getId() const;
    BDD getBDD() const;
    // get Variable that
    Variable newVarAtSameLevel();
    int getLevel();
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