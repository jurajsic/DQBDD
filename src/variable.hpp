#ifndef DQBF_BDD_VARIABLE_HPP
#define DQBF_BDD_VARIABLE_HPP

#include <unordered_set>
#include "cuddObj.hh"

class Variable {
private:
    // each variable has unique ID based on Cudd manager
    unsigned int id;
    BDD representation;
    Cudd mgr;

public:
    operator BDD() const;
    Variable() = delete;
    Variable(int id, Cudd &mgr);
    //Variable(BDD repr);
    unsigned int getId() const;
    BDD getBDD() const;
    // get a new Variable that is at the same BDD level as this variable
    Variable newVarAtSameLevel();
    int getLevel() const;
    bool operator==(const Variable &anotherVariable) const;
    //BDD operator&(const Variable& other) const;
    //BDD operator|(const Variable& other) const;
    BDD operator&(const BDD& other) const;
    BDD operator|(const BDD& other) const;
    BDD operator!() const;
};

std::ostream& operator<<(std::ostream& os, const Variable& obj);

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