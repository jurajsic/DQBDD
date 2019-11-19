#ifndef CLAUSE_HPP
#define CLAUSE_HPP

#include <unordered_set>

typedef int CNFVariable;
typedef std::unordered_set<int> CNFVariableSet;

class Clause {
private:
    CNFVariableSet vars;
public:
    Clause() = default;
    Clause(CNFVariable var);
    Clause(CNFVariableSet vars);

    void addVar(CNFVariable var);
    void addVars(CNFVariableSet vars);

    bool hasVar(CNFVariable var);
    

};

#endif