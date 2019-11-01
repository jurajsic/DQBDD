#ifndef FORMULA_HPP
#define FORMULA_HPP

#include <unordered_map>
#include <unordered_set>
#include "Variable.hpp"
#include "cuddObj.hh"

typedef std::unordered_set<Variable> VariableSet;
typedef std::unordered_map<Variable, VariableSet> DependencyMap;

class Formula {
private:
    VariableSet univVars;
    VariableSet existVars;
    
    // maps universal variable x to the set of existential variables that depend on x
    DependencyMap univVarsDependencies;
    // maps existential variable y to the set of universal variables that depend on which y depends
    DependencyMap existVarsDependencies;

    // propositional predicate inside DQBF formula in BDD form
    BDD matrix;

    Cudd mgr;

public:
    Formula() = delete;
    Formula(const Cudd &mgr);
    Formula(const Formula &f) = default;
    VariableSet getUnivVars() const;
    VariableSet getExistVars() const;
    BDD getMatrix() const;
    void setMatrix(const BDD &matrix);
    void addUnivVar(Variable uVar);
    void addExistVar(Variable eVar);
    void addExistVar(Variable eVar, VariableSet dependencies);
    void addDependency(Variable eVar, Variable uVar);
    void addDependency(Variable eVar, VariableSet dependencies);
    void removeDependency(Variable eVar, Variable uVar);
    void removeDependency(Variable eVar, VariableSet dependencies);
    void removeUnivVar(Variable uVar);
    void removeExistVar(Variable eVar);
    void removeVar(Variable var);
    VariableSet getExistVarDependencies(Variable eVar);
    VariableSet getUnivVarDependencies(Variable uVar);
    bool dependsOnEverything(Variable eVar);
    void removeUnusedVars();
};

#endif