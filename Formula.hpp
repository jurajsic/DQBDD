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

public:
    Formula() = default;
    Formula(const Formula &f) = default;
    VariableSet getUnivVars() const;
    VariableSet getExistVars() const;
    BDD getMatrix() const;
    void setMatrix(BDD matrix);
    void addUnivVar(Variable uVar);
    void addExistVar(Variable eVar);
    void addExistVar(Variable eVar, VariableSet dependencies);
    void addDependency(Variable eVar, Variable uVar);
    void addDependency(Variable eVar, VariableSet dependencies);
    void removeDependency(Variable eVar, Variable uVar);
    void removeDependency(Variable eVar, VariableSet dependencies);
    void removeUnivVar(Variable uVar);
    void removeExistVar(Variable eVar);
    VariableSet getExistVarDependencies(Variable eVar);
    VariableSet getUnivVarDependencies(Variable uVar);
    bool dependsOnEverything(Variable eVar);
    // TODO isMatrixOne/Zero nahradit za toto (treba to vobec, viem asi to zistit aj cez get)
    /*
    bool isTrue();
    bool isFalse();
    bool isSimple();
    */
    // TODO implement this:
    void cleanUnusedVars();
};

#endif