#ifndef FORMULA_HPP
#define FORMULA_HPP

#include <unordered_map>
#include <unordered_set>
#include "Variable.hpp"
#include "bdd.h"

class Formula {
private:
    std::unordered_set<Variable> univVars;
    std::unordered_set<Variable> existVars;
    
    // maps universal variable x to the set of existential variables that depend on x
    std::unordered_map<Variable, std::unordered_set<Variable>> univVarsDependencies;
    // maps existential variable y to the set of universal variables that depend on which y depends
    std::unordered_map<Variable, std::unordered_set<Variable>> existVarsDependencies;

    // propositional predicate inside DQBF formula in BDD form
    bdd matrix;

public:
    Formula() = default;
    Formula(const Formula &f) = default;
    std::unordered_set<Variable> getUnivVars() const;
    std::unordered_set<Variable> getExistVars() const;
    bdd getMatrix() const;
    void setMatrix(bdd matrix);
    void addUnivVar(Variable uVar);
    void addExistVar(Variable eVar);
    void addDependency(Variable eVar, Variable uVar);
    void addDependency(Variable eVar, std::unordered_set<Variable> dependencies);
    void removeDependency(Variable eVar, Variable uVar);
    void removeDependency(Variable eVar, std::unordered_set<Variable> dependencies);
    void removeUnivVar(Variable uVar);
    void removeExistVar(Variable eVar);
    std::unordered_set<Variable> getExistVarDependencies(Variable eVar);
    std::unordered_set<Variable> getUnivVarDependencies(Variable uVar);
    bool dependsOnEverything(Variable eVar);
};

#endif