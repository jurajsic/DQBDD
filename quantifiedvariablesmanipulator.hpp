#ifndef QUANTIFIEDVARIABLESMANIPULATOR_HPP
#define QUANTIFIEDVARIABLESMANIPULATOR_HPP

#include <unordered_map>
#include "variable.hpp"


class VariableSet : public std::unordered_set<Variable> {
public:
    bool contains(Variable const &var) {
        return (this->count(var) != 0);
    }
};

//typedef std::unordered_set<Variable> VariableSet;
typedef std::unordered_map<Variable, VariableSet> DependencyMap;

class QuantifiedVariablesManipulator {
private:
    VariableSet univVars;
    VariableSet existVars;
    VariableSet vars;
    
    // maps universal variable x to the set of existential variables that depend on x
    DependencyMap &univVarsDependencies;
    // maps existential variable y to the set of universal variables that depend on which y depends
    DependencyMap &existVarsDependencies;

    //VariableSet &allUnivVars;
    //VariableSet &allExistVars;
public:
    QuantifiedVariablesManipulator();

    VariableSet getUnivVars() const;
    VariableSet getExistVars() const;
    VariableSet getVars() const;

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

    bool isUnivVar(Variable var);
    bool isExistVar(Variable var);

    bool dependsOnEverything(Variable eVar);
};

#endif