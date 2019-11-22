#ifndef QUANTIFIEDVARIABLESMANIPULATOR_HPP
#define QUANTIFIEDVARIABLESMANIPULATOR_HPP

#include <unordered_map>
#include "variable.hpp"

class VariableSet : public std::unordered_set<Variable> {
public:
    // TODO check if this works 
    using std::unordered_set<Variable>::unordered_set;

    bool contains(Variable const &var) {
        return (this->count(var) != 0);
    }
};

//typedef std::unordered_set<Variable> VariableSet;
typedef std::unordered_map<Variable, VariableSet> DependencyMap;

class QuantifiedVariablesManager {
private:
    // maps universal variable x to the set of existential variables that depend on x
    DependencyMap univVarsDependencies;
    // maps existential variable y to the set of universal variables that depend on which y depends
    DependencyMap existVarsDependencies;

    std::unordered_map<Variable,int> numberOfUsedExistVars;
    std::unordered_map<Variable,int> numberOfUsedUnivVars;
public:
    void addExistVarInstance(Variable eVar);
    void removeExistVarInstance(Variable eVar);

    void addUnivVarInstance(Variable uVar);
    void removeUnivVarInstance(Variable uVar);

    
    void addDependency(Variable eVar, VariableSet dependencies);
    void removeDependency(Variable eVar, VariableSet dependencies);

    VariableSet const &getExistVarDependencies(Variable eVar);
    VariableSet const &getUnivVarDependencies(Variable uVar);
    
    bool isVarUniv(Variable var);
    bool isVarExist(Variable var);
};

class QuantifiedVariablesManipulator {
private:
    VariableSet univVars;
    VariableSet existVars;
    
    // internal QVManager that is used if external is not supplied in constructor
    //QuantifiedVariablesManager internalQVManager;

protected:
    // this either refers to external manager or internal if external is not supplied -> always to external??
    QuantifiedVariablesManager *qvMgr;

public:
    QuantifiedVariablesManipulator() = delete;
    QuantifiedVariablesManipulator(QuantifiedVariablesManager &qvMgr);
    ~QuantifiedVariablesManipulator();
    //QuantifiedVariablesManipulator(const QuantifiedVariablesManipulator &qvm);

    //TODO!!!!!!!!!!!!!
    //QuantifiedVariablesManipulator& operator=(QuantifiedVariablesManipulator &qvm);

    VariableSet const &getUnivVars() const;
    VariableSet const &getExistVars() const;
    //VariableSet getVars() const;

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

    VariableSet const &getExistVarDependencies(Variable eVar);
    VariableSet const &getUnivVarDependencies(Variable uVar);

    bool isVarHere(Variable var);
    bool isVarUniv(Variable var);
    bool isVarExist(Variable var);

    void clear();

    //bool dependsOnEverything(Variable eVar);
};

#endif