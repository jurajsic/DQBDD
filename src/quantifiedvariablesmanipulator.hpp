#ifndef DQBDD_QUANTIFIEDVARIABLESMANIPULATOR_HPP
#define DQBDD_QUANTIFIEDVARIABLESMANIPULATOR_HPP

#include <unordered_map>

#include "DQBDDvariable.hpp"

class VariableSet : public std::unordered_set<Variable> {
public:
    using std::unordered_set<Variable>::unordered_set;

    bool contains(Variable const &var) const;
};


typedef std::unordered_map<Variable, VariableSet> DependencyMap;

/**
 * @brief Manager of quantified variables.
 * 
 * This manager saves the dependencies of existential and universal variables and has
 * methods for adding and removing instances of variables and for adding, removing and
 * retrieving the dependencies. 
 * It should be used only with QuantifiedVariablesManipulator, where multiple manipulators can
 * be under the same QuantifiedVariablesManager.
 */
class QuantifiedVariablesManager {
private:
    // maps universal variable x to the set of existential variables that depend on x
    DependencyMap univVarsDependencies;
    // maps existential variable y to the set of universal variables that depend on which y depends
    DependencyMap existVarsDependencies;

    std::unordered_map<Variable,int> numberOfUsedExistVars;
    std::unordered_map<Variable,int> numberOfUsedUnivVars;

    unsigned numberOfUnivVars;
    unsigned numberOfExistVars;
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

    unsigned getNumberOfUnivVars();
    unsigned getNumberOfExistVars();
    //VariableSet getAllUnivVars();
    //VariableSet getAllExistVars();
};

// TODO: solve supportSet, ci treba mat tu ulozene a ako to funguje diamond pri virtual metode getSupportSet()

/**
 * @brief Quantified variables manipulator
 * 
 * Manipulator of quantified variables to which universal and existential variables can be
 * saved. Each manipulator must be constructed with a pointer to QuantifiedVariablesManager
 * which takes care of dependencies between variables and this manager cannot be deleted
 * while we work with this manipulator.
 * It should be used as a base class for representing subformulas of DQBF, where this represents
 * the prefix of subformula, while the matrix (which can be created from subformulas) needs to be
 * represented in subclass. The sublass needs to also define method getSupportSet() which should 
 * return the variables that occur in the matrix.
 */
class QuantifiedVariablesManipulator {
private:
    // the sets of universal and existential variables which are saved in this manipulator
    VariableSet univVars;
    VariableSet existVars;
    
    // internal QVManager that is used if external is not supplied in constructor
    //QuantifiedVariablesManager internalQVManager;

    // used for printing the matrix
    virtual std::ostream& print(std::ostream& out) const;
    friend std::ostream& operator<<(std::ostream& out, const QuantifiedVariablesManipulator& qvm);

protected:
    // this either refers to external manager or internal if external is not supplied -> always to external??
    QuantifiedVariablesManager *qvMgr;
    // TODO: check if this is needed here, maybe should be deleted, if it is used in formula or tree, should be inside, it doesnt make sense to keep it here
    // the variables in support set do not need to be the same which are in univVars + existVars,
    // depends on the 
    VariableSet supportSet;

public:
    QuantifiedVariablesManipulator() = delete;
    QuantifiedVariablesManipulator(QuantifiedVariablesManager &qvMgr);
    ~QuantifiedVariablesManipulator();
    QuantifiedVariablesManipulator(const QuantifiedVariablesManipulator &qvm);

    //TODO!!!!!!!!!!!!!
    //QuantifiedVariablesManipulator& operator=(QuantifiedVariablesManipulator &qvm);

    VariableSet const &getUnivVars() const;
    VariableSet const &getExistVars() const;
    //VariableSet getVars() const;

    QuantifiedVariablesManager* getManager();

    void addUnivVar(Variable uVar);
    void addExistVar(Variable eVar);
    void addExistVar(Variable eVar, VariableSet dependencies);

    void addDependency(Variable eVar, Variable uVar);
    void addDependency(Variable eVar, VariableSet dependencies);

    void removeDependency(Variable eVar, Variable uVar);
    void removeDependency(Variable eVar, VariableSet dependencies);

    // does nothing if uVar is not saved here
    void removeUnivVar(Variable uVar);
    // does nothing if eVar is not saved here
    void removeExistVar(Variable eVar);
    // removes (univ or exist) var from this manipulator if it is saved here, otherwise does nothing
    void removeVar(Variable var);

    VariableSet const &getExistVarDependencies(Variable eVar) const;
    VariableSet const &getUnivVarDependencies(Variable uVar) const;

    // checks if var is saved in this manipulator 
    bool isVarHereQuantified(Variable var) const;
    // check if var is universal in this manipulator's vars manager
    bool isVarUniv(Variable var) const;
    // check if var is existential in this manipulator's vars manager
    bool isVarExist(Variable var) const;

    // removes all saved variables from this manipulator
    void clear();

    // removes all variables not in support set
    void removeUnusedVars();
    virtual VariableSet const &getSupportSet();

    //bool dependsOnEverything(Variable eVar);
};

std::ostream& operator<<(std::ostream& os, const QuantifiedVariablesManipulator& qvm);

#endif