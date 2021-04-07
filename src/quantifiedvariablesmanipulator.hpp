/*
 * This file is part of DQBDD.
 *
 * Copyright 2020 Juraj Síč
 *
 * DQBDD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * DQBDD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with DQBDD. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DQBDD_QUANTIFIEDVARIABLESMANIPULATOR_HPP
#define DQBDD_QUANTIFIEDVARIABLESMANIPULATOR_HPP

#include <unordered_map>
#include <ostream>

#include "dqbddvariable.hpp"
#include "dqbddoptions.hpp"

namespace dqbdd {

class VariableSet : public std::unordered_set<Variable> {
public:
    using std::unordered_set<Variable>::unordered_set;

    bool contains(Variable const &var) const;
    /**
     * @brief returns the intersection of this set with vs
     */
    VariableSet intersect(const VariableSet &vs) const;
    /**
     * @brief returns the union of this set with vs
     */
    VariableSet unite(const VariableSet &vs) const;
    /**
     * @brief returns this set minus vs
     */
    VariableSet minus(const VariableSet &vs) const;
    /**
     * @brief Returns the set of variables occuring in BDD b
     */
    static VariableSet getSupportSetOfBDD(const BDD &b, Cudd &mgr);
};

std::ostream& operator<<(std::ostream& os, const VariableSet& variableSet);


using DependencyMap = std::unordered_map<Variable, VariableSet>;

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

    unsigned numberOfUnivVars = 0;
    unsigned numberOfExistVars = 0;

public:
    Options options;

    QuantifiedVariablesManager() = default;
    QuantifiedVariablesManager(Options options);

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

/**
 * @brief Quantified variables manipulator
 * 
 * Manipulator of quantified variables to which universal and existential variables can be
 * saved. Each manipulator must be constructed with a pointer to QuantifiedVariablesManager
 * which takes care of dependencies between variables and this manager cannot be deleted
 * while we work with this manipulator.
 * It should be used as a base class for representing subformulas of DQBF, where this represents
 * the prefix of subformula, while the matrix (which can be created from subformulas) needs to be
 * represented in subclass. The sublass should also redefine the method getSupportSet() which should 
 * update the set supportSet so it contains the variables that occur in the matrix and return it.
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
    // QVManager that takes care of dependencies for variables occuring in this manipulator
    QuantifiedVariablesManager *qvMgr;

    // the set of variables occuring in the matrix (NOT those that occur in this manipulator)
    VariableSet supportSet;

public:
    QuantifiedVariablesManipulator() = delete;
    QuantifiedVariablesManipulator(QuantifiedVariablesManager &qvMgr);
    ~QuantifiedVariablesManipulator();
    QuantifiedVariablesManipulator(const QuantifiedVariablesManipulator &qvm);

    QuantifiedVariablesManipulator& operator=(QuantifiedVariablesManipulator &qvm) = delete;

    VariableSet const &getUnivVars() const;
    VariableSet const &getExistVars() const;

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
};

std::ostream& operator<<(std::ostream& os, const QuantifiedVariablesManipulator& qvm);

} // namespace dqbdd

#endif