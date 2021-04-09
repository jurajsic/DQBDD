/*
 * This file is part of DQBDD.
 *
 * Copyright 2020, 2021 Juraj Síč
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

#ifndef DQBDD_FORMULA_HPP
#define DQBDD_FORMULA_HPP

#include <unordered_map>
#include <unordered_set>
#include <ostream>

#include <cuddObj.hh>

#include "dqbddvariable.hpp"
#include "quantifiedvariablesmanipulator.hpp"

namespace dqbdd {

class Formula : public virtual QuantifiedVariablesManipulator {
private:
    // the Cudd manager used for BDD manipulation
    Cudd mgr;
    // the matrix of DQBF formula in BDD form
    BDD matrix = mgr.bddZero();

    std::ostream& print(std::ostream& out) const override;

    bool needToRecomputeSupportSet = true;

    // used for saving the order of universal variables to eliminate (the last variable in vector is eliminated as first)    
    std::vector<Variable> univVarsOrderToRemove;
    // used for NumOfLeftoverVarsInConjuncts, we save the conjuncts for minimal variable here
    BDD minf1,minf2;
    void initializeUnivVarEliminationOrder();
    Variable getUnivVarToEliminate();

protected:
    bool needsToRecomputeSupportSet() const;

public:
    Formula() = delete;
    Formula(const Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    Formula(const Cudd &mgr, QuantifiedVariablesManipulator &qvManipulator);
    Formula(const Formula &f) = delete;

    virtual ~Formula() = default;

    VariableSet const &getSupportSet() override;

    BDD getMatrix() const;
    void setMatrix(const BDD &matrix);

    /**
     * @brief Eliminates the universal variable uVarToEliminate
     * Based on universal expansion:
     * \forall x \psi = \psi[0/x] /\ \psi[y1',...,yn'/y1,...,yn][1/x]
     * where y1,...,yn are all existential variables dependent on x
     * and y1',...,yn' are new variables.
     */
    void eliminateUnivVar(Variable uVarToEliminate);
    /**
     * @brief Eliminates the universal variable uVarToEliminate
     * Based on universal expansion:
     * \forall x \psi = \psi[0/x] /\ \psi[y1',...,yn'/y1,...,yn][1/x]
     * where y1,...,yn are all existential variables dependent on x
     * and y1',...,yn' are new variables.
     * @param useAlreadyComputedf1f2 If true, use minf1 and minf2 as 
     * \psi[0/x] and \psi[y1',...,yn'/y1,...,yn][1/x]
     */
    void eliminateUnivVar(Variable uVarToEliminate, bool useAlreadyComputedf1f2);

    /**
     * @brief Eliminates the existential variable existVarToEliminate
     * Based on existential elimination:
     * \exists y \psi = \psi[0/y] \/ \psi[1/y]
     * It is assumed that universal variables in matrix are only those
     * that existVarToEliminate depends on and all existential variables
     * in matrix depend only on those variables that existVarToEliminate 
     * depends on.
     */
    void eliminateExistVar(Variable existVarToEliminate);
    /**
     * @brief Eliminates the existential variables existVarsToEliminate
     * Based on existential elimination:
     * \exists y \psi = \psi[0/y] \/ \psi[1/y]
     * It is assumed that universal variables in matrix are only those
     * that all existVarsToEliminate depend on and all existential variables
     * in matrix depend only on those variables that existVarsToEliminate 
     * depend on.
     */
    void eliminateExistVars(VariableSet existVarsToEliminate);
    /**
     * @brief Get all existential variables that can be eliminated
     * 
     * @return existential variables y1,...,yn that depend on all universal
     * variables in matrix and all existential variables in matrix
     * depend only on those universal variables that y1,..,yn depend on 
     */
    VariableSet getPossibleExistVarsToEliminate();

    /**
     * @brief Eliminate all possible existential variables
     * 
     * @return true if something was eliminated
     */
    bool eliminatePossibleExistVars();

    /**
     * @brief Eliminate all universal variables that do not create new copies after elimination
     * 
     * @return true if something was eliminated
     */
    bool eliminateSimpleUnivVars();

    /**
     * @brief Eliminates all possible vars (both universal and existential)
     */
    void eliminatePossibleVars();

    void printStats();
};

} // namespace dqbdd

#endif