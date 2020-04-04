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

#ifndef DQBDD_FORMULA_HPP
#define DQBDD_FORMULA_HPP

#include <unordered_map>
#include <unordered_set>

#include <cuddObj.hh>

#include "DQBDDvariable.hpp"
#include "quantifiedvariablesmanipulator.hpp"

class Formula : public virtual QuantifiedVariablesManipulator {
private:
    // the Cudd manager used for BDD manipulation
    Cudd mgr;
    // the matrix of DQBF formula in BDD form
    BDD matrix = mgr.bddZero();

    std::ostream& print(std::ostream& out) const override;

    bool needToRecomputeSupportSet = true;

    // used for saving the order of universal variables to eliminate    
    std::vector<Variable> univVarsOrderToRemove;
    void initializeUnivVarEliminationOrder();
    Variable getUnivVarToEliminate();

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
     * @brief Eliminates all possible vars (both universal and existential)
     */
    void eliminatePossibleVars();

    void printStats();
};

#endif