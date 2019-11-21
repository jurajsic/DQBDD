#ifndef FORMULA_HPP
#define FORMULA_HPP

#include <unordered_map>
#include <unordered_set>
#include "variable.hpp"
#include "quantifiedvariablesmanipulator.hpp"
#include "cuddObj.hh"

class Formula : public QuantifiedVariablesManipulator {
private:

    // propositional predicate inside DQBF formula in BDD form
    BDD matrix;

    Cudd mgr;

public:
    Formula() = delete;
    //Formula(const Cudd &mgr);
    Formula(const Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    Formula(const Formula &f) = default;

    VariableSet getSupportSet();

    BDD getMatrix() const;
    void setMatrix(const BDD &matrix);

    void removeUnusedVars();

    void eliminateUnivVar(Variable uVarToEliminate);

    void eliminateExistVar(Variable existVarToEliminate);
    void eliminateExistVars(VariableSet existVarsToEliminate);
    // eliminates all existential variables that are possible to eliminate based on Theorem 5 from DQBF localization paper
    // returns number of eliminated existential variables
    int eliminatePossibleExistVars();
};

#endif