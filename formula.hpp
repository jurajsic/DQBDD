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
    Formula(const Cudd &mgr);
    Formula(const Formula &f) = default;


    BDD getMatrix() const;
    void setMatrix(const BDD &matrix);

    void removeUnusedVars();

    void eliminateUnivVar(Variable uVarToEliminate);

    void eliminateExistVar(Variable existVarToEliminate);
    void eliminateExistVars(std::vector<Variable> existVarsToEliminate);
    void eliminatePossibleExistVars();
};

#endif