#ifndef DQBDD_FORMULA_HPP
#define DQBDD_FORMULA_HPP

#include <unordered_map>
#include <unordered_set>

#include <cuddObj.hh>

#include "DQBDDvariable.hpp"
#include "quantifiedvariablesmanipulator.hpp"

class Formula : public virtual QuantifiedVariablesManipulator {
private:

    Cudd mgr;
    // propositional predicate inside DQBF formula in BDD form
    BDD matrix = mgr.bddZero();

    std::ostream& print(std::ostream& out) const override;

    bool needToRecomputeSupportSet = true;

public:
    Formula() = delete;
    //Formula(const Cudd &mgr);
    Formula(const Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    Formula(const Cudd &mgr, QuantifiedVariablesManipulator &qvManipulator);
    Formula(const Formula &f) = delete;

    virtual ~Formula() = default;

    VariableSet const &getSupportSet() override;

    BDD getMatrix() const;
    void setMatrix(const BDD &matrix);

    void eliminateUnivVar(Variable uVarToEliminate);

    void eliminateExistVar(Variable existVarToEliminate);
    void eliminateExistVars(VariableSet existVarsToEliminate);
    VariableSet getPossibleExistVarsToEliminate();
    // eliminates all existential variables that are possible to eliminate based on Theorem 5 from DQBF localization paper
    // returns number of eliminated existential variables
    //int eliminatePossibleExistVars();

    // TODO implemenet -> should eliminate all universal variables and all possible exist (can add new exist)
    void eliminatePossibleVars();

    void printStats();
};

#endif