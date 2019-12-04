#ifndef DQBF_BDD_SIMPLESOLVER_HPP
#define DQBF_BDD_SIMPLESOLVER_HPP

#include "formula.hpp"
#include "solver.hpp"

class SimpleSolver : public Solver {
private:
    QuantifiedVariablesManager qvMgr;
    Formula formula;
    Variable getSomeUnivVar(int choice=0);
    void reorder();

    std::vector<Variable> univVarsOrderToRemove;
    void setUnivVarsOrder();
public:
    SimpleSolver() = delete;
    SimpleSolver(const Cudd &mgr);
    void readFile(std::ifstream& file);
    void setTest1Formula();
    void setTest2Formula();
    //void setFormula(Formula formula);
    bool solve();
    void runTests();
};

#endif