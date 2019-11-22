#ifndef SIMPLESOLVER_HPP
#define SIMPLESOLVER_HPP

#include <fstream>
#include "formula.hpp"
#include "solver.hpp"

class SimpleSolver : Solver {
private:
    QuantifiedVariablesManager qvMgr;
    Formula formula;
    Variable getSomeUnivVar(int choice=0);
    void printFormulaStats();
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
};

#endif