#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <fstream>
#include "cuddObj.hh"
#include "formula.hpp"
//#include "BDDProcessor.hpp"

class Solver {
private:
    Cudd mgr;
    Formula formula;
    Variable getSomeUnivVar(int choice=0);
    void printFormulaStats();
    void reorder();

    std::vector<Variable> univVarsOrderToRemove;
    void setUnivVarsOrder();
public:
    Solver() = delete;
    Solver(const Cudd &mgr);
    void readFile(std::ifstream& file);
    void setFormula(Formula formula);
    bool solve();
};

#endif