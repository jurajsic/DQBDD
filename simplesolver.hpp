#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <fstream>
#include "cuddObj.hh"
#include "formula.hpp"
//#include "BDDProcessor.hpp"

class Solver {
protected:
    Cudd mgr;
public:
    Solver() = delete;
    Solver(const Cudd &mgr);
    virtual void readFile(std::ifstream& file) = 0;
    virtual bool solve() = 0;
    virtual void runTests() = 0;
};

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