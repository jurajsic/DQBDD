#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <fstream>
#include "cuddObj.hh"
#include "Formula.hpp"
//#include "BDDProcessor.hpp"

class Solver {
private:
    Cudd mgr;
    Formula formula;
    Variable getSomeUnivVar();
public:
    Solver() = delete;
    Solver(const Cudd &mgr);
    void readFile(std::ifstream& file);
    void setFormula(Formula formula);
    bool solve();
};

#endif