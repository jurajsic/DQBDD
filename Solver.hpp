#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <fstream>
#include "cuddObj.hh"
#include "Formula.hpp"
//#include "BDDProcessor.hpp"

class Solver {
private:
    Formula formula;
    //BDDProcessor& bddProcessor = BDDProcessor::getInstance();
    Cudd mgr;
    Variable getSomeUnivVar();
    BDD getVarBDDFromStr(std::string strVar);
public:
    Solver() = delete;
    Solver(Cudd mgr);
    Solver(std::ifstream& file);
    Solver(Formula formula);
    void readFile(std::ifstream& file);
    void setFormula(Formula formula);
    bool solve();
};

#endif