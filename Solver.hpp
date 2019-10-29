#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <fstream>
#include "bdd.h"
#include "Formula.hpp"
#include "BDDProcessor.hpp"

class Solver {
private:
    Formula formula;
    BDDProcessor& bddProcessor = BDDProcessor::getInstance();
    Variable getSomeUnivVar();
    bdd getVarBDDFromStr(std::string strVar);
public:
    Solver() = default;
    Solver(std::ifstream& file);
    Solver(Formula formula);
    void readFile(std::ifstream& file);
    void setFormula(Formula formula);
    bool solve();
};

#endif