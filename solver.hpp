#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <fstream>
#include "cuddObj.hh"

class Solver {
protected:
    Cudd mgr;
public:
    Solver() = delete;
    Solver(const Cudd &mgr);
    virtual void readFile(std::ifstream& file) = 0;
    virtual bool solve() = 0;
    //virtual void runTests() = 0;
};

#endif