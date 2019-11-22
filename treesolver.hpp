#ifndef TREESOLVER_HPP
#define TREESOLVER_HPP

#include "solver.hpp"
#include "quantifiertree.hpp"

class TreeSolver : public Solver {
private:
    QuantifiedVariablesManager qvMgr;

    QuantifierTree  *qt = nullptr;
public:
    TreeSolver(const Cudd &mgr);
    void readFile(std::ifstream& file);
    bool solve();
    //void runTests();
};

#endif