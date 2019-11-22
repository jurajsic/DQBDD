#ifndef TREESOLVER_HPP
#define TREESOLVER_HPP

#include "solver.hpp"
#include "quantifiertree.hpp"

class TreeSolver : public Solver {
private:
    QuantifiedVariablesManager qvMgr;

    QuantifierTreeNode *root = nullptr;
public:
    TreeSolver(const Cudd &mgr);
    ~TreeSolver();
    void readFile(std::ifstream& file);
    bool solve();
    //void runTests();
    void setTest1Formula();
};

#endif