#include <fstream>
#include <iostream>
//#include "BDDProcessor.hpp"
#include "Solver.hpp"

Formula getTestFormula1(Cudd &mgr) {
    Formula f(mgr);
    Variable y1 = mgr.bddVar();
    Variable y2 = mgr.bddVar();
    Variable x1 = mgr.bddVar();
    Variable x2 = mgr.bddVar();
    f.addExistVar(y1);
    f.addExistVar(y2);
    f.addUnivVar(x1);
    f.addDependency(y1, x1);
    f.addUnivVar(x2);
    f.addDependency(y2, x2);
    BDD m = x1 & x2;
    m = m.Xnor(y1.getRepr().Xnor(y2));
    f.setMatrix(m);
    return f;
}

int main(int argc, char **argw)
{
    Cudd mgr;
    //mgr.AutodynDisable();
    std::ifstream input_file;
    Solver solver(mgr);
    if (argc > 1) {
        input_file.open(argw[1]);
        if (!input_file.is_open()) {
            std::cerr << "Could not open input file." << std::endl;
            return -1;
        }
        solver.readFile(input_file);
    } else {
        solver.setFormula(getTestFormula1(mgr));
    }
    
    if (solver.solve()) {
        std::cout << "SAT" << std::endl;
    } else {
        std::cout << "UNSAT" << std::endl;
    }
}