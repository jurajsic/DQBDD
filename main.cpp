#include <fstream>
#include <iostream>
//#include "BDDProcessor.hpp"
#include "solver.hpp"

Formula getTestFormula1(Cudd &mgr) {
    Formula f(mgr);
    // TODO constructor variable(BDD) was deleted, can be fixed?
    Variable y1(0,mgr);
    Variable y2(1,mgr);
    Variable x1(2,mgr);
    Variable x2(3,mgr);
    f.addExistVar(y1);
    f.addExistVar(y2);
    f.addUnivVar(x1);
    f.addDependency(y1, x1);
    f.addUnivVar(x2);
    f.addDependency(y2, x2);
    BDD m = x1 & x2;
    m = m.Xnor(y1.getBDD().Xnor(y2));
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