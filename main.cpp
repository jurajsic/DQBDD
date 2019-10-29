#include <fstream>
#include "BDDProcessor.hpp"
#include "Solver.hpp"

Formula getTestFormula1() {
    BDDProcessor& BDDPr = BDDProcessor::getInstance();
    BDDPr.initialize(100,100);
    BDDPr.setNumOfVars(4);
    Formula f;
    Variable y1 = Variable(0);
    Variable y2 = Variable(1);
    Variable x1 = Variable(2);
    Variable x2 = Variable(3);
    f.addExistVar(y1);
    f.addExistVar(y2);
    f.addUnivVar(x1);
    f.addDependency(y1, x1);
    f.addUnivVar(x2);
    f.addDependency(y2, x2);
    f.setMatrix(bdd_biimp(bdd_ithvar(0) & bdd_ithvar(1), bdd_biimp(bdd_ithvar(2), bdd_ithvar(3))));
    return f;
}

int main(int argc, char **argw)
{
    try {
    std::ifstream input_file;
    Solver solver;
    if (argc > 1) {
        input_file.open(argw[1]);
        if (!input_file.is_open()) {
            std::cerr << "Could not open input file." << std::endl;
            return -1;
        }
        solver.readFile(input_file);
    } else {
        solver.setFormula(getTestFormula1());
    }
    
    if (solver.solve()) {
        std::cout << "SAT" << std::endl;
    } else {
        std::cout << "UNSAT" << std::endl;
    }
    } catch (char const* e) {
        std::cout << e << std::endl;
    }
}