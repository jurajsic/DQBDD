#include <fstream>
#include "BDDProcessor.hpp"
#include "Solver.hpp"

Formula getTestFormula1() {
    BDDProcessor& BDDPr = BDDProcessor::getInstance();
    BDDPr.initialize(100,100);
    BDDPr.setNumOfVars(10);
    Formula f;
    f.addExistVar(Variable(0));
    f.addExistVar(Variable(1));
    f.addUnivVar(Variable(2));
    f.addUnivVar(Variable(3));
    f.setMatrix(bdd_biimp(bdd_ithvar(0) & bdd_ithvar(1), bdd_biimp(bdd_ithvar(2), bdd_ithvar(3))));
    return f;
}

int main(int argc, char **argw)
{
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

}