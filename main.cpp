#include <fstream>
#include <iostream>
#include "simplesolver.hpp"
#include "treesolver.hpp"

// TODO add some sort of checking (maybe creating dependency graph -> in manager but what 
// to do with tree) of whether formula is in QBF form??

int main(int argc, char **argw)
{
    Cudd mgr;
    //mgr.AutodynDisable();
    std::ifstream input_file;
    Solver *solver;
    if (argc > 2) {
        if (std::stoi(argw[1]) == 0) {
            solver = new SimpleSolver(mgr);
        } else {
            solver = new TreeSolver(mgr);
        }
        input_file.open(argw[2]);
        if (!input_file.is_open()) {
            std::cerr << "Could not open input file." << std::endl;
            return -1;
        }
        solver->readFile(input_file);
        if (solver->solve()) {
            std::cout << "SAT" << std::endl;
        } else {
            std::cout << "UNSAT" << std::endl;
        }
    } else {
        solver = new TreeSolver(mgr);
        solver->runTests();
        //solver.setTest2Formula();
        //std::cout << solver.solve() << std::endl;
    }
    delete solver;
    return 0;
}