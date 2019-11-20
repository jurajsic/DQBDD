#include <fstream>
#include <iostream>
//#include "BDDProcessor.hpp"
#include "solver.hpp"

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
        solver.setTest1Formula();
    }
    
    if (solver.solve()) {
        std::cout << "SAT" << std::endl;
    } else {
        std::cout << "UNSAT" << std::endl;
    }
}