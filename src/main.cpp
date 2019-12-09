#include <fstream>
#include <iostream>
#include "simplesolver.hpp"
#include "treesolver.hpp"
#include "HQSpreinterface.hpp"

// TODO add some sort of checking (maybe creating dependency graph -> in manager but what 
// to do with tree) of whether formula is in QBF form??

int main(int argc, char **argw)
{
    Cudd mgr;
    //mgr.AutodynDisable();
    std::ifstream input_file;
    Solver *solver;
    if (argc <= 1) {
        throw "wrong args";
        return -1;
    } else {
        if (std::stoi(argw[1]) == 0) {
            solver = new SimpleSolver(mgr);
        } else {
            solver = new TreeSolver(mgr);
        }
    }

    if (argc > 2) {

/**************************/
        delete solver;
        QuantifiedVariablesManager qvMgr;
        Parser *parser = new HQSPreInterface(mgr, qvMgr);
        parser->parse(argw[2]);
        std::cout << "Parsing finished" << std::endl;
        Formula *f = parser->getFormula();
        std::cout << "FOrmula got" << std::endl;
        delete parser;
        f->eliminatePossibleVars();
        if (f->getMatrix().IsOne()) {
            std::cout << "SAT" << std::endl;
        } else {
            std::cout << "UNSAT" << std::endl;
        }
        delete f;
        return 0;
/**********************/

        input_file.open(argw[2]);
        if (!input_file.is_open()) {
            std::cerr << "Could not open input file." << std::endl;
            delete solver;
            return -1;
        }
        solver->readFile(input_file);
        bool isSat =  solver->solve();
        if (isSat) {
            std::cout << "SAT" << std::endl;
            delete solver;
            return 10;
        } else {
            std::cout << "UNSAT" << std::endl;
            delete solver;
            return 20;
        }
    } else {
        solver->runTests();
        //solver.setTest2Formula();
        //std::cout << solver.solve() << std::endl;
    }
    delete solver;
    return 0;
}