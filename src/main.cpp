/*
 * This file is part of DQBDD.
 *
 * Copyright 2020 Juraj Síč
 *
 * DQBDD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * DQBDD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with DQBDD. If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <iostream>

#include <cxxopts.hpp>

//#include "simplesolver.hpp"
//#include "treesolver.hpp"

#include "HQSpreinterface.hpp"

// TODO add some sort of checking (maybe creating dependency graph -> in manager but what 
// to do with tree) of whether formula is in QBF form??

enum ReturnCode {
    SAT = 10,
    UNSAT = 20,
    SATPRE = 30, // solved by HQSpre
    UNSATPRE = 40, // solved by HQSpre
    UNKNOWN = 50
};

int main(int argc, char **argv)
{
    std::cout << "This is DQBDD version 0.1" << std::endl;

    // argument parsing
    cxxopts::Options options("DQBDD", "A DQBF solver using BDDs.");
    options.add_options()
        ("h,help", "Print usage")
        ("l,localise", "Push quantifiers inside the formula")
        ("r,removal", "The heuristics to decide what to remove in each subformula", cxxopts::value<int>()->default_value("1"))
        ("e,expansion", "The heuristics to decide how to choose the next universal variable for expansion", cxxopts::value<int>()->default_value("1"))
        ("f,file","DQDIMACS file to solve",cxxopts::value<std::string>())
        ;
    options.parse_positional({"file"});
    options.positional_help("<input file>");
    
    cxxopts::ParseResult *result;
    try {
        result = new cxxopts::ParseResult(options.parse(argc, argv));
    } catch (const cxxopts::OptionParseException& e) {
        std::cerr << "Parsing error: " << e.what() << std::endl;
        return -1;
    }

    if (result->count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }
    
    if (!result->count("file")) {
        std::cerr << "No file specified, try 'DQBDD --help' for more info" << std::endl;
        return -1;
    }

    std::string fileName = (*result)["file"].as<std::string>();
    bool localise = result->count("localise");
    // for now, this is doing nothing
    int removeHeuristic = (*result)["removal"].as<int>();
    // for noww this is doing nothing
    int varToExpansionHeuristic = (*result)["expansion"].as<int>();


    Cudd mgr;
    //mgr.AutodynDisable();
    std::ifstream input_file;
    //Solver *solver;
    /*if (argc <= 1) {
        throw "wrong args";
        return -1;
    }*/ /*else {
        if (std::stoi(argw[1]) == 0) {
            solver = new SimpleSolver(mgr);
        } else {
            solver = new TreeSolver(mgr);
        }
    }*/

    

   // if (argc > 2) {

/**************************/
        //delete solver;
        QuantifiedVariablesManager qvMgr;
        Parser *parser = new HQSPreInterface(mgr, qvMgr);
        Formula *f;
        bool preprocessorSolved = false;
        /*std::ofstream statusOutput;
        if (argc > 3) {
            statusOutput.open(argw[3], std::ios::out | std::ios::app);
            statusOutput << argw[2] << ',';
        }*/
        try {
            preprocessorSolved = parser->parse(fileName);
            std::cout << "Parsing finished" << std::endl;
            //statusOutput << 'P';
            if (!localise) {
                f = parser->getFormula();
            } else {
                auto qtroot = parser->getQuantifierTree();
                std::cout << "Created quantifier tree" << std::endl;
                //statusOutput << 'Q';
                f = qtroot->changeToFormula(mgr);
            }
        } catch(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > s) {
            std::cerr << s << std::endl;
            delete parser;
            return -1; 
        }

        
        std::cout << "Created BDD formula" << std::endl;
        //statusOutput << 'F';
        ReturnCode rc;
        delete parser;

        if (!preprocessorSolved) {
            f->eliminatePossibleVars();
            //statusOutput << 'E';
        } else {
            std::cout << "Solved by preprocessor" << std::endl;
        }

        //statusOutput << 'S';
        
        if (f->getMatrix().IsOne()) {
            std::cout << "SAT" << std::endl;
            rc = preprocessorSolved ? ReturnCode::SATPRE : ReturnCode::SAT;
        } else if (f->getMatrix().IsZero()) {
            std::cout << "UNSAT" << std::endl;
            rc = preprocessorSolved ? ReturnCode::UNSATPRE : ReturnCode::UNSAT;
        } else { // this should not be reachable
            std::cout << "UNKNOWN" << std::endl;
            rc = ReturnCode::UNKNOWN;
        }
        
        delete f;
        return rc;
   // }
/**********************/
/*
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
    */
}