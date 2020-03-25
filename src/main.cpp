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

#include "HQSpreinterface.hpp"
#include "DQDIMACSparser.hpp"

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
    // argument parsing
    cxxopts::Options options("DQBDD", "A DQBF solver using BDDs.");
    options.add_options()
        ("h,help", "Print usage")
        ("v,version", "Print the version number")
        ("l,localise", "Push quantifiers inside the formula")
        ("p,preprocess", "Use preprocessing")
        ("r,removal", "The heuristics to decide what to remove in each subformula", cxxopts::value<int>()->default_value("0"))
        ("e,expansion", "The heuristics to decide how to choose the next universal variable for expansion", cxxopts::value<int>()->default_value("1"))
        ("f,file","DQDIMACS file to solve",cxxopts::value<std::string>())
        ;
    options.parse_positional({"file"});
    options.positional_help("<input file>");
    
    std::unique_ptr<cxxopts::ParseResult> result;
    try {
        result.reset(new cxxopts::ParseResult(options.parse(argc, argv)));
    } catch (const cxxopts::OptionParseException& e) {
        std::cerr << "Parsing error: " << e.what() << std::endl;
        return -1;
    }


    if (result->count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    if (result->count("version")) {
        std::cout << "DQBDD 0.1" << std::endl;
        return 0;
    }
    
    if (!result->count("file")) {
        std::cerr << "No file specified, try 'DQBDD --help' for more info" << std::endl;
        return -1;
    }

    std::string fileName = (*result)["file"].as<std::string>();
    bool localise = result->count("localise");
    bool preprocess = result->count("preprocess");
    // for now, this is doing nothing
    // TODO check if it is not out of range 
    UnivVarElimHeuristic removeHeuristic = static_cast<UnivVarElimHeuristic>((*result)["removal"].as<int>());
    // for noww this is doing nothing
    int varToExpansionHeuristic = (*result)["expansion"].as<int>();


    Cudd mgr;
    //mgr.AutodynDisable();
    QuantifiedVariablesManager qvMgr(removeHeuristic);
    Formula *f = nullptr;
    bool preprocessorSolved = false;
    Parser *parser;
    if (preprocess) {
        parser = new HQSPreInterface(mgr, qvMgr);
    } else {
        parser = new DQDIMACSParser(mgr,qvMgr);
    }

    try {
        preprocessorSolved = parser->parse(fileName);
        std::cout << "Parsing finished" << std::endl;
        if (!localise) {
            f = parser->getFormula();
        } else {
            auto qtroot = parser->getQuantifierTree();
            if (!preprocessorSolved) {
                std::cout << "Created quantifier tree" << std::endl;
                qtroot->localise();
                std::cout << "Pushed quantifiers inside" << std::endl;
            }
            f = qtroot->changeToFormula(mgr);
        }
    } catch(const std::exception &e) {
        std::cerr << e.what() << std::endl;
        delete parser;
        return -1; 
    }

    delete parser;
    
    ReturnCode rc;

    if (!preprocessorSolved) {
        std::cout << "Created BDD formula" << std::endl;
        try {
            f->eliminatePossibleVars();
        } catch(const std::exception &e) {
            std::cerr << e.what() << std::endl;
            delete f;
            return -1;
        }
    } else {
        std::cout << "Solved by preprocessor" << std::endl;
    }

    
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
}