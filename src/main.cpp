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
    SATPRE = 10, // solved by HQSpre
    UNSATPRE = 20, // solved by HQSpre
    UNKNOWN = 0
};

int main(int argc, char **argv)
{
    // argument parsing
    cxxopts::Options optionsParser("DQBDD", "A DQBF solver using BDDs.");
    optionsParser.add_options()
        ("h,help", "Print usage")
        ("v,version", "Print the version number")
        ("l,localise", "Force pushing quantifiers inside the tree")
        ("p,preprocess", "Use preprocessing")
        ("e,elimination-choice", "Decide what to eliminate on each level of quantifier tree during transformation to formula"//\
                0: Do not eliminate anything\
                1: Eliminate only existential variables\
              2: Eliminate all universal variables and possible existential variables"
              , cxxopts::value<int>()->default_value("1"))
        ("u,uvar-choice", "The heuristics by which the next universal variable for elimination is chosen", cxxopts::value<int>()->default_value("0"))
        ("f,file","DQDIMACS file to solve",cxxopts::value<std::string>())
        ;
    optionsParser.parse_positional({"file"});
    optionsParser.positional_help("<input file>");
    
    std::unique_ptr<cxxopts::ParseResult> result;
    try {
        result.reset(new cxxopts::ParseResult(optionsParser.parse(argc, argv)));
    } catch (const cxxopts::OptionParseException& e) {
        std::cerr << "Parsing error: " << e.what() << std::endl;
        return -1;
    }


    if (result->count("help")) {
        std::cout << optionsParser.help() << std::endl;
        return 0;
    }

    if (result->count("version")) {
        std::cout << "DQBDD 0.3" << std::endl;
        return 0;
    }
    
    if (!result->count("file")) {
        std::cerr << "No file specified, try 'DQBDD --help' for more info" << std::endl;
        return -1;
    }

    std::string fileName = (*result)["file"].as<std::string>();
    bool localise = result->count("localise");
    bool preprocess = result->count("preprocess");
    Options options;
    options.treeElimChoice = static_cast<TreeElimChoice>((*result)["elimination-choice"].as<int>());
    options.uVarElimChoice = static_cast<UnivVarElimChoice>((*result)["uvar-choice"].as<int>());


    Cudd mgr;
    //mgr.AutodynDisable();
    QuantifiedVariablesManager qvMgr(options);
    Formula *f = nullptr;
    bool preprocessorSolved = false;

    try {
        std::unique_ptr<Parser> parser;
        std::cout << "Parsing" << std::endl;
        if (preprocess) {
            parser = std::make_unique<HQSPreInterface>(mgr, qvMgr);
            std::cout << "Starting HQSpre" << std::endl;
        } else {
            parser = std::make_unique<DQDIMACSParser>(mgr,qvMgr);
        }
        preprocessorSolved = parser->parse(fileName);
        std::cout << "Parsing finished" << std::endl;
        if (!localise) {
            if (!preprocessorSolved) {
                std::cout << "Creating BDD formula" << std::endl;
            }
            f = parser->getFormula();
        } else {
            if (!preprocessorSolved) {
                std::cout << "Creating quantifier tree" << std::endl;
            }
            auto qtroot = parser->getQuantifierTree();
            if (!preprocessorSolved) {
                std::cout << "Quantifier tree created" << std::endl
                          //<< *qtroot << std::endl
                          << "Pushing quantifiers inside" << std::endl;
                //std::cout << qtroot->getUnivVars() << std::endl
                //          << qtroot->getExistVars() << std::endl;
                qtroot->localise();
                std::cout << "Quantifiers pushed inside" << std::endl
                          //<< *qtroot << std::endl
                          << "Creating BDD formula" << std::endl;
            }
            f = qtroot->changeToFormula(mgr);
            //std::cout << *f <<std::endl;
        }
    } catch(const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return -1; 
    }

    
    ReturnCode rc;

    if (!preprocessorSolved) {
        f->removeUnusedVars();
        std::cout << "BDD formula created" << std::endl;
        f->printStats();
        //std::cout << "Universal variables: " << f->getUnivVars() << std::endl;
        //std::cout << "Existential variables: " << f->getExistVars() << std::endl;
        std::cout << "Eliminating universal variables in the created formula" << std::endl;
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