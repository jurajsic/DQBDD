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
    Cudd mgr1;
    //mgr.AutodynDisable();
    QuantifiedVariablesManager qvMgr1;
    Variable x1(1, mgr1);
    Variable x2(2, mgr1);
    Variable x3(3, mgr1);
    Variable y1(11, mgr1);
    Variable y2(12, mgr1);
    QuantifiedVariablesManipulator prefix(qvMgr1);
    prefix.addUnivVar(x1);
    prefix.addUnivVar(x2);
    prefix.addUnivVar(x3);
    prefix.addExistVar(y1, VariableSet{x1});
    prefix.addExistVar(y2);
    
    QuantifierTreeFormula *qtx1 = new QuantifierTreeFormula(mgr1, qvMgr1);
    qtx1->setMatrix(x1);
    QuantifierTreeFormula *qtx21 = new QuantifierTreeFormula(mgr1, qvMgr1);
    qtx21->setMatrix(x2);
    QuantifierTreeFormula *qtx22 = new QuantifierTreeFormula(mgr1, qvMgr1);
    qtx22->setMatrix(x2);
    QuantifierTreeFormula *qtx23 = new QuantifierTreeFormula(mgr1, qvMgr1);
    qtx23->setMatrix(!x2);
    QuantifierTreeFormula *qtx3 = new QuantifierTreeFormula(mgr1, qvMgr1);
    qtx3->setMatrix(x3);
    QuantifierTreeFormula *qty11 = new QuantifierTreeFormula(mgr1, qvMgr1);
    qty11->setMatrix(y1);
    QuantifierTreeFormula *qty12 = new QuantifierTreeFormula(mgr1, qvMgr1);
    qty12->setMatrix(y1);
    QuantifierTreeFormula *qty21 = new QuantifierTreeFormula(mgr1, qvMgr1);
    qty21->setMatrix(y2);
    QuantifierTreeFormula *qty22 = new QuantifierTreeFormula(mgr1, qvMgr1);
    qty22->setMatrix(!y2);
    auto qtttt = new QuantifierTree(false, std::list<QuantifierTreeNode*>{
        new QuantifierTree(true, std::list<QuantifierTreeNode*>{qtx1,qty11},qvMgr1),
        new QuantifierTree(true, std::list<QuantifierTreeNode*>{qtx21,qty12},qvMgr1),
        new QuantifierTree(true, std::list<QuantifierTreeNode*>{qtx22,qty21},qvMgr1),
        new QuantifierTree(true, std::list<QuantifierTreeNode*>{qtx23,qty22,qtx3},qvMgr1),
    },prefix);
    prefix.clear();
    
    std::cout << qtttt->getSupportSet() << std::endl;
    std::cout << qtttt->getUVarsOutsideThisSubtree() << std::endl;
    std::cout << qtttt->getUVarsSupportSet() << std::endl;
    std::cout << *qtttt << std::endl;
    qtttt->localise();
    std::cout << *qtttt << std::endl;
    auto fffff = qtttt->changeToFormula(mgr1);
    std::cout << *fffff << std::endl;
    return 0;

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
                          << *qtroot << std::endl
                          << "Pushing quantifiers inside" << std::endl;
                qtroot->localise();
                std::cout << "Quantifiers pushed inside" << std::endl
                          << *qtroot << std::endl
                          << "Creating BDD formula" << std::endl;
            }
            f = qtroot->changeToFormula(mgr);
            std::cout << *f <<std::endl;
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
        std::cout << "Universal variables: " << f->getUnivVars() << std::endl;
        std::cout << "Existential variables: " << f->getExistVars() << std::endl;
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