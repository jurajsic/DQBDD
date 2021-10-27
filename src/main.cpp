/*
 * This file is part of DQBDD.
 *
 * Copyright 2020, 2021 Juraj Síč
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
#include <algorithm>
#include <stdexcept>

#include <cxxopts.hpp>

#include "hqspreinterface.hpp"
#include "dqdimacsparser.hpp"
#include "prenexdqcirparser.hpp"

enum ReturnCode {
    SAT = 10,
    UNSAT = 20,
    SATPRE = 10, // solved by HQSpre
    UNSATPRE = 20, // solved by HQSpre
    UNKNOWN = 0
};

std::string getLowercaseFileType(std::string path) {
    auto sep = path.find_last_of("\\/");
    if (sep != std::string::npos) {
        path = path.substr(sep + 1, path.size() - sep - 1);
    }

    auto dot = path.find_last_of(".");
    if (dot != std::string::npos)
    {
        path = path.substr(dot + 1, path.size() - dot - 1);
        transform(path.begin(), path.end(), path.begin(), [](unsigned char c){ return std::tolower(c); });
        return path;
    }
    else
    {
        return "";
    }
}

void checkAndPrintDQCIR(dqbdd::GateParser &parser, const std::unique_ptr<cxxopts::ParseResult> &result) {
    if (result->count("dqcir-output")) {
        std::cout << "Turning into DQCIR format" << std::endl;
        std::string outputFileName = (*result)["dqcir-output"].as<std::string>();
        std::ofstream outputFile(outputFileName);
        if (outputFile.is_open()) {
            parser.printPrenexDQCIR(outputFile);
            // std::cout << "Transforming to NNF" << std::endl;
            // parser.transformToNNF();
            // std::cout << "Transformed to NNF" << std::endl;
            // outputFile << "After NNF" << std::endl;
            // parser.printGates(outputFile);
            // std::cout << "Collapsing gates" << std::endl;
            // parser.collapseGates();
            // std::cout << "Gates collapsed" << std::endl;
            // outputFile << "After collapsing gates" << std::endl;
            // parser.printGates(outputFile);
            outputFile.close();
        } else {
            throw std::runtime_error("Could not open output DQCIR file");
            return;
        }
        std::cout << "Parsed and possibly preprocessed input file was succesfully written into DQCIR output file " << outputFileName << std::endl;
    }
    if (result->count("dqcir-output-cleansed")) {
        std::cout << "Turning into cleansed DQCIR format" << std::endl;
        std::string outputFileName = (*result)["dqcir-output-cleansed"].as<std::string>();
        std::ofstream outputFile(outputFileName);
        if (outputFile.is_open()) {
            parser.printPrenexCleansedDQCIR(outputFile);
            outputFile.close();
        } else {
            throw std::runtime_error("Could not open output DQCIR file");
            return;
        }
        std::cout << "Parsed and possibly preprocessed input file was succesfully written into cleansed DQCIR output file " << outputFileName << std::endl;
    }
}

int main(int argc, char **argv)
{
    // argument parsing
    cxxopts::Options optionsParser("DQBDD", "A DQBF solver using BDDs.");
    optionsParser.add_options()
        ("h,help", "Print usage")
        ("v,version", "Print the version number")
        ("l,localise", "Use quantifier tree with localisation of quantifiers", cxxopts::value<int>()->default_value("1"))
        ("p,preprocess", "Use preprocessing (only for (DQ)DIMACS files)", cxxopts::value<int>()->default_value("1"))
        ("e,elimination-choice", "Decide what to eliminate on each level of quantifier tree during transformation to formula", cxxopts::value<int>()->default_value("1"))
        ("u,uvar-choice", "The heuristics by which the next universal variable for elimination is chosen", cxxopts::value<int>()->default_value("0"))
        ("d,dyn-reordering", "Allow dynamic reordering of variables in BDDs", cxxopts::value<int>()->default_value("1"))
        ("force-filetype", "Forces the filetype (0 - (DQ)DIMACS, 1 - (D)QCIR)", cxxopts::value<int>())
        ("dqcir-output", "Writes parsed (and possibly preprocessed) input file into given output file in DQCIR format without solving", cxxopts::value<std::string>())
        ("dqcir-output-cleansed", "Writes parsed (and possibly preprocessed) input file into given output file in cleansed DQCIR format without solving", cxxopts::value<std::string>())
        ("f,file","(DQ)DIMACS/(D)QCIR file to solve", cxxopts::value<std::string>())
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
        std::cout << "DQBDD 1.2" << std::endl;
        return 0;
    }
    
    if (!result->count("file")) {
        std::cerr << "No file specified, try 'DQBDD --help' for more info" << std::endl;
        return -1;
    }

    std::string fileName = (*result)["file"].as<std::string>();
    bool localise = (*result)["localise"].as<int>();
    bool preprocess = (*result)["preprocess"].as<int>();
    bool dynReorder = (*result)["dyn-reordering"].as<int>();
    dqbdd::Options options;
    options.treeElimChoice = static_cast<dqbdd::TreeElimChoice>((*result)["elimination-choice"].as<int>());
    options.uVarElimChoice = static_cast<dqbdd::UnivVarElimChoice>((*result)["uvar-choice"].as<int>());
    
    // fileType = 0 - (DQ)DIMACS, filetype = 1 - (D)QCIR
    int fileType = 0;
    if (!result->count("force-filetype")) { // try to guess filetype automatically
        std::string fileExt = getLowercaseFileType(fileName);
        if (fileExt == "dimacs" || fileExt == "qdimacs" || fileExt == "dqdimacs") {
            fileType = 0;
        } else if (fileExt == "qcir" || fileExt == "dqcir") {
            fileType = 1;
        } else {
            std::cout << "The filetype could not be determined, defaulting to (DQ)DIMACS" << std::endl;
            fileType = 0;
        }
    } else {
        fileType = (*result)["force-filetype"].as<int>();
    }

    Cudd mgr;
    if (dynReorder) {
        mgr.AutodynEnable(Cudd_ReorderingType::CUDD_REORDER_SIFT);
    } else {
        mgr.AutodynDisable();
    }

    dqbdd::QuantifiedVariablesManager qvMgr(options);

    if (result->count("hqspre-dqcir-output")) {
        dqbdd::HQSPreInterface hqspreParser(mgr, qvMgr);
        std::cout << "Starting HQSpre" << std::endl;
        hqspreParser.parse(fileName);
        std::cout << "Turning into DQCIR format" << std::endl;
        auto outputFileName = fileName.substr(0, fileName.size()-9) + ".dqcir";
        std::ofstream outputFile(outputFileName);
        if (outputFile.is_open()) {
            hqspreParser.printPrenexDQCIR(outputFile);
            outputFile.close();
        } else {
            std::cerr << "Could not open output file" << std::endl;
            return -1;
        }
        std::cout << "Input file was successfully preprocessed and tranformed into DQCIR file: " << outputFileName << std::endl;
        return 0;
    }

    dqbdd::Formula *f = nullptr;

    try {
        std::unique_ptr<dqbdd::GateParser> parser;
        std::cout << "Parsing" << std::endl;
        if (preprocess && fileType == 0) {
            std::unique_ptr<dqbdd::HQSPreInterface> hqspreParser(new dqbdd::HQSPreInterface(mgr, qvMgr));
            std::cout << "Starting HQSpre" << std::endl;
            hqspreParser->parse(fileName);
            if (hqspreParser->getPreprocessorResult() != dqbdd::HQSPreResult::UNKNOWN) {
                std::cout << "Solved by preprocessor" << std::endl;
                checkAndPrintDQCIR(*hqspreParser, result);
                if (hqspreParser->getPreprocessorResult() == dqbdd::HQSPreResult::SAT) {
                    std::cout << "SAT" << std::endl;
                    return ReturnCode::SATPRE;
                } else {
                    std::cout << "UNSAT" << std::endl;
                    return ReturnCode::UNSATPRE;
                }
            } else {
                parser = std::move(hqspreParser);
            }
        } else {
            if (fileType == 0) {
                parser = std::make_unique<dqbdd::DQDIMACSParser>(mgr,qvMgr);
            } else if (fileType == 1) {
                parser = std::make_unique<dqbdd::PrenexDQCIRParser>(mgr,qvMgr);
            } else {
                std::cerr << "Unknown filetype (this should not happen)" << std::endl;
                return -1;
            }
            parser->parse(fileName);
        }
        std::cout << "Parsing finished" << std::endl;

        if (result->count("dqcir-output") || result->count("dqcir-output-cleansed")) {
            checkAndPrintDQCIR(*parser, result);
            return 0;
        }

        if (!localise) {
            std::cout << "Creating BDD formula" << std::endl;
            f = parser->getFormula();
        } else {
            std::cout << "Creating quantifier tree" << std::endl;

            auto qtroot = parser->getQuantifierTree();

            std::cout << "Quantifier tree created with " 
                        << qtroot->getUnivVars().size() << " universal and "
                        << qtroot->getExistVars().size() << " existential variables quantified in it." << std::endl
                        //<< *qtroot << std::endl
                        << "Pushing quantifiers inside" << std::endl;
            //std::cout << qtroot->getUnivVars() << std::endl
            //          << qtroot->getExistVars() << std::endl;

            qtroot->localise( dqbdd::VariableSet{ } );

            std::cout << "Quantifiers pushed inside" << std::endl
                        //<< *qtroot << std::endl
                        << "Creating BDD formula" << std::endl;
            
            f = qtroot->changeToFormula(mgr);
            //std::cout << *f <<std::endl;
        }
    } catch(const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return -1; 
    }

    f->removeUnusedVars();
    std::cout << "BDD formula created" << std::endl;
    f->printStats();
    //std::cout << "Universal variables: " << f->getUnivVars() << std::endl;
    //std::cout << "Existential variables: " << f->getExistVars() << std::endl;
    std::cout << "Eliminating variables in the created formula" << std::endl;
    try {
        f->eliminatePossibleVars();
    } catch(const std::exception &e) {
        std::cerr << e.what() << std::endl;
        delete f;
        return -1;
    }

    ReturnCode rc;
    
    if (f->getMatrix().IsOne()) {
        std::cout << "SAT" << std::endl;
        rc = ReturnCode::SAT;
    } else if (f->getMatrix().IsZero()) {
        std::cout << "UNSAT" << std::endl;
        rc = ReturnCode::UNSAT;
    } else { // this should not be reachable
        std::cout << "UNKNOWN" << std::endl;
        rc = ReturnCode::UNKNOWN;
    }
    
    delete f;
    return rc;
}