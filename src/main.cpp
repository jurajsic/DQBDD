/*
 * This file is part of DQBDD.
 *
 * Copyright 2020-2022 Juraj Síč
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

#include <easylogging++.hpp>
#include <cxxopts.hpp>

#include "hqspreinterface.hpp"
#include "dqdimacsparser.hpp"
#include "prenexdqcirparser.hpp"
#include "version.hpp"

INITIALIZE_EASYLOGGINGPP

enum ReturnCode {
    SAT = 10,
    UNSAT = 20,
    SATPRE = 10, // solved by HQSpre
    UNSATPRE = 20, // solved by HQSpre
    UNKNOWN = 0,
    ERROR = -1,
    NOSOLVINGSUCCESS = 0, // finished succesfully without solving (print help message, print DQCIR, etc.)
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
        VLOG(1) << "Turning into DQCIR format";
        std::string outputFileName = (*result)["dqcir-output"].as<std::string>();
        std::ofstream outputFile(outputFileName);
        if (outputFile.is_open()) {
            parser.printPrenexDQCIR(outputFile);
            outputFile.close();
        } else {
            throw std::runtime_error("Could not open output DQCIR file");
            return;
        }
        VLOG(1) << "Parsed and possibly preprocessed input file was succesfully written into DQCIR output file " << outputFileName;
    }
    if (result->count("dqcir-output-cleansed")) {
        VLOG(1) << "Turning into cleansed DQCIR format";
        std::string outputFileName = (*result)["dqcir-output-cleansed"].as<std::string>();
        std::ofstream outputFile(outputFileName);
        if (outputFile.is_open()) {
            parser.printPrenexCleansedDQCIR(outputFile);
            outputFile.close();
        } else {
            throw std::runtime_error("Could not open output DQCIR file");
            return;
        }
        VLOG(1) << "Parsed and possibly preprocessed input file was succesfully written into cleansed DQCIR output file " << outputFileName;
    }
}

int main(int argc, char **argv)
{
    cxxopts::Options optionsParser("DQBDD", "A DQBF solver using BDDs.");
    optionsParser.add_options()
        ("h,help", "Print usage")
        ("V,version", "Print the version number")
        ("v,verbose", "Use verbose output (use multiple -v options, or set --verbose=1-9 for the desired verbose level)", cxxopts::value<int>()->default_value("0")->implicit_value("1"))
        ("l,localise", "Use quantifier tree with localisation of quantifiers", cxxopts::value<int>()->default_value("1"))
        ("p,preprocess", "Use preprocessing (only for (DQ)DIMACS files)", cxxopts::value<int>()->default_value("1"))
        ("e,elimination-choice", "Decide what to eliminate on each level of quantifier tree during transformation to formula", cxxopts::value<int>()->default_value("1"))
        ("u,uvar-choice", "The heuristics by which the next universal variable for elimination is chosen", cxxopts::value<int>()->default_value("0"))
        ("d,dyn-reordering", "Allow dynamic reordering of variables in BDDs", cxxopts::value<int>()->default_value("1"))
        ("force-filetype", "Forces the filetype (0 - (DQ)DIMACS, 1 - (D)QCIR)", cxxopts::value<int>())
        ("dqcir-output", "Writes parsed (and possibly preprocessed) input file into given output file in DQCIR format without solving", cxxopts::value<std::string>())
        ("dqcir-output-cleansed", "Writes parsed (and possibly preprocessed) input file into given output file in cleansed DQCIR format without solving", cxxopts::value<std::string>())
        ("print-to-dot", "Prints created quantifier tree into given output file in DOT format", cxxopts::value<std::string>())
        ("print-to-dot-localise", "Prints quantifier tree after localisation into given output file in DOT format", cxxopts::value<std::string>())
        ("initial-ordering", "Decides the initial ordering of variables", cxxopts::value<int>()->default_value("0")) // TODO write better, right now 0 means initial ordering is kinda random given by parsing (for DQDIMACS, if HQSpre -> first univ vars, then exist vars, if not HQSpre -> based on the number in file, for DQCIR for cleansed -> based on the number in file, for not cleansed -> based on the first occurence), 1 means based on expected elimination, the first to eliminate are on the lowest level
        ("f,file","(DQ)DIMACS/(D)QCIR file to solve", cxxopts::value<std::string>())
        ;
    optionsParser.parse_positional({"file"});
    optionsParser.positional_help("<input file>");
    
    std::unique_ptr<cxxopts::ParseResult> result;
    try {
        result.reset(new cxxopts::ParseResult(optionsParser.parse(argc, argv)));
    } catch (const cxxopts::OptionParseException& e) {
        std::cerr << "ERROR: " << e.what() << ", try 'dqbdd --help' for more info" << std::endl;
        return ReturnCode::ERROR;
    }

    if (result->count("help")) {
        std::cout << optionsParser.help() << std::endl;
        return ReturnCode::NOSOLVINGSUCCESS;
    }

    if (result->count("version")) {
        std::cout << "DQBDD " << dqbdd_VERSION_MAJOR << "." << dqbdd_VERSION_MINOR << std::endl;
        return ReturnCode::NOSOLVINGSUCCESS;
    }
    
    if (!result->count("file")) {
        std::cerr << "ERROR: No file specified, try 'dqbdd --help' for more info" << std::endl;
        return ReturnCode::ERROR;
    }

    // TODO create custom logger for DQBDD
    // configure default configuration for logging
    el::Configurations defConf;
    defConf.setToDefault();
    defConf.setGlobally(el::ConfigurationType::ToFile, std::string("false"));
    defConf.setGlobally(el::ConfigurationType::Format, "%datetime %level [DQBDD] %msg");
    defConf.set(el::Level::Verbose, el::ConfigurationType::Format, "%datetime %level-%vlevel [DQBDD] %msg");
    el::Loggers::setDefaultConfigurations(defConf, true);
    // set the verbosity level
    int verboseLevel = static_cast<dqbdd::TreeElimChoice>((*result)["verbose"].as<int>());
    if (verboseLevel == 1) {
        el::Loggers::setVerboseLevel(result->count("verbose"));
    } else {
        el::Loggers::setVerboseLevel(verboseLevel);
    }
    VLOG(1) << "Verbose level set to " << el::Loggers::verboseLevel();

    for (const cxxopts::KeyValue &optionWithValue : result->arguments()) {
        if (optionWithValue.key() != "verbose") {
            VLOG(2) << "Option " << optionWithValue.key() << " is set to " << optionWithValue.value();  
        }
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
            LOG(WARNING) << "The filetype could not be determined, defaulting to (DQ)DIMACS";
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

    dqbdd::Formula *f = nullptr;

    try {
        std::unique_ptr<dqbdd::GateParser> parser;
        VLOG(1) << "Parsing";
        if (preprocess && fileType == 0) {
            std::unique_ptr<dqbdd::HQSPreInterface> hqspreParser(new dqbdd::HQSPreInterface(mgr, qvMgr));
            VLOG(1) << "Starting HQSpre";
            hqspreParser->parse(fileName);
            if (hqspreParser->getPreprocessorResult() != dqbdd::HQSPreResult::UNKNOWN) {
                VLOG(1) << "Solved by preprocessor";
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
                std::cerr << "ERROR: Unknown filetype (this should not happen)" << std::endl;
                return ReturnCode::ERROR;
            }
            parser->parse(fileName);
        }
        VLOG(1) << "Parsing finished";

        if (result->count("dqcir-output") || result->count("dqcir-output-cleansed")) {
            checkAndPrintDQCIR(*parser, result);
            return ReturnCode::NOSOLVINGSUCCESS;
        }

        if (!localise) {
            VLOG(1) << "Creating BDD formula";
            f = parser->getFormula();
        } else {
            VLOG(1) << "Creating quantifier tree";
            auto qtroot = parser->getQuantifierTree();
            VLOG(1) << "Quantifier tree created"; 
            VLOG(2) << "There are " << qtroot->getUnivVars().size() << " universal and " << qtroot->getExistVars().size() << " existential variables quantified in the prefix of the quantifier tree";
            VLOG(3) << "Univ vars in the prefix of the root: " << qtroot->getUnivVars();
            VLOG(3) << "Exist vars in the prefix of the root: " << qtroot->getExistVars();
            VLOG(9) << "Current formula: " << *qtroot; // for crazy level of verbosity, we print the whole tree (as a formula)

            if (result->count("print-to-dot") > 0) {
                std::string outputFileName = (*result)["print-to-dot"].as<std::string>();
                std::ofstream outputFile(outputFileName);
                if (outputFile.is_open()) {
                    qtroot->printToDot(outputFile);
                } else {
                    std::cerr << "ERROR: Could not open output DOT file" << std::endl;
                    delete qtroot;
                    return ReturnCode::ERROR;
                }
            }

            VLOG(1) << "Pushing quantifiers inside";
            qtroot->localise( dqbdd::VariableSet{ } );
            VLOG(1) << "Quantifiers pushed inside";
            VLOG(2) << "There are " << qtroot->getUnivVars().size() << " universal and " << qtroot->getExistVars().size() << " existential variables quantified in the prefix of the quantifier tree";
            VLOG(3) << "Univ vars in the prefix of the root: " << qtroot->getUnivVars();
            VLOG(3) << "Exist vars in the prefix of the root: " << qtroot->getExistVars();
            VLOG(9) << "Current formula: " << *qtroot; // for crazy level of verbosity, we print the whole tree (as a formula)

            if (result->count("print-to-dot-localise") > 0) {
                std::string outputFileName = (*result)["print-to-dot-localise"].as<std::string>();
                std::ofstream outputFile(outputFileName);
                if (outputFile.is_open()) {
                    qtroot->printToDot(outputFile);
                } else {
                    std::cerr << "ERROR: Could not open output DOT file" << std::endl;
                    delete qtroot;
                    return ReturnCode::ERROR;
                }
            }

            if ((*result)["initial-ordering"].as<int>() == 1) { // TODO create enum for this option, we will probably implement more initial orderings
                VLOG(1) << "Set initial ordering of variables based on the expected elimination order";
                qtroot->getManager()->reorderVars(mgr); // TODO this reordering kinda assumes that uvar-choice==0
            }

            VLOG(1) << "Creating BDD formula from the quantifier tree";
            f = qtroot->changeToFormula(mgr);
        }
    } catch(const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return ReturnCode::ERROR; 
    }

    f->removeUnusedVars();
    VLOG(1) << "Main BDD formula created";
    if (el::Loggers::verboseLevel() > 0) {
        f->printStats();
    }
    VLOG(3) << "Universal variables: " << f->getUnivVars();
    VLOG(3) << "Existential variables: " << f->getExistVars();
    VLOG(9) << "Current formula: " << *f; // for crazy level of verbosity, we print the whole formula
    VLOG(1) << "Eliminating variables in the created formula";
    try {
        f->eliminatePossibleVars();
    } catch(const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        delete f;
        return ReturnCode::ERROR;
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