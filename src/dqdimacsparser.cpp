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

#include <sstream>
#include <fstream>

// for warnings
#include <iostream>

#include "dqdimacsparser.hpp"
#include "dqbddexceptions.hpp"

namespace dqbdd {

DQDIMACSParser::DQDIMACSParser(Cudd &mgr, QuantifiedVariablesManager &qvmgr) : GateParser(mgr, qvmgr) {}

void DQDIMACSParser::parse(std::string fileName) {
    std::string line;
    std::ifstream inputFile(fileName);

    if (!inputFile.is_open()) {
        std::string errorMes = "Could not open file '";
        errorMes += fileName + "'.";
        throw dqbddException(errorMes);
    }

    bool pLineProcessed = false;
    bool prefixFinished = false;
    unsigned long expectedNumOfClauses = 0;
    unsigned long maximumVariable = 0;
    std::string lastToken = "";
    std::vector<GateLiteral> clauses;

    auto varStringToUnsignedLongWithMaximumVariableCheck = [&maximumVariable](std::string varString)->unsigned long {
        unsigned long varID = std::stoul(varString);
        if (varID > maximumVariable) {
            std::cerr << "WARNING: Variable with ID " << varID << " was found during parsing (DQ)DIMACS file, which is larger than the allowed maximum from problem line" << std::endl;
        }
        return varID;
    };

    auto getLiteralFromStr = [&varStringToUnsignedLongWithMaximumVariableCheck](std::string dqdimacsLiteral) {
        bool isNegated = (dqdimacsLiteral[0] == '-');
        if (isNegated) {
            dqdimacsLiteral = dqdimacsLiteral.substr(1);
        }
        unsigned long i = varStringToUnsignedLongWithMaximumVariableCheck(dqdimacsLiteral);
        return GateLiteral(!isNegated, i);
    };

    while(std::getline(inputFile, line)) {
        if (line == "") { // TODO whitespaces maybe also can be on empty line???
            // TODO can we have empty lines? if not warning
            continue;
        }
        std::istringstream streamline(line);
        std::string token;
        streamline >> token;
        std::string currentFirstToken = token;
        if (token == "c") {
            continue;
        } else if (token == "p") {
            if (pLineProcessed) {
                throw dqbddException("There are multiple problem lines (lines starting with 'p') in input (DQ)DIMACS.");
            }
            streamline >> token;
            if (token != "cnf") {
                std::cerr << "WARNING: The problem line (i.e. the first line after comments) in input (DQ)DIMACS should have the form 'p cnf <num> <num>'" << std::endl;
            }
            streamline >> token;
            maximumVariable = std::stoul(token);
            streamline >> token;
            expectedNumOfClauses = std::stoul(token);
            pLineProcessed = true;
        } else if (!pLineProcessed && (token == "a" || token == "e" || token == "d")) {
            throw dqbddException("Input DQDIMACS file is missing the problem line (line starting with 'p').");
        } else if (prefixFinished && (token == "a" || token == "e" || token == "d")) {
            throw dqbddException("Prefix in input DQDIMACS file cannot be between the definition of clauses.");
        } else if (token == "a") {
            if (lastToken == "a") {
                std::cerr << "WARNING: Multiple 'a' lines in input after each other." << std::endl;
            }
            while (streamline >> token) {
                if (token != "0") {
                    addUnivVar(varStringToUnsignedLongWithMaximumVariableCheck(token));
                }
            }
        } else if (token == "e") {
            if (lastToken == "e") {
                std::cerr << "WARNING: Multiple 'e' lines in input after each other." << std::endl;
            }
            while (streamline >> token) {
                if (token != "0") {
                    addExistVar(varStringToUnsignedLongWithMaximumVariableCheck(token), true);
                }
            }    
        } else if (token == "d") {
            streamline >> token;
            unsigned long existVarID = varStringToUnsignedLongWithMaximumVariableCheck(token);
            std::vector<unsigned long> dependenciesID;
            while (streamline >> token) {
                if (token != "0") {
                    dependenciesID.push_back(varStringToUnsignedLongWithMaximumVariableCheck(token));
                }
            }
            addExistVar(existVarID, dependenciesID);
        } else { // parse clause (disjunction of literals)
            prefixFinished = true;
            std::vector<GateLiteral> literals;
            literals.push_back(getLiteralFromStr(token));
            while (streamline >> token) {
                if (token == "0") {
                    continue;
                }
                literals.push_back(getLiteralFromStr(token));
            }
            // the clause if the disjunction of literals
            unsigned long clauseGateID = addGate(GateType::OR, literals);
            clauses.push_back(GateLiteral(true, clauseGateID));
        }
        lastToken = currentFirstToken;
    }

    if (expectedNumOfClauses != clauses.size()) {
        std::cout << "WARNING: Expected number of clauses is different from the real number of clauses in input (DQ)DIMACS file." << std::endl;
    }

    // we add the conjunction of all clauses
    finishedParsing(true, addGate(GateType::AND, clauses));

    return;
}

} // namespace dqbdd