/*
 * This file is part of DQBDD.
 *
 * Copyright 2021 Juraj Síč
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
#include <algorithm>

// for warnings
#include <iostream>

#include "prenexdqcirparser.hpp"
#include "dqbddexceptions.hpp"

namespace dqbdd {

PrenexDQCIRParser::PrenexDQCIRParser(Cudd &mgr, QuantifiedVariablesManager &qvmgr) : GateParser(mgr, qvmgr) {}

void PrenexDQCIRParser::parse(std::string fileName) {
    std::ifstream inputFile(fileName);

    if (!inputFile.is_open()) {
        std::string errorMes = "Could not open file '";
        errorMes += fileName + "'.";
        throw dqbddException(errorMes);
    }

    std::string line;
    bool prefixFinished = false;
    GateLiteral outputGate;
    bool firstLineParsed = false;
    bool isCleansed = false;
    unsigned long maximumAllowedGateID;

    unsigned long maxGateID = 0;
    std::unordered_map<std::string, unsigned long> gateStringToID;
    auto getIDFromGateString = [&isCleansed, &maximumAllowedGateID, &gateStringToID, &maxGateID](std::string gateString)->unsigned long {
        if (isCleansed) { // for cleansed (D)QCIR we just return the number on the input
            unsigned long id = std::stoul(gateString);
            if (id > maximumAllowedGateID) {
                std::cerr << "WARNING: Variable or gate " << id << " was found during parsing (D)QCIR file, which is larger than the allowed maximum from the first line" << std::endl;
            }
            return id;
        } else if (gateStringToID.count(gateString) > 0) {
            return gateStringToID[gateString];
        } else {
            ++maxGateID;
            gateStringToID[gateString] = maxGateID;
            return maxGateID;
        }
    };

    auto getLiteralFromString = [&getIDFromGateString](std::string LiteralStr)->GateLiteral {
        if (LiteralStr[0] == '-') {
            return GateLiteral(false, getIDFromGateString(LiteralStr.substr(1))); 
        } else {
            return GateLiteral(true, getIDFromGateString(LiteralStr));
        }
        // long litNum = std::stol(LiteralStr);
        // if (litNum < 0) {
        //     return GateLiteral(false, -litNum);
        // } else {
        //     return GateLiteral(true, litNum);
        // }
    };

    while(std::getline(inputFile, line)) {
        if (firstLineParsed && (line[0] == '#' || line == "")) {
            continue;
        }

        std::replace_if(line.begin(), line.end(), [](char c) {return (c == '(' || c == ')' || c == ',');}, ' ');
        std::istringstream streamLine(line);
        std::string token;
        streamLine >> token;

        if (!firstLineParsed) {
            if (token != "#QCIR-14") {
                throw dqbddException("First line of (D)QCIR file should start with '#QCIR-14'");
            } else {
                if (streamLine >> token) {
                    maximumAllowedGateID = std::stoul(token);
                    isCleansed = true;
                }
            }
            firstLineParsed = true;
            continue;
        } else if (!prefixFinished) {  // processing quantifier prefix
            std::transform(token.begin(), token.end(), token.begin(), [](unsigned char c){ return std::tolower(c); });
            if (token == "exists") {
                while (streamLine >> token) {
                    addExistVar(getIDFromGateString(token), true);
                }
            } else if (token == "free") {
                while (streamLine >> token) {
                    addExistVar(getIDFromGateString(token), false);
                }
            } else if (token == "forall") {
                while (streamLine >> token) {
                    addUnivVar(getIDFromGateString(token));
                }
            } else if (token == "depend") {
                streamLine >> token;
                unsigned long existVarID = getIDFromGateString(token);
                std::vector<unsigned long> dependenciesID;
                while (streamLine >> token) {
                    if (token != "0") {
                        dependenciesID.push_back(getIDFromGateString(token));
                    }
                }
                addExistVar(existVarID, dependenciesID);
            } else if (token == "output") {
                streamLine >> token;
                outputGate = getLiteralFromString(token);
                prefixFinished = true;
            } else {
                throw dqbddException("Unexpected token found in the quantifier prefix of the input (D)QCIR file.");
            }
        } else {
            // processing gates
            unsigned long inputGateID = getIDFromGateString(token);
            GateType inputGateType;
            std::vector<GateLiteral> operands;
            
            streamLine >> token;
            if (token != "=") {
                throw dqbddException("Unexpected token found in the gates part of the input (D)QCIR file");
            }
            
            streamLine >> token;
            std::transform(token.begin(), token.end(), token.begin(), [](unsigned char c){ return std::tolower(c); });

            if (token == "and") {
                inputGateType = GateType::AND;
            } else if (token == "or") {
                inputGateType = GateType::OR;
            } else if (token == "ite") {
                inputGateType = GateType::MUX;
            } else if (token == "xor") {
                inputGateType = GateType::XOR;
            } else {
                throw dqbddException("Unexpected token in the gates part of the input (D)QCIR file");
            }
            
            while (streamLine >> token) {
                operands.push_back(getLiteralFromString(token));
            }

            addGate(inputGateID, inputGateType, operands);
        }
    }

    finishedParsing(outputGate.first, outputGate.second);
}

} // namespace dqbdd