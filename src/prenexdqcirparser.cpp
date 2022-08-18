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

#include <easylogging++.hpp>

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
    std::string outputToken;
    bool firstLineParsed = false;
    bool isCleansed = false;
    unsigned long maximumAllowedGateID;

    std::unordered_map<std::string, unsigned long> gateStringToID;
    auto getIDFromGateString = [this, &isCleansed, &maximumAllowedGateID, &gateStringToID](std::string gateString)->unsigned long {
        if (isCleansed) { // for cleansed (D)QCIR we just return the number on the input
            unsigned long id = std::stoul(gateString);
            if (id > maximumAllowedGateID) {
                LOG(WARNING) << "Variable or gate " << id << " was found during parsing (D)QCIR file, which is larger than the allowed maximum from the first line";
            }
            return id;
        } else { 
            if (gateStringToID.count(gateString) == 0) {
                gateStringToID[gateString] = getMaxGateID()+1;
            }
            return gateStringToID[gateString];
        }
    };

    auto getLiteralFromString = [&getIDFromGateString](std::string LiteralStr)->GateLiteral {
        if (LiteralStr[0] == '-') {
            return GateLiteral(false, getIDFromGateString(LiteralStr.substr(1))); 
        } else {
            return GateLiteral(true, getIDFromGateString(LiteralStr));
        }
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
            firstLineParsed = true;
            if (token != "#QCIR-G14") {
                 LOG(WARNING) << "First line of (D)QCIR file should start with '#QCIR-G14'";
                 // if the first line does not correctly start with "#QCIR-G14", then if it is not comment, we should probably process it
                 if (line[0] == '#') {
                    continue;
                 }
            } else {
                if (streamLine >> token) { // if #QCIR-G14 is followed by something...
                    maximumAllowedGateID = std::stoul(token); // ...it has to be maximum number of (normal/gate) variables and...
                    isCleansed = true; // ...it also means that (D)QCIR is in cleansed form
                }
                continue;
            }
        }
        
        if (!prefixFinished) {  // processing quantifier prefix
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
                streamLine >> outputToken;
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

    finishedParsing(getLiteralFromString(outputToken));
}

} // namespace dqbdd