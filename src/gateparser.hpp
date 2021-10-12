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

#ifndef DQBDD_GATEPARSER_HPP
#define DQBDD_GATEPARSER_HPP

#include <utility>
#include <vector>
#include <unordered_map>

#include "parser.hpp"

namespace dqbdd {

enum class GateType {
    AND,
    OR,
    MUX,
    XOR,
    VAR,
};

// bool represents if the gate is negated (if false, then it is negated), num is gate id
using GateLiteral = std::pair<bool, unsigned long>;

struct Gate {
    unsigned long ID;
    GateType type;
    std::vector<GateLiteral> operands;
};

class GateParser : public Parser {
private:
    Cudd &mgr;
    QuantifiedVariablesManipulator DQBFPrefix;

    GateLiteral outputGateLiteral;
    std::unordered_map<unsigned long, Gate> gateIDtoGate;

    // true if parse(fileName) was called and we did not transformed formula either with getForumla() or getQuantifierTree()
    bool isFormulaParsed;

    void clearParser();

    // TODO implement
    void transformToNNF();

protected:
    // TODO implement
    void addGate(unsigned long gateID, GateType type);
    // TODO implement
    void addGate(unsigned long gateID, GateType type, std::vector<unsigned long> &operands);
    // TODO implement
    void finishedParsing(unsigned long outputGateID);

public:
    GateParser(Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    // TODO return should be void, HQSpreparser can have its own member bool preprocessorSolved
    // TODO explain that this function should parse formula from some file by adding gates using addGate and at the end finishedParsing should be called with outputGate
    virtual bool parse(std::string fileName) = 0;
    /**
     * @brief Get DQBF from formula saved as gates in this parser
     * 
     * Can be only used after parsing a formula saved in some file with parse(fileName) and by using
     * this method, the parsed formula is deleted from this parser
     */
    Formula* getFormula();
    /**
     * @brief Get DQBF as a quantifier tree from formula saved as gates in this parser
     * 
     * Can be only used after parsing a formula saved in some file with parse(fileName) and by using
     * this method, the parsed formula is deleted from this parser
     */
    QuantifierTreeNode* getQuantifierTree();

    // TODO implement
    void printPrenexDQCIR(std::ostream &output);
    // TODO implement
    void printPrenexCleansedDQCIR(std::ostream &output);
};

} // namespace dqbdd

#endif