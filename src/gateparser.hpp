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
#include <list>

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

    GateLiteral outputGateLiteral;
    std::unordered_map<unsigned long, Gate> gateIDToGate;

    // ordered list of gates ID based on when they were added to the parser
    std::list<unsigned long> gateInputOrder;

    unsigned long maxGateID = 0;

    // true if parse(fileName) was called and we have not transformed formula either with getForumla() or getQuantifierTree()
    bool isFormulaParsed = false;

    // TODO description
    unsigned long addNewGateAtPositionWithoutChecks(GateType type, const std::vector<GateLiteral> &operands, const std::list<unsigned long>::iterator &position);

    // TODO description
    void removeMUXAndXORGates();

    // TODO describe - fucks up order of gateIDToGate, can only be used if we are not plannign to print after
    void transformToNNF();
    std::unordered_map<unsigned long, unsigned long> gateIDToNegatedGateID;
    // TODO describe - fucks up order of gateIDToGate, can only be used if we are not plannign to print after
    unsigned long getNegatedGateID(unsigned long gateIDToNegate);
    // TODO describe - fucks up order of gateIDToGate, can only be used if we are not plannign to print after
    void pushNegation(GateLiteral &gateLiteralToPushNegation);

    // TODO description
    void printPrefixAndGates(std::ostream &output);

protected:
    QuantifiedVariablesManipulator DQBFPrefix;

    // TODO description
    void addGate(unsigned long gateID, GateType type);
    // TODO description - operands should be already added gates, or they should be var (implicitly existential without dependencies, if they are not in DQBFprefix)
    void addGate(unsigned long gateID, GateType type, const std::vector<GateLiteral> &operands);
    // TODO description
    void finishedParsing(bool outputGateNegation, unsigned long outputGateID);

public:
    // TODO description
    GateParser(Cudd &mgr, QuantifiedVariablesManager &qvmgr);

    virtual ~GateParser() = default;

    // TODO return should be void, HQSpreparser can have its own member bool preprocessorSolved
    // TODO explain that this function should parse formula from some file by adding gates using addGate and at the end finishedParsing should be called with outputGate
    virtual bool parse(std::string fileName) = 0;
    /**
     * @brief Get DQBF from formula saved as gates in this parser
     * 
     * Can be only used after parsing a formula saved in some file with parse(fileName) and by using
     * this method, the parsed formula is deleted from this parser
     */
    Formula* getFormula() override;
    /**
     * @brief Get DQBF as a quantifier tree from formula saved as gates in this parser
     * 
     * Can be only used after parsing a formula saved in some file with parse(fileName) and by using
     * this method, the parsed formula is deleted from this parser
     */
    QuantifierTreeNode* getQuantifierTree() override;

    // TODO description
    void printPrenexDQCIR(std::ostream &output);
    // TODO description - important: removes MUX and XOR gates!
    void printPrenexCleansedDQCIR(std::ostream &output);

    // TODO description
    void clearParser();
};

} // namespace dqbdd

#endif