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
    AND, // can have multiple operands, zero operands represent true
    OR,  // can have multiple operands, zero operands represent false
    MUX, // = ITE, must have 3 operands, MUX(A,B,C) = (A AND B) OR (!A AND C) 
    XOR, // must have 2 operands, A XOR B = (A AND !B) OR (!A AND B)
    VAR, // must have no operands, represents variable
};

// bool represents if the gate is negated (if false, then it is negated), unsigned long is gate id
using GateLiteral = std::pair<bool, unsigned long>;

struct Gate {
    unsigned long ID;
    GateType type;
    std::vector<GateLiteral> operands;
};

/**
 * @brief Base class of DQBF parser with gates
 * Derived parsers must implement the method parse
 */
class GateParser : public Parser {
private:
    Cudd &mgr;
    QuantifiedVariablesManipulator DQBFPrefix;

    GateLiteral outputGateLiteral;
    std::unordered_map<unsigned long, Gate> gateIDToGate;

    // ordered list of gates ID based on when they were added to the parser (used for printing to DQCIR)
    std::list<unsigned long> gateInputOrder;

    unsigned long maxGateID = 0;

    // true if parse(fileName) was called and we have not transformed formula either with getForumla() or getQuantifierTree()
    bool isFormulaParsed = false;

    /**
     * @brief Adds new gate with given type and operands and saves it at given position of gateIDToGate
     * 
     * @param position The position in gateIDToGate before which the newly added gate will be saved
     * @return the ID of the newly added gate 
     */
    unsigned long addNewGateAtPositionWithoutChecks(GateType type, const std::vector<GateLiteral> &operands, const std::list<unsigned long>::iterator &position);

    /**
     * @brief Replaces MUX and XOR gates with combination of AND and OR gates 
     */
    void removeMUXAndXORGates();

    /**
     * @brief Transforms the gates into negation normal form (contains only AND and OR gates and negation is only before VAR gates) 
     * Important: The order of gateInputOrder will be wrong after calling this function, can be therefore used only if we are 
     * "using up" the parser (i.e. in the getQuantifierTree() method)
     */
    void transformToNNF();
    // helping map for transformToNNF(), maps gate IDs to their negated gates' ID
    std::unordered_map<unsigned long, unsigned long> gateIDToNegatedGateID;
    // helping function for transformToNNF(), creates or gets already created negated gate ID of given gate ID (gateInputOrder is wrong after this)
    unsigned long getNegatedGateID(unsigned long gateIDToNegate);
    /* helping function for transformToNNF(), recursively pushes the negations in the given gate literal 
     * (calling this on outputGate will transform the gates into negation normal form, assuming there are no MUX and XOR gates)
     */
    void pushNegation(GateLiteral &gateLiteralToPushNegation);

    /**
     * @brief Collapses multiple AND and OR gates into one
     * The method finds all AND (OR) gates that have AND (OR) gates as operands and collapses them into one,
     * i.e. if we have AND(AND(x1, OR(x2, OR(x1, x3)), AND(x3, x4)), AND(x5)), it collapses it to 
     * AND(x1, OR(x2, x1, x3), x3, x4, x5). It collapses constant values (i.e. AND(x1, AND()) turns to
     * AND(x1) as AND() represents true so it has no impact).
     * 
     * It is assumed that gates are in NNF
     */
    void collapseGates();

    // prints prefix of DQCIR (without output gate)
    void printPrefix(std::ostream &output);

    // prints gates for DQCIR (with output gate of prefix)
    void printGates(std::ostream &output);

protected:
    /**
     * @brief Adds existential variable to prefix
     * Adds existential variable with ID existVarID to the DQBF prefix and a VAR gate with the same ID. Adding VAR gate 
     * for this variable after calling this function will result in thrown exception! Also, adding another variable with
     * the same ID will throw an exception (even if we add it as an existential again). 
     * If dependsOnAllDefinedUnivVars is true, this existential variable will depend on all already added universal variables,
     * otherwise it has an empty dependency set.
     */
    void addExistVar(unsigned long existVarID, bool dependsOnAllDefinedUnivVars = false);
    /**
     * @brief Adds existential variable to prefix with given dependency set
     * Adds existential variable with ID existVarID to the DQBF prefix and a VAR gate with the same ID. Adding VAR gate 
     * for this variable after calling this function will result in thrown exception! Also, adding another variable with
     * the same ID will throw an exception (even if we add it as an existential again). 
     * Newly added existential variable will depend on all universal variables with IDs from dependencySetVarIDs. The
     * variables from the dependency set must have been already added as universal variables with addUnivVar.
     */
    void addExistVar(unsigned long existVarID, const std::vector<unsigned long> &dependencySetVarIDs);
    /**
     * @brief Adds universal variable to prefix
     * Adds universal variable with ID univVarID to the DQBF prefix and a VAR gate with the same ID. Adding VAR gate 
     * for this variable after calling this function will result in thrown exception! Also, adding another variable with
     * the same ID will throw an exception (even if we add it as an universal again).
     */
    void addUnivVar(unsigned long univVarID);

    /**
     * @brief Adds new gate with given ID and type and without operands
     */
    void addGate(unsigned long gateID, GateType type);
    /**
     * @brief Adds new gate with given ID, type and operands
     * The operands should be gate literals with gates which were already added into this parser (except for VAR gates).
     * If an operand is a VAR gate, then it will be added as a new existential variable without dependencies.
     */
    void addGate(unsigned long gateID, GateType type, const std::vector<GateLiteral> &operands);
    /**
     * @brief Adds a new gate with a given type and without operands
     * 
     * @return the ID of the newly added gate 
     */
    unsigned long addGate(GateType type);
    /**
     * @brief Adds new gate with given type and operands
     * The operands should be gate literals with gates which were already added into this parser (except for VAR gates).
     * If an operand is a VAR gate, then it will be added as a new existential variable without dependencies.
     * 
     * @return the ID of the newly added gate 
     */
    unsigned long addGate(GateType type, const std::vector<GateLiteral> &operands);
    /**
     * @brief Method which should be called after parsing is finished with the output gate
     * 
     * @param isOutputGatePositive if false, the output from output gate is negated
     */
    void finishedParsing(bool isOutputGatePositive, unsigned long outputGateID);

public:
    GateParser(Cudd &mgr, QuantifiedVariablesManager &qvmgr);

    /**
     * @brief Parses DQBF formula from file
     * Derived class should implement this method where it should add existential and universal variables to DQBF prefix
     * using addExistVar and addUnivVar methods. It should also add gates using any of the addGate methods. Finally,
     * after parsing is finished, the method finishedParsing should be called by which the output gate can be defined.
     */
    virtual void parse(std::string fileName) override = 0;
    /**
     * @brief Get DQBF from formula saved as gates in this parser
     * 
     * Can be only used after parsing a formula saved in some file with parse(fileName) and by using
     * this method, the parsed formula is deleted from this parser.
     * 
     * @return pointer to resulting formula, needs to be later deleted
     */
    Formula* getFormula() override;
    /**
     * @brief Get DQBF as a quantifier tree from formula saved as gates in this parser
     * Can be only used after parsing a formula saved in some file with parse(fileName) and by using
     * this method, the parsed formula is deleted from this parser. The result should be also later deleted!
     * 
     * @return pointer to resulting quantifier tree, needs to be later deleted
     */
    QuantifierTreeNode* getQuantifierTree() override;

    /**
     * @brief Prints the DQBF saved in this parser as prenex DQCIR format
     */
    void printPrenexDQCIR(std::ostream &output);
    /**
     * @brief Prints the DQBF saved in this parser as prenex cleansed DQCIR format
     * Important: MUX and XOR gates will be replaced with combination of AND and XOR gates.
     */
    void printPrenexCleansedDQCIR(std::ostream &output);

    /**
     * @brief Clears this parser, making it possible to use it for parsing again
     */
    void clearParser();
};

} // namespace dqbdd

#endif