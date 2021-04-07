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

#ifndef DQBDD_PRENEXCLEANSEDQCIRPARSER_HPP
#define DQBDD_PRENEXCLEANSEDQCIRPARSER_HPP

#include <unordered_map>
#include <string>
#include <vector>

#include "dqbddformula.hpp"
#include "quantifiertree.hpp"
#include "parser.hpp"

namespace dqbdd {

/**
 * @brief Parser of formulas in prenex cleansed (D)QCIR format
 * 
 * Prenex cleansed QCIR format is defined at qbflib website, we define prenex cleansed DQCIR format as QCIR 
 * where quantifier of type henkin(v, v1, ..., vn) can be used in any place that exists or forall quantifier 
 * can be used in QCIR. It represents existential variable v with dependency set Dv = {v1, ..., vn}. It is 
 * assumed that v1, ..., vn were already defined as universal variables.
 */
class PrenexCleansedQCIRParser : public Parser {
    Cudd &mgr;
    QuantifiedVariablesManipulator DQBFPrefix;

    // bool says if it is negated or not (true is not, false it is), unsigned long is the gate/var number
    using Literal = std::pair<bool, unsigned long>;
    // bool - true=and, false=or; vector of Literals = operands
    using OperationAndOperands = std::pair<bool, std::vector<Literal>>;

    Literal outputGate;
    std::unordered_map<unsigned long, OperationAndOperands> gates;

    Literal getLiteralFromString(std::string LiteralStr);

    BDD getBDDFromGate(unsigned long gate);
    QuantifierTreeNode* getQTFromGate(unsigned long gate);

public:
    PrenexCleansedQCIRParser(Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    // returns true if resulting formula is trivial (equal to TRUE or FALSE) - can also mean that preprocessor solved
    bool parse(std::string fileName) override;
    Formula* getFormula() override;
    QuantifierTreeNode* getQuantifierTree() override;
};

} // namespace dqbdd

#endif