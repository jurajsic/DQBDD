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

#ifndef DQBDD_HQSPRE_INTERFACE_HPP
#define DQBDD_HQSPRE_INTERFACE_HPP

#include <memory>

#include "parser.hpp"

class HQSPreInterface : public Parser {
private:
    // using pimpl idiom to hide implementation of HQSPre
    class HQSPreFormulaWrapper;
    std::unique_ptr<HQSPreFormulaWrapper> formulaPtr;

    Cudd &mgr;
    QuantifiedVariablesManipulator DQBFPrefix;
public:
    HQSPreInterface(Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    /**
     * @brief Parses a file in DQDIMACS format and runs HQSpre preprocessor
     * 
     * @param fileName name of the file to parse
     */
    bool parse(std::string fileName);
    Formula* getFormula();
    QuantifierTreeNode* getQuantifierTree();
};

#endif