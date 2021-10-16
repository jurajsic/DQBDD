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

#ifndef DQBDD_HQSPRE_INTERFACE_HPP
#define DQBDD_HQSPRE_INTERFACE_HPP

// #include <memory>
// #include <ostream>

#include "gateparser.hpp"

namespace dqbdd {

enum class HQSPreResult {
    UNKNOWN,
    SAT,
    UNSAT,
};

/**
 * @brief Parser which uses HQSpre preprocessor to also preprocess the parsed formula
 */
class HQSPreInterface : public GateParser {
private:
    // using pimpl idiom to hide the implementation of HQSPre
    // class HQSPreFormulaWrapper;
    // std::unique_ptr<HQSPreFormulaWrapper> formulaPtr;

    // Cudd &mgr;
    // QuantifiedVariablesManipulator DQBFPrefix;

    HQSPreResult result = HQSPreResult::UNKNOWN;
public:
    HQSPreInterface(Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    /**
     * @brief Parses a file in DQDIMACS format and runs HQSpre preprocessor
     * 
     * @param fileName name of the file to parse
     * @return true if formula was solved by preprocessor
     */
    // TODO change bool to void and add some member to check if solved by preprocessor
    void parse(std::string fileName) override;

    HQSPreResult getPreprocessorResult();
    // Formula* getFormula() override;
    // QuantifierTreeNode* getQuantifierTree() override;
    // ~HQSPreInterface();

    // void turnIntoDQCIR(std::ostream &output);
};

} // namespace dqbdd

#endif