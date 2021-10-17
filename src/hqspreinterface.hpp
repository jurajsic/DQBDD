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

#include "gateparser.hpp"

namespace dqbdd {

enum class HQSPreResult {
    UNKNOWN,
    SAT,
    UNSAT,
};

/**
 * @brief Parser of formulas in (DQ)DIMACS format that also uses HQSpre preprocessor
 */
class HQSPreInterface : public GateParser {
private:
    HQSPreResult result = HQSPreResult::UNKNOWN;
public:
    HQSPreInterface(Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    /**
     * @brief Parses a file in DQDIMACS format and runs HQSpre preprocessor
     */
    void parse(std::string fileName) override;
    /**
     * @brief Agter parsing, returns the result of HQSpre preprocessing
     */
    HQSPreResult getPreprocessorResult();
};

} // namespace dqbdd

#endif