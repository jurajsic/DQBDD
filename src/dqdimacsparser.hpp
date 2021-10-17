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

#ifndef DQBDD_DQDIMACSPARSER_HPP
#define DQBDD_DQDIMACSPARSER_HPP

#include "gateparser.hpp"

namespace dqbdd {

/**
 * @brief Parser of formulas in (DQ)DIMACS format
 */
class DQDIMACSParser : public GateParser {
public:
    DQDIMACSParser(Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    /**
     * @brief Parses DQBF from file in (DQ)DIMACS format
     */
    void parse(std::string fileName) override;
};

} // namespace dqbdd

#endif