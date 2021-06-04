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

#include <pair>
#include <vector>

#include "parser.hpp"

namespace dqbdd {

enum class GateType {
    AND,
    OR,
    MUX,
    XOR,
}

class Gate {
private:
    GateType type;
    std::vector<std::shared_ptr<Gate>> operands;
}

using GateLiteral = std::pair<bool, Gate>;

class GateParser : public Parser {
private:
    // TODO gates
}

} // namespace dqbdd

#endif