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

#ifndef DQBDD_PARSER_HPP
#define DQBDD_PARSER_HPP

#include "dqbddformula.hpp"
#include "quantifiertree.hpp"

namespace dqbdd {

/**
 * @brief Base class for formula parsers to inherit from
 */
class Parser {
public:
    Parser() = default;
    virtual ~Parser() = default;
    virtual void parse(std::string fileName) = 0;
    virtual Formula* getFormula() = 0;
    virtual QuantifierTreeNode* getQuantifierTree() = 0;
};

} // namespace dqbdd

#endif