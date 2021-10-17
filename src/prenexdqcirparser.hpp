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

#ifndef DQBDD_PRENEXDQCIRPARSER_HPP
#define DQBDD_PRENEXDQCIRPARSER_HPP

#include <unordered_map>
#include <string>
#include <vector>

#include "dqbddformula.hpp"
#include "quantifiertree.hpp"
#include "gateparser.hpp"

namespace dqbdd {

/**
 * @brief Parser of formulas in prenex (D)QCIR format
 * 
 * Prenex QCIR format is defined at http://www.qbflib.org/qcir.pdf, we define prenex DQCIR format as QCIR 
 * where quantifier of type depend(v, v1, ..., vn) can be used in any place that exists or forall quantifier 
 * can be used in QCIR. It represents existential variable v with dependency set Dv = {v1, ..., vn}. It is 
 * assumed that v1, ..., vn were already defined as universal variables (i.e. with forall quantifier).
 */
class PrenexDQCIRParser : public GateParser {
public:
    PrenexDQCIRParser(Cudd &mgr, QuantifiedVariablesManager &qvmgr);
    /**
     * @brief Parses DQBF from file in prenex DQCIR format
     */
    void parse(std::string fileName) override;
};

} // namespace dqbdd

#endif