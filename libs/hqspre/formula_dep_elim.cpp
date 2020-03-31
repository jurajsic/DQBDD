/*
 * This file is part of HQSpre.
 *
 * Copyright 2016/17 Ralf Wimmer, Sven Reimer, Paolo Marin, Bernd Becker
 * Albert-Ludwigs-Universitaet Freiburg, Freiburg im Breisgau, Germany
 *
 * HQSpre is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HQSpre is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with HQSpre. If not, see <http://www.gnu.org/licenses/>.
 */

#include <utility>

#include "auxil.hpp"
#include "formula.hpp"
#include "gate.hpp"
#include "literal.hpp"
#include "prefix.hpp"

/**
 * \file formula_dep_elim.cpp
 * \brief Implementation of dependency elimination
 * \author Ralf Wimmer
 * \date 02/2016
 */

namespace hqspre {

/**
 * \brief Eliminates a single dependency from the formula.
 *
 * If the dependency \f$x\in D_y\f$ is to be eliminated, this is done
 * by introducing two additional variables \f$y^0\f$ and \f$y^1\f$ and
 * adding clauses for \f$y \equiv((x\land y^1)\lor(\neg x\land y^0))\f$.
 * Then the new variables can be independent of \f$x\f$.
 * In case that there is an implication \f$(\neg)x\to (\neg)y\f$, we can omit
 * one of the new variables.
 * Afterwards, \f$y\f$ is a Tseitin variable (for a multiplexer output)
 * and depends on all universal variables.
 *
 * \return the pair <y0, y1> (if one of these variables did not need to be
 * created it is replaced by 0). \todo Use a SAT call to test for an implication
 * between \f$x\f$ and \f$y\f$. \todo Check if we can add blocked implications
 * to reduce the number of variable copies.
 */
std::pair<Variable, Variable>
Formula::elimDependency(const Variable x, const Variable y)
{
    if (!depends(x, y)) return std::make_pair(0u, 0u);

    val_assert(_prefix);
    val_assert(_prefix->type() == PrefixType::DQBF);
    val_assert(_dqbf_prefix);
    val_assert(minVarIndex() <= x && x <= maxVarIndex());
    val_assert(minVarIndex() <= y && y <= maxVarIndex());
    val_assert(isUniversal(x));
    val_assert(isExistential(y));
    val_assert_msg(!_gates.isGateOutput(y), "You may not remove a dependency from a gate output!");

    // ITE(a,b,c)=d
    //    is the same as
    // (a,c,~d), (~a,b,~d), (~a,~b,d), (a,~c,d)

    // If there are implications between a literal of x and a literal of y,
    // we can save one copy of y and simplify the clauses.
    // Therefore we check for implications here:
    auto       implied      = hasImplicationTransitive(var2lit(x, false), y);
    const bool x_implies_y  = implied.first;
    const bool x_implies_ny = implied.second;

    implied                  = hasImplicationTransitive(var2lit(x, true), y);
    const bool nx_implies_y  = implied.first;
    const bool nx_implies_ny = implied.second;

    removeDependency(x, y);

    if (x_implies_y && x_implies_ny) {
        // x => y and x => ~y means that ~x is unit (and actually the formula
        // unsatisfiable as x is universal).
        pushUnit(var2lit(x, true), PureStatus::UNIT);
        return std::make_pair(0u, 0u);
    }

    if (nx_implies_y && nx_implies_ny) {
        // ~x => y and ~x => ~y means that x is unit (and actually the formula
        // unsatisfiable as x is universal).
        pushUnit(var2lit(x, false), PureStatus::UNIT);
        return std::make_pair(0u, 0u);
    }

    if (x_implies_y && nx_implies_y) {
        // x => y && ~x => y means that y is always implied.
        pushUnit(var2lit(y, false), PureStatus::UNIT);
        return std::make_pair(0u, 0u);
    }

    if (x_implies_ny && nx_implies_ny) {
        // x => ~y && ~x => ~y means that ~y is always implied.
        pushUnit(var2lit(y, true), PureStatus::UNIT);
        return std::make_pair(0u, 0u);
    }

    if (x_implies_y && nx_implies_ny) {
        // x => y and ~x => ~y means x = y. So replace y by x
        replaceLiteral(var2lit(y, false), var2lit(x, false));
        return std::make_pair(0u, 0u);
    }

    if (x_implies_ny && nx_implies_y) {
        // x => y and ~x => ~y means ~x = y. So replace y by ~x.
        replaceLiteral(var2lit(y, false), var2lit(x, true));
        return std::make_pair(0u, 0u);
    }

    if (x_implies_y) {
        // x -> y. That means, y1 = 1 and the equation simplifies to y = x + ~x*y0
        // Clauses: (x, y0, ~y) (~x, y) (x, y, ~y0)
        const Variable y0 = copyVar(y);

        // make y depend on all universal variables
        _prefix->moveToRMB(y);

        const int c_nr1 = addClause({var2lit(x, false), var2lit(y0, false), var2lit(y, true)});
        const int c_nr2 = addClause({var2lit(x, true), var2lit(y, false)});
        const int c_nr3 = addClause({var2lit(x, false), var2lit(y, false), var2lit(y0, true)});

        /*
                Gate g(GateType::MUX_GATE);
                g._output_literal = var2lit(y, false);
        //TODO: Gate with constant input
                g._input_literals = { var2lit(x, false), 1, var2lit(y0, false) }; //
        select, in1, in0 g._encoding_clauses.reserve(3); if (c_nr1 >= 0)
        g._encoding_clauses.push_back(c_nr1); if (c_nr2 >= 0)
        g._encoding_clauses.push_back(c_nr2); if (c_nr3 >= 0)
        g._encoding_clauses.push_back(c_nr3); _gates.addGate(std::move(g));
        */

        return std::make_pair(y0, 0u);

    } else if (x_implies_ny) {
        // x -> ~y. That means, y1 = 0 and the equation simplifies to y = ~x*y0,
        // i.e., an AND-gate Clauses: (~x, ~y) (y0, ~y) (x, ~y0, y)
        const Variable y0 = copyVar(y);

        // make y depend on all universal variables
        _prefix->moveToRMB(y);

        const int c_nr1 = addClause({var2lit(x, true), var2lit(y, true)});
        const int c_nr2 = addClause({var2lit(y0, false), var2lit(y, true)});
        const int c_nr3 = addClause({var2lit(x, false), var2lit(y0, true), var2lit(y, false)});

        if (c_nr1 >= 0 && c_nr2 >= 0 && c_nr3 >= 0) {
            Gate g(GateType::AND_GATE);
            g._output_literal = var2lit(y, false);
            g._input_literals = {var2lit(y0, false), var2lit(x, true)};
            g._encoding_clauses.reserve(3);
            g._encoding_clauses.push_back(c_nr1);
            g._encoding_clauses.push_back(c_nr2);
            g._encoding_clauses.push_back(c_nr3);
            _gates.addGate(std::move(g));
        }

        return std::make_pair(y0, 0u);

    } else if (nx_implies_y) {
        // ~x -> y. That means, y0 = 1 and the equation simplifies to y = x*y1 + ~x
        // Clauses: (~x, y1, ~y) (x, y) (~x, y, ~y1)
        const Variable y1 = copyVar(y);

        // make y depend on all universal variables
        _prefix->moveToRMB(y);

        const int c_nr1 = addClause({var2lit(x, true), var2lit(y1, false), var2lit(y, true)});
        const int c_nr2 = addClause({var2lit(x, false), var2lit(y, false)});
        const int c_nr3 = addClause({var2lit(x, true), var2lit(y, false), var2lit(y1, true)});

        /*
                Gate g(GateType::MUX_GATE);
                g._output_literal = var2lit(y, false);
                //TODO: Gate with constant input
                g._input_literals = { var2lit(x, false), var2lit(y1, false), 1 };
                g._encoding_clauses.reserve(3);
                if (c_nr1 >= 0) g._encoding_clauses.push_back(c_nr1);
                if (c_nr2 >= 0) g._encoding_clauses.push_back(c_nr2);
                if (c_nr3 >= 0) g._encoding_clauses.push_back(c_nr3);
                _gates.addGate(std::move(g));
        */

        return std::make_pair(0u, y1);

    } else if (nx_implies_ny) {
        // ~x -> ~y. That means, y0 = 0 and the quation simplifies to y = x*y1,
        // i.e., an AND-gate Clauses: (x, ~y) (y1, ~y) (~x, ~y1, y)
        const Variable y1 = copyVar(y);

        // make y depend on all universal variables
        _prefix->moveToRMB(y);

        const int c_nr1 = addClause({var2lit(x, false), var2lit(y, true)});
        const int c_nr2 = addClause({var2lit(y1, false), var2lit(y, true)});
        const int c_nr3 = addClause({var2lit(x, true), var2lit(y1, true), var2lit(y, false)});

        if (c_nr1 >= 0 && c_nr2 >= 0 && c_nr3 >= 0) {
            Gate g(GateType::AND_GATE);
            g._output_literal = var2lit(y, false);
            g._input_literals = {var2lit(y1, false), var2lit(x, false)};
            g._encoding_clauses.reserve(3);
            g._encoding_clauses.push_back(c_nr1);
            g._encoding_clauses.push_back(c_nr2);
            g._encoding_clauses.push_back(c_nr3);
            _gates.addGate(std::move(g));
        }

        return std::make_pair(0u, y1);

    } else {
        // no implications -> proper multiplexer y = x*y1 + ~x*y0
        const Variable y0 = copyVar(y);
        const Variable y1 = copyVar(y);

        // make y depend on all universal variables
        _prefix->moveToRMB(y);

        const int c_nr1 = addClause({var2lit(x, false), var2lit(y0, false), var2lit(y, true)});
        const int c_nr2 = addClause({var2lit(x, false), var2lit(y0, true), var2lit(y, false)});
        const int c_nr3 = addClause({var2lit(x, true), var2lit(y1, false), var2lit(y, true)});
        const int c_nr4 = addClause({var2lit(x, true), var2lit(y1, true), var2lit(y, false)});

        if (c_nr1 >= 0 && c_nr2 >= 0 && c_nr3 >= 0 && c_nr4 >= 0) {
            Gate g(GateType::MUX_GATE);
            g._output_literal = var2lit(y, false);
            g._input_literals = {var2lit(x, false), var2lit(y1, false), var2lit(y0, false)};  // select, input1, input0
            g._encoding_clauses.reserve(4);
            g._encoding_clauses.push_back(c_nr1);
            g._encoding_clauses.push_back(c_nr2);
            g._encoding_clauses.push_back(c_nr3);
            g._encoding_clauses.push_back(c_nr4);
            _gates.addGate(std::move(g));
        }

        return std::make_pair(y0, y1);
    }

    return std::make_pair(0u, 0u);
}

}  // end namespace hqspre
