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

 * HQSpre is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with HQSpre. If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <deque>
#include <utility>
#include <vector>

#include "formula.hpp"
#include "gate.hpp"
#include "literal.hpp"
#include "prefix.hpp"

namespace hqspre {

Literal Formula::addAndGate(const Literal input1, const Literal input2)
{
    const Variable output_var = addEVar();
    _prefix->moveToRMB(output_var);

    Clause::ClauseData clause(2, 0);

    // input1, ~output
    clause[0] = input1;
    clause[1] = var2lit(output_var, true);
    const ClauseID c_nr1 = addClause(clause);

    // input2, ~output
    clause[0] = input2;
    const ClauseID c_nr2 = addClause(clause);

    // ~input1, ~input2, output
    clause[0] = negate(input1);
    clause[1] = negate(input2);
    clause.push_back( var2lit(output_var, false) );
    const ClauseID c_nr3 = addClause(clause);

    Gate g(GateType::AND_GATE);
    g._output_literal = var2lit(output_var,false);
    g._input_literals = {input1, input2};
    g._encoding_clauses = {c_nr1, c_nr2, c_nr3};
    _gates.addGate(std::move(g));

    return var2lit(output_var, false);
}


Literal Formula::addNandGate(const Literal input1, const Literal input2)
{
    return negate(addAndGate(input1, input2));
}


Literal Formula::addOrGate(const Literal input1, const Literal input2)
{
    return negate(addAndGate(negate(input1), negate(input2)));
}


Literal Formula::addNorGate(const Literal input1, const Literal input2)
{
    return addAndGate(negate(input1), negate(input2));
}


Literal Formula::addXorGate(const Literal input1, const Literal input2)
{
    const Variable output_var = addEVar();
    _prefix->moveToRMB(output_var);

    Clause::ClauseData clause(3, 0);

    // ~input1, input2, output
    clause[0] = negate(input1);
    clause[1] = input2;
    clause[2] = var2lit(output_var, false);
    const ClauseID c_nr1 = addClause(clause);

    // input1, ~input2, output
    clause[0] = input1;
    clause[1] = negate(input2);
    // output did not change
    const ClauseID c_nr2 = addClause(clause);

    // input1, input2, ~output
    clause[0] = input1;
    clause[1] = input2;
    clause[2] = var2lit(output_var, true);
    const ClauseID c_nr3 = addClause(clause);

    // ~input1, ~input2, ~output
    clause[0] = negate(input1);
    clause[1] = negate(input2);
    // output did not change
    const ClauseID c_nr4 = addClause(clause);

    Gate g(GateType::XOR_GATE);
    g._output_literal = var2lit(output_var, false);
    g._input_literals = { input1, input2 };
    g._encoding_clauses = { c_nr1, c_nr2, c_nr3, c_nr4 };
    _gates.addGate(std::move(g));

    return var2lit(output_var, false);
}


Literal Formula::addMux(const Literal x0, const Literal x1, const Literal select)
{
    const Variable output_var = addEVar();
    _prefix->moveToRMB(output_var);

    Clause::ClauseData clause(3, 0);

    // x0, select, ~output
    clause[0] = x0;
    clause[1] = negate(select);
    clause[2] = var2lit(output_var, true);
    const ClauseID c_nr1 = addClause(clause);

    // ~x0, select, output
    clause[0] = negate(x0);
    clause[1] = select;
    clause[2] = var2lit(output_var, false);
    const ClauseID c_nr2 = addClause(clause);

    // x1, ~select, ~output
    clause[0] = x1;
    clause[1] = negate(select);
    clause[2] = var2lit(output_var, true);
    const ClauseID c_nr3 = addClause(clause);

    // ~x1, ~select, output
    clause[0] = negate(x0);
    clause[1] = negate(select);
    clause[2] = var2lit(output_var, false);
    const ClauseID c_nr4 = addClause(clause);

    Gate g(GateType::MUX_GATE);
    g._output_literal = var2lit(output_var, false);
    g._input_literals = {select, x1, x0};
    g._encoding_clauses = { c_nr1, c_nr2, c_nr3, c_nr4 };
    _gates.addGate(std::move(g));

    return var2lit(output_var, false);
}

Literal Formula::addAndGate(const std::vector<Literal>& inputs)
{
    const Variable output_var = addEVar();
    _prefix->moveToRMB(output_var);

    Clause::ClauseData clause(2, 0);
    Gate g(GateType::AND_GATE);
    g._output_literal = var2lit(output_var, false);
    g._input_literals = inputs;
    g._encoding_clauses.reserve(inputs.size() + 1);

    // input i, negate(output)
    clause[1] = var2lit(output_var, true);
    for (const Literal lit: inputs) {
        clause[0] = lit;
        const ClauseID c_nr = addClause(clause);
        g._encoding_clauses.push_back(c_nr);
    }

    clause.clear();

    // negate(all inputs), output
    for (const Literal lit: inputs) {
        clause.push_back(negate(lit));
    }
    clause.push_back(var2lit(output_var, false));

    const ClauseID c_nr = addClause(clause);
    g._encoding_clauses.push_back(c_nr);

    _gates.addGate(std::move(g));

    return var2lit(output_var, false);
}

Literal Formula::addNandGate(const std::vector<Literal>& inputs)
{
    return negate(addAndGate(inputs));
}


Literal Formula::addOrGate(const std::vector<Literal>& inputs) {
    const Variable output_var = addEVar();
    _prefix->moveToRMB(output_var);

    Clause::ClauseData clause(2, 0);
    Gate g(GateType::AND_GATE); // or-gates are just and-gates with negated in- and outputs!
    g._output_literal = var2lit(output_var, false);
    g._input_literals = inputs;
    g._encoding_clauses.reserve(inputs.size() + 1);

    //input i, negate(output)
    clause[1] = var2lit(output_var, false);
    for (const Literal lit: inputs) {
        clause[0] = lit;
        const ClauseID c_nr = addClause(clause);
        g._encoding_clauses.push_back(c_nr);
    }

    clause.clear();

    //negate(all inputs), output
    for (const Literal lit: inputs) {
        clause.push_back(negate(lit));
    }
    clause.push_back(var2lit(output_var, true));

    const ClauseID c_nr = addClause(clause);
    g._encoding_clauses.push_back(c_nr);

    _gates.addGate(std::move(g));

    return var2lit(output_var, true);
}


Literal Formula::addNorGate(const std::vector<Literal>& inputs)
{
    return negate(addOrGate(inputs));
}

Literal Formula::addXorGate(const std::vector<Literal>& inputs)
{
    std::deque<Literal> queue(inputs.begin(), inputs.end());

    while (queue.size() >= 2) {
        const Literal x1 = queue.front();
        queue.pop_front();
        const Literal x2 = queue.front();
        queue.pop_front();
        const Literal output = addXorGate(x1, x2);
        queue.push_back(output);
    }
    return queue.front();
}

} // end namespace hqspre
