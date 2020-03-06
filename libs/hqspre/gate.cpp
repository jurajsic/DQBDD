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

#include "gate.hpp"

#include <algorithm>
#include <functional>
#include <ostream>

#include <easylogging++.hpp>

/**
 * \file gate.cpp
 * \brief Contains the implementation of gate-related functions
 * \author Ralf Wimmer
 * \date 04/2017
 */

namespace hqspre {

/**
 * \brief Prints a gate to the given output stream.
 */
std::ostream&
operator<<(std::ostream& stream, const Gate& gate)
{
    stream << lit2dimacs(gate._output_literal) << " = ";
    switch (gate._type) {
        case GateType::AND_GATE:
            stream << "AND( ";
            break;
        case GateType::XOR_GATE:
            stream << "XOR( ";
            break;
        case GateType::MUX_GATE:
            stream << "MUX( ";
            break;
        default:
            stream << "UNKNOWN( ";
            break;
    }

    for (const Literal lit : gate._input_literals) {
        stream << lit2dimacs(lit) << " ";
    }
    stream << ") clauses";

    for (const auto c_nr : gate._encoding_clauses) {
        stream << " " << c_nr;
    }

    return stream;
}

/**
 * \brief Adds a gate to the database and updates the data structures.
 *
 * \param gate the gate to add
 * \note Gates with no encoding clauses are considered invalid and therefore
 * ignored.
 */
void
GateInfo::addGate(const Gate& gate)
{
    val_assert(!gate._encoding_clauses.empty());
    val_assert(gate._output_literal >= 2);
    val_assert(_gate_output[lit2var(gate._output_literal)] == false);

    _gates.emplace_back(gate);
    _gate_output[lit2var(gate._output_literal)] = true;
    for (const Literal lit : gate._input_literals) {
        ++_gate_input[lit2var(lit)];
    }

    for (const ClauseID c_nr : gate._encoding_clauses) {
        val_assert(c_nr < _clause_modified.size());
        val_assert(c_nr < _gate_clause.size());
        val_assert(!_gate_clause[c_nr]);

        _clause_modified[c_nr] = false;
        _gate_clause[c_nr]     = true;
    }
}

/**
 * \brief Adds a gate to the database and updates the data structures.
 *
 * \param gate the gate to add
 * \note Gates with no encoding clauses are considered invalid and therefore
 * ignored.
 */
void
GateInfo::addGate(Gate&& gate)
{
    val_assert(!gate._encoding_clauses.empty());
    val_assert(gate._output_literal >= 2);
    val_assert(_gate_output[lit2var(gate._output_literal)] == false);

    _gates.push_back(std::move(gate));
    _gate_output[lit2var(gate._output_literal)] = true;
    for (const Literal lit : gate._input_literals) {
        ++_gate_input[lit2var(lit)];
        val_assert(_gate_input[lit2var(lit)] >= 1);
    }

    for (const ClauseID c_nr : gate._encoding_clauses) {
        val_assert(c_nr < _clause_modified.size());
        val_assert(c_nr < _gate_clause.size());
        val_assert(!_gate_clause[c_nr]);

        _clause_modified[c_nr] = false;
        _gate_clause[c_nr]     = true;
    }
}

/**
 * \brief Returns the number of valid gates in the database.
 *
 * A gate is valid, if _encoding_clauses is not empty, none of its encoding
 * clauses has been modified, and the output literal is different from 0.
 */
std::size_t
GateInfo::numValidGates() const
{
    return std::count_if(_gates.cbegin(), _gates.cend(), [this](const Gate& g) -> bool { return gateValid(g); });
}

/**
 * \brief Makes a gate in the database invalid.
 *
 * It is still present in the list of gates, but not considered anymore
 * for the remaining data structures.
 */
void
GateInfo::invalidateGate(Gate& g)
{
    if (!gateValid(g)) return;

    _gate_output[lit2var(g._output_literal)] = false;

    for (const Literal lit : g._input_literals) {
        val_assert(_gate_input[lit2var(lit)] >= 1);
        --_gate_input[lit2var(lit)];
    }

    for (const auto c_nr : g._encoding_clauses) {
        _gate_clause[c_nr] = false;
    }

    g._encoding_clauses.clear();
    g._input_literals.clear();
    g._output_literal = 0;
}

/**
 * \brief Checks if a gate is still valid.
 *
 * In order to be valid, all encoding claues must be untouched.
 * Additionally the vector of encoding clauses must not be empty
 * and the output literal not 0.
 * \param g the gate to be checked. Note that the gate must be in the gate
 * database. \return true iff the gate is valid
 */
bool
GateInfo::gateValid(const Gate& g) const
{
    if (g._encoding_clauses.empty() || g._output_literal == 0) return false;

    for (const ClauseID c_nr : g._encoding_clauses) {
        val_assert(c_nr < _clause_modified.size());
        if (_clause_modified[c_nr]) return false;
    }

    return true;
}

/**
 * \brief Updates the database by removing all invalid gates.
 */
void
GateInfo::update()
{
    // Remove invalid gates
    std::size_t finished = 0;
    for (std::size_t curr_pos = 0; curr_pos < _gates.size(); ++curr_pos) {
        if (gateValid(_gates[curr_pos])) {
            if (curr_pos > finished) {
                _gates[finished] = std::move(_gates[curr_pos]);
            }
            ++finished;
        }
    }
    _gates.resize(finished);

    // Update the vectors with gate input and output information
    std::fill(_gate_input.begin(), _gate_input.end(), 0u);
    std::fill(_gate_output.begin(), _gate_output.end(), false);
    std::fill(_gate_clause.begin(), _gate_clause.end(), false);
    std::fill(_clause_modified.begin(), _clause_modified.end(), false);

    for (const auto& gate : _gates) {
        _gate_output[lit2var(gate._output_literal)] = true;
        for (const Literal input_lit : gate._input_literals) {
            ++_gate_input[lit2var(input_lit)];
        }
        for (const ClauseID c_nr : gate._encoding_clauses) {
            _gate_clause[c_nr] = true;
        }
    }
}

/**
 * \brief Brings the gates in the database into a topological order from the
 * inputs to the ouputs.
 *
 * \note It does not remove the invalid gates. This must be done by the user by
 * calling GateInfo::update() before sorting (if necessary).
 */
void
GateInfo::sort()
{
    if (_gates.empty()) return;

    // create dependency graph
    const std::size_t num_vars = _gate_input.size();

    std::vector<std::vector<Variable>> gateDep(num_vars);
    for (const auto& gate : _gates) {
        if (!gateValid(gate)) continue;
        const Variable out_var = lit2var(gate._output_literal);
        for (const Literal in_lit : gate._input_literals) {
            gateDep[out_var].push_back(lit2var(in_lit));
        }
    }

    std::vector<int> which_gate(num_vars, -1);
    for (std::size_t i = 0; i < _gates.size(); ++i) {
        which_gate[lit2var(_gates[i]._output_literal)] = static_cast<int>(i);
    }

    std::vector<bool> processed(num_vars, false);
    for (Variable var = 1; var < num_vars; ++var) {
        processed[var] = !_gate_output[var];
    }

    // dfs is a recursive function which performs depth-first search on the
    // graph defined by the gate dependencies.
    std::vector<Variable> order;
    order.reserve(_gates.size());

    std::function<void(Variable)> dfs = [&gateDep, &processed, &dfs, &order](Variable var) -> void {
        processed[var] = true;
        for (Variable dep : gateDep[var]) {
            if (!processed[dep]) dfs(dep);
        }
        order.push_back(var);
    };  // end lambda function

    for (Variable var = 1; var < num_vars; ++var) {
        if (!processed[var] && _gate_input[var] == 0 && _gate_output[var]) dfs(var);
    }

    std::vector<Gate> result;
    result.reserve(_gates.size());
    for (const Variable var : order) {
        val_assert(which_gate[var] != -1);
        result.push_back(std::move(_gates[which_gate[var]]));
    }

    val_assert(_gates.size() == result.size());

    _gates = std::move(result);
}

/**
 * \brief Removes all information from the database.
 */
void
GateInfo::clear()
{
    _gates.clear();
    _gate_input.clear();
    _gate_output.clear();
    _gate_clause.clear();
    _clause_modified.clear();
}

/**
 * \brief Checks -- as far as possible -- the consistency of the GateInfo data
 * structure. \return true iff the object seems consistent
 */
bool
GateInfo::checkConsistency() const
{
    if (_gate_input.size() != _gate_output.size()) {
        LOG(ERROR) << "gate_input and gate_output have different sizes!";
        return false;
    }

    if (_clause_modified.size() != _gate_clause.size()) {
        LOG(ERROR) << "gate_clause and clause_modified have different sizes!";
        return false;
    }

    for (const auto& g : _gates) {
        for (const ClauseID c_nr : g._encoding_clauses) {
            if (c_nr >= _gate_clause.size()) {
                LOG(ERROR) << "gate_clause too small for encoding clause " << c_nr << '!';
                return false;
            }
        }
    }

    for (std::size_t var = 1; var < _gate_input.size(); ++var) {
        if (_gate_input[var] > _gates.size()) {
            LOG(ERROR) << "gate_input[" << var << "] = " << _gate_input[var] << " is invalid!";
            return false;
        }
    }

    std::vector<unsigned int> go(_gate_output.size(), 0);
    std::vector<unsigned int> cu(_clause_modified.size(), 0);
    for (const auto& g : _gates) {
        if (gateValid(g)) {
            go[lit2var(g._output_literal)] += 1;
            for (const ClauseID c_nr : g._encoding_clauses) {
                cu[c_nr] += 1;
            }
        }
    }

    for (std::size_t var = 1; var < go.size(); ++var) {
        if (go[var] > 1) {
            LOG(ERROR) << "Variable " << var << " is used " << go[var] << " times as gate output!";
            return false;
        }
    }

    for (ClauseID c_nr = 0; c_nr < cu.size(); ++c_nr) {
        if (cu[c_nr] > 1) {
            LOG(ERROR) << "Clause " << c_nr << " is used " << cu[c_nr] << " times as encoding clause!";
            return false;
        }
    }

    return true;
}

}  // end namespace hqspre
