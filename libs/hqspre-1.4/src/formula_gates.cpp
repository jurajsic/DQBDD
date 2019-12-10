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
#include <cstdint>
#include <iterator>
#include <map>
#include <stack>
#include <unordered_map>
#include <utility>
#include <vector>

#include <easylogging++.hpp>

#include "aux.hpp"
#include "clause.hpp"
#include "gate.hpp"
#include "literal.hpp"
#include "timer.hpp"

extern "C" {
#include <picosat.h>
}

#include "formula.hpp"
#include "sat_solver.hpp"

/**
 * \file formula_gates.cpp
 * \brief Implementation of gate detection and structural hashing
 * \author Ralf Wimmer
 * \date 01/2016
 */

namespace hqspre {


/**
 * \brief Checks if a given variable is suited to be a gate output.
 *
 * A variable can only be a gate output if it is existential and depends
 * at least on the variable the gate inputs depend on. This function does
 * not check if all necessary clauses are present.
 * \param output_var the variable that is assumed to be a gate output
 * \param clause a clause that contains the gate inputs and the gate output
 * \return true if the described conditions are satisfied
 * \sa Formula::determineGates()
 */
bool Formula::checkGatePrecondition(const Variable output_var, const Clause& clause) const
{
    if (isUniversal(output_var)) return false;
    for (const auto lit: clause) {
        const auto current_var = lit2var(lit);
        if (current_var == output_var) continue;
        if (!dependenciesSubset(current_var, output_var)) return false;
    }
    return true;
}


namespace internal {

/**
 * \brief Adds the dependencies of a single gate to gateDep
 *
 * This function is used by the gate detection for analyzing
 * whether a gate would create a cycle in the circuit.
 * \sa cyclicDep
 * \sa Formula::determineGates()
 */
template <typename Allocator>
void addDep(std::vector<std::vector<Variable>>& gateDep, const Variable base_var, const std::vector<Variable, Allocator>& parents)
{
#ifndef NDEBUG
    val_assert(gateDep.size() > base_var);
    for (const Variable var: parents) {
        val_assert(gateDep.size() > var);
    }
#endif

    gateDep[base_var].reserve(gateDep[base_var].size() + parents.size());
    std::copy(parents.cbegin(), parents.cend(), std::back_inserter(gateDep[base_var]));
}


/**
 * \brief Returns true if adding another gate would create a cyclic dependency.
 *
 * This function is used by the gate detection for analyzing whether a gate
 * would create a cycle in the circuit.
 * \sa addDep
 * \sa Formula::determineGates
 */
template <typename Allocator>
bool cyclicDep(const std::vector<std::vector<Variable>>& gateDep, const Variable base_var, const std::vector<Variable, Allocator>& parents)
{
#ifndef NDEBUG
    val_assert(gateDep.size() > base_var);
    for (const Variable var: parents) {
        val_assert(gateDep.size() > var);
    }
#endif

    std::stack<Variable> pending;
    std::vector<bool> processed(gateDep.size(), false);

    for (const Variable p: parents) {
        val_assert(p != base_var);
        pending.push(p);
    }

    while (!pending.empty()) {
        const Variable pt = pending.top();
        pending.pop();

        if (!processed[pt]) {
            const auto& dep = gateDep[pt];
            processed[pt] = true;

            for (const Variable p: dep) {
                if (p == base_var) return true;
                else if (!processed[p]) pending.push(p);
            }
        }
    }
    return false;
}


/**
 * \brief Checks if two clauses contain the same variables (with different polarities).
 * \pre Requires c1 and c2 to be sorted.
 * \pre Clauses must not be tautological or contain duplicate literals.
 */
bool sameVars(const Clause& c1, const Clause& c2)
{
    if (c1.size() != c2.size()) return false;

    return std::equal(c1.cbegin(), c1.cend(), c2.cbegin(),
                [](const Literal a, const Literal b) -> bool { return lit2var(a) == lit2var(b); }
            );
}

/**
 * \brief Counts the number of positive literals in a clause.
 */
unsigned long int posCount(const Clause& clause)
{
    return std::count_if(clause.cbegin(), clause.cend(), isPositive);
}

} // end namespace internal


/**
 * \brief Performs gate detection on a quantified formula in CNF.
 *
 * We can detect 2-input XOR and multi-input AND-gates with arbitrarily
 * negated inputs and outputs. The gates are returned in topological order
 * (from inputs to outputs).
 * \sa Formula::substituteGate(const Gate&)
 * \sa Formula::gateDependencies()
 */
bool Formula::determineGates(const bool and_gates, const bool xor_gates, const bool mux_gates, const bool semantic_gates)
{
    using namespace internal;

    if (!_unit_stack.empty()) unitPropagation();

    ScopeTimer gate_detection(getTimer(WhichTimer::GATE_DETECTION));

    // Clear the vector that contains the information which clauses
    // encode gates and which variables are gate in- and outputs.
    _gates.update();

    // Dependency graph of the gate outputs on the gate inputs.
    // Used to avoid cyclic dependencies, which are not sound in a circuit.
    std::vector<std::vector<Variable>> gateDep(maxVarIndex() + 1);

    // Fill in the dependencies of the gates we have already detected.
    for (const auto& g: _gates) {
        for (const Literal lit: g._input_literals) {
            gateDep[lit2var(g._output_literal)].push_back(lit2var(lit));
        }
    }

    std::size_t count = 0;
    if (and_gates) { count += findANDGates(gateDep, _settings.extend_gates); }
    if (xor_gates) { count += findXORGates(gateDep); }
    if (mux_gates) { count += findMUXGates(gateDep); }
    if (semantic_gates) { count += findSemanticGates(gateDep); }

    // Create a topological ordering of the variables in gateDep.
    _gates.sort();
    VLOG(2) << "We currently have " << _gates.numGates() << " gates (" << count << " new ones).";

    return count > 0;
}


/**
 * \brief Tries to identify equivalent gates (a.k.a. structural hashing).
 *
 * If two gates of the same type have the same inputs, then the gate outputs
 * are equivalent.
 * \return true iff the formula has been changed
 */
bool Formula::findEquivalentGates()
{
    ScopeTimer st(getTimer(WhichTimer::EQUIV_GATES));

    _gates.update();
    _gates.sort();

    const auto& my_gates = _gates.getGates();
    if (my_gates.empty()) return false;

    VLOG(1) << __FUNCTION__ << " called with " << my_gates.size() << " gates.";

    std::size_t count = 0;

    std::map<std::vector<Literal>, Literal> and_map;
    std::map<std::vector<Literal>, Literal> xor_map;
    std::map<std::vector<Literal>, Literal> mux_map;
    std::unordered_map<Literal, Literal> lit_map;

    const auto mapLiteral = [&lit_map](const Literal in) -> Literal
    {
        const auto found = lit_map.find(in);
        if (found == lit_map.cend()) {
            lit_map[in] = in;
            lit_map[negate(in)] = negate(in);
            return in;
        } else {
            return found->second;
        }
    };

    std::vector<Literal> mapped_inputs;

    for (const auto& gate: my_gates) {
        if (!_gates.gateValid(gate)) continue;

        mapped_inputs.clear();

        if (gate._type == GateType::AND_GATE || gate._type == GateType::XOR_GATE || gate._type == GateType::MUX_GATE) {

            std::map<std::vector<Literal>, Literal>* current_map = nullptr;
            switch (gate._type) {
            case GateType::AND_GATE:
                current_map = &and_map;
                break;
            case GateType::XOR_GATE:
                current_map = &xor_map;
                break;
            case GateType::MUX_GATE:
                current_map = &mux_map;
                break;
            default:
                current_map = nullptr;
                break;
            }

            val_assert(current_map);

            for (const Literal lit: gate._input_literals) {
                mapped_inputs.push_back( mapLiteral(lit) );
            }

            // Sort the inputs if the gate has only symmetric inputs.
            // This is the case for AND and XOR, but not for MUX.
            if (gate._type == GateType::AND_GATE || gate._type == GateType::XOR_GATE) {
                std::sort(mapped_inputs.begin(), mapped_inputs.end());
            }

            const auto inputs_found = current_map->find(mapped_inputs);
            if (inputs_found == current_map->cend()) {
                (*current_map)[mapped_inputs] = gate._output_literal;
            } else {
                lit_map[gate._output_literal] = inputs_found->second;
                lit_map[negate(gate._output_literal)] = negate(inputs_found->second);
            }
        }
    }

    for (const auto& entry: lit_map) {
        if (isPositive(entry.first) && entry.first != entry.second) {
            VLOG(3) << "Gate outputs " << lit2dimacs(entry.first) << " and " << lit2dimacs(entry.second) << " are equivalent.";
            replaceLiteral(entry.first, entry.second);

#ifdef SKOLEM
            _skolem_data.push_back(std::make_unique<SkolemEquiv>(lit2var(entry.second), negate_if(entry.first, isNegative(entry.second))));
#endif
            ++count;
        }
    }

    VLOG(2) << __FUNCTION__ << " replaced " << count << " variables.";
    stat(Statistics::EQUIV_GATES) += count;
    return count > 0;
}


/**
 * \brief Only used internally by Formula::findMUXGates.
 */
static bool find_c1(const Clause& c, const Literal out, Literal& sel, Literal& in1)
{
    // Check if c1 is of the form { ~sel, in1, ~out }
    sel = 0;
    in1 = 0;
    for (std::size_t i = 0; i < 3; ++i) {
        if (c[i] == negate(out)) continue;
        else { in1 = sel; sel = c[i]; }
    }
    sel = negate(sel);

    return sel != 0 && in1 != 0;
}

/**
 * \brief Only used internally by Formula::findMUXGates.
 */
static bool find_c3(const Clause& c, const Literal sel, const Literal in1, const Literal out, Literal& in0)
{

    bool found_out = false;
    bool found_sel = false;
    in0 = 0;
    for (std::size_t i = 0; i < 3; ++i) {
        if (c[i] == sel) { found_sel = true; }
        else if (c[i] == negate(out)) { found_out = true; }
        else if (c[i] == in1 || c[i] == negate(in1)) { return false; }
        else { in0 = c[i]; }
    }

    return in0 != 0 && found_sel && found_out;
}


/**
 * \brief Tries to find Tseitin-encoded multiplexer gates in the formula.
 *
 * The implementation roughly follows the following literature:
 * Harald Seltner: <i>Extracting Hardware Circuits from CNF Formulas</i>,
 * Master's Thesis, University of Linz, Austria, July 2014.
 *
 * \return the number of multiplexers found.
 */
std::size_t Formula::findMUXGates(std::vector<std::vector<Variable>>& gateDep)
{
    ScopeTimer sc(getTimer(WhichTimer::GATE_MUX_DETECTION));

    val_assert(_unit_stack.empty());

    Clause::ClauseData tmp_clause(3, 0);
    std::size_t count = 0;

    for (Literal output = minLitIndex(); output <= maxLitIndex(); ++output)
    {
        const Variable output_var = lit2var(output);
        if (!isExistential(output_var)) continue;
        if (_gates.isGateOutput(output_var)) continue;

        bool gate_found = false;
        Literal sel = 0;
        Literal in0 = 0;
        Literal in1 = 0;

        const Literal neg_output = negate(output);
        for (const ClauseID c_nr1: _occ_list[neg_output]) {
            val_assert(!clauseDeleted(c_nr1));

            if (_gates.isGateClause(c_nr1)) continue;
            const Clause& c1 = _clauses[c_nr1];

            // now, sel and i1 are candidates for the select signal and the 1-input.
            // If that does not work, we have to swap them and try again.
            if (c1.size() != 3) continue;
            if (!checkGatePrecondition(lit2var(output), c1)) continue;
            if (!find_c1(c1, output, sel, in1)) continue;

            auto rest_check = [this, &gateDep, &tmp_clause, c_nr1](Literal o, Literal select, Literal i1, Literal& i0) -> bool
                {
                    // check if there is the clause c2 = {~sel, ~in1, out}
                    tmp_clause[0] = negate(select);
                    tmp_clause[1] = negate(i1);
                    tmp_clause[2] = o;
                    std::sort(tmp_clause.begin(), tmp_clause.end());
                    const int c_nr2 = findClause(tmp_clause);
                    val_assert(c_nr2 < 0 || !clauseDeleted(c_nr2));

                    if (c_nr2 >= 0 && !_gates.isGateClause(c_nr2)) {
                        for (const ClauseID c_nr3: _occ_list[select]) {
                            val_assert(!clauseDeleted(c_nr3));

                            if (_gates.isGateClause(c_nr3)) continue;
                            const Clause& c3 = _clauses[c_nr3];

                            // Check if c3 is of the form {sel, in0, ~out}
                            if (c3.size() != 3) return false;
                            if (!checkGatePrecondition(lit2var(o), c3)) return false;
                            if (!find_c3(c3, select, i1, o, i0)) continue;

                            // check if there is the clause c4 = {sel, ~in0, out}
                            tmp_clause[0] = select;
                            tmp_clause[1] = negate(i0);
                            tmp_clause[2] = o;
                            std::sort(tmp_clause.begin(), tmp_clause.end());
                            const int c_nr4 = findClause(tmp_clause);

                            val_assert(c_nr4 < 0 || !clauseDeleted(c_nr4));

                            if (c_nr4 >= 0 && !_gates.isGateClause(c_nr4)) {
                                // check for cyclic dependencies
                                tmp_clause[0] = lit2var(select);
                                tmp_clause[1] = lit2var(i1);
                                tmp_clause[2] = lit2var(i0);
                                if (internal::cyclicDep(gateDep, lit2var(o), tmp_clause)) continue;

                                internal::addDep(gateDep, lit2var(o), tmp_clause);
                                Gate g(GateType::MUX_GATE,  // Type
                                        o,                  // output
                                        { select, i1, i0 }, // inputs
                                        { c_nr1, static_cast<ClauseID>(c_nr2), c_nr3, static_cast<ClauseID>(c_nr4) } // Clauses
                                      );

                                // Check if there are the optional clauses:
                                // (a) { ~in1, ~in0, out }
                                tmp_clause[0] = negate(i1);
                                tmp_clause[1] = negate(i0);
                                tmp_clause[2] = o;
                                std::sort(tmp_clause.begin(), tmp_clause.end());
                                const int c_nr5 = findClause(tmp_clause);
                                if (c_nr5 >= 0) {
                                    val_assert(!clauseDeleted(c_nr5));
                                    g._encoding_clauses.push_back(c_nr5);
                                }

                                // (b) { in1, in0, ~out }
                                tmp_clause[0] = i1;
                                tmp_clause[1] = i0;
                                tmp_clause[2] = negate(o);
                                std::sort(tmp_clause.begin(), tmp_clause.end());
                                const int c_nr6 = findClause(tmp_clause);
                                if (c_nr6 >= 0) {
                                    val_assert(!clauseDeleted(c_nr6));
                                    g._encoding_clauses.push_back(c_nr6);
                                }

                                _gates.addGate(std::move(g));
                                for (auto c_nr: g._encoding_clauses) {
                                    _clauses[c_nr].setStatus(ClauseStatus::MANDATORY);
                                }

                                return true;
                            }
                        }
                    }
                    return false;
            }; // end lambda function

            gate_found = rest_check(output, sel, in1, in0);

            if (!gate_found) {
                std::swap(sel, in1);
                in1 = negate(in1);
                sel = negate(sel);
                gate_found = rest_check(output, sel, in1, in0);
            }

            if (gate_found) {
                ++count;
                break;
            }
        }
    }

    VLOG(2) << __FUNCTION__ << " found " << count << " MUX gates.";
    stat(Statistics::GATES_MUX) += count;
    return count;
}


/**
 * \brief Tries to find Tseitin-encoded multi-input AND gates in the formula.
 *
 * \return the number of AND gates found.
 */
std::size_t Formula::findANDGates(std::vector<std::vector<Variable>>& gateDep, const bool extend)
{
    using namespace internal;

    ScopeTimer sc(getTimer(WhichTimer::GATE_AND_DETECTION));

    val_assert(_unit_stack.empty());

    std::size_t count = 0;

    std::vector<Variable> parents;
    parents.reserve(6);

    // Used for AND-gates to store the clauses which encode the implications.
    std::vector<ClauseID> clause_ids;
    clause_ids.reserve(6);


    // Looking for (multi-)AND-gates with arbitrarily negated inputs and output
    std::vector<Clause::ClauseData> to_add;
    Clause::ClauseData tmp_bin(2,0);

    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (clauseDeleted(c_nr) || _clauses[c_nr].size() <= 2 || _gates.isGateClause(c_nr)) continue;

        for (std::size_t l_nr = 0; l_nr < _clauses[c_nr].size(); ++l_nr)
        {
            const Literal base_lit = _clauses[c_nr][l_nr];
            const Variable base_var = lit2var(base_lit);
            if (_gates.isGateOutput(base_var)) continue;
            if (!checkGatePrecondition(base_var, _clauses[c_nr])) continue;

            parents.clear();
            clause_ids.clear();
            to_add.clear();

            bool has_imp = true;
            for (const Literal lit: _clauses[c_nr]) {
                if (lit == base_lit) continue;
                const int bin_id = getImplicationClause(base_lit, negate(lit));

                if (bin_id < 0) {
                    if (!extend) { has_imp = false; break; }

                    // Try to get the missing implication transitively
                    // or as blocked implication. If both fails,
                    // base_lit is not the output of an AND gate.
                    tmp_bin[0] = negate(base_lit);
                    tmp_bin[1] = lit;

                    if (clauseBlocked(tmp_bin)) {
                        to_add.push_back(tmp_bin);
                        continue;
                    }

                    const auto trans = hasImplicationTransitive(base_lit, lit2var(lit));
                    if ( (isPositive(lit) && trans.first) || (isNegative(lit) && trans.second) ) {
                        to_add.push_back(tmp_bin);
                        continue;
                    }

                    has_imp = false;
                    break;
                } else if (_gates.isGateClause(bin_id)) {
                    // The implication is available, but already used for another gate.
                    has_imp = false;
                    break;
                }
                parents.push_back(lit2var(lit));
                clause_ids.push_back(bin_id);
            } // end for

            if (!has_imp) continue;

            if (!to_add.empty()) {
                VLOG(3) << "If we add " << to_add.size() << " clause(s), we get an additional AND gate.";

                for (auto& cl: to_add) {
                    const auto impl_nr = addClause(std::move(cl));
                    if (impl_nr >= 0) clause_ids.push_back(impl_nr);
                }
            }

            // Check for cyclic dependencies
            if (cyclicDep(gateDep, base_var, parents)) continue;

            addDep(gateDep, base_var, parents);

            // We have found an AND-gate
            Gate g(GateType::AND_GATE, base_lit);
            g._encoding_clauses.reserve(_clauses[c_nr].size());
            g._input_literals.reserve(_clauses[c_nr].size() - 1);
            g._encoding_clauses.push_back(c_nr);
            _clauses[c_nr].setStatus(ClauseStatus::MANDATORY);
            for (const ClauseID bin_id: clause_ids) {
                g._encoding_clauses.push_back(bin_id);
                _clauses[bin_id].setStatus(ClauseStatus::MANDATORY);
            }
            for (const Literal lit: _clauses[c_nr]) {
                if (lit != base_lit) g._input_literals.push_back(negate(lit));
            }

            _gates.addGate(std::move(g));
            ++count;
            break;
        }
    }

    VLOG(2) << __FUNCTION__ << " found " << count << " AND gates.";
    stat(Statistics::GATES_AND) += count;
    return count;
}

/**
 * \brief Tries to find Tseitin-encoded 2-input XOR gates in the formula.
 *
 * \return the number of XOR gates found.
 */
std::size_t Formula::findXORGates(std::vector<std::vector<Variable>>& gateDep)
{
    using namespace internal;

    ScopeTimer sc(getTimer(WhichTimer::GATE_XOR_DETECTION));

    val_assert(_unit_stack.empty());

    std::size_t count = 0;

    std::vector<Variable> parents;
    parents.reserve(6);

    // Looking for XOR-gates with arbitrarily negated inputs and output
    for (Variable base_var = 1; base_var <= maxVarIndex(); ++base_var)
    {
        if (isUniversal(base_var)) continue;
        if (_gates.isGateOutput(base_var)) continue;
        if (_occ_list[var2lit(base_var, false)].size() < 2 || _occ_list[var2lit(base_var, true)].size() < 2) continue;

        const Literal base_lit_pos = var2lit(base_var, false);
        const Literal base_lit_neg = var2lit(base_var, true);
        bool found_gate = false;

        for (std::size_t p1 = 0; p1 < _occ_list[base_lit_pos].size(); ++p1) {
            const ClauseID i1 = _occ_list[base_lit_pos][p1];
            val_assert(!clauseDeleted(i1));
            const auto& c1 = _clauses[i1];
            if (c1.size() != 3) continue;
            if (_gates.isGateClause(i1)) continue;
            if (!checkGatePrecondition(base_var, c1)) continue;

            // Check for cyclic dependencies.
            parents.clear();
            for (const Literal lit: c1) {
                if (lit == base_lit_pos) continue;
                parents.push_back(lit2var(lit));
            }
            if (cyclicDep(gateDep, base_var, parents)) continue;

            for (std::size_t p2 = p1 + 1; p2 < _occ_list[base_lit_pos].size(); ++p2) {
                const ClauseID i2 = _occ_list[base_lit_pos][p2];
                val_assert(!clauseDeleted(i2));
                const auto& c2 = _clauses[i2];
                if (_gates.isGateClause(i2)) continue;
                if (!sameVars(c1, c2)) continue;

                for (std::size_t p3 = 0; p3 < _occ_list[base_lit_neg].size(); ++p3) {
                    const ClauseID i3 = _occ_list[base_lit_neg][p3];
                    val_assert(!clauseDeleted(i3));
                    const auto& c3 = _clauses[i3];
                    if (_gates.isGateClause(i3)) continue;
                    if (!sameVars(c1, c3)) continue;

                    for (std::size_t p4 = p3 + 1;  p4 < _occ_list[base_lit_neg].size(); ++p4) {
                        const ClauseID i4 = _occ_list[base_lit_neg][p4];
                        val_assert(!clauseDeleted(i4));
                        const auto& c4 = _clauses[i4];
                        if (_gates.isGateClause(i4)) continue;
                        if (!sameVars(c1, c4)) continue;

                        if (c1 == c2 || c1 == c3 || c1 == c4 || c2 == c3 || c2 == c4 || c3 == c4) continue;

                        const auto pc1 = posCount(c1);
                        const auto pc2 = posCount(c2);
                        const auto pc3 = posCount(c3);
                        const auto pc4 = posCount(c4);

                        if( !( (  ( pc1 & 1 ) &&  ( pc2 & 1 ) &&  ( pc3 & 1 ) &&  ( pc4 & 1 ) ) ||
                               ( !( pc1 & 1 ) && !( pc2 & 1 ) && !( pc3 & 1 ) && !( pc4 & 1 ) ) ) )
                        {
                            // mixed polarities -> skipping
                            continue;
                        }

                        // Now 'base_var' is the output of an XOR gate.
                        addDep(gateDep, base_var, parents);
                        Gate g(GateType::XOR_GATE,
                                var2lit(base_var, (pc1 & 1) == 1),
                                {var2lit(parents[0]), var2lit(parents[1])},
                                { i1, i2, i3, i4 }
                              );
                        _clauses[i1].setStatus(ClauseStatus::MANDATORY);
                        _clauses[i2].setStatus(ClauseStatus::MANDATORY);
                        _clauses[i3].setStatus(ClauseStatus::MANDATORY);
                        _clauses[i4].setStatus(ClauseStatus::MANDATORY);
                        _gates.addGate(std::move(g));
                        ++count;
                        found_gate = true;
                        break;
                    } // for c4
                    if (found_gate) break;
                } // for c3
                if (found_gate) break;
            } // for c2
            if (found_gate) break;
        } // for c1
    } // for base_var
    // Detection of XOR gates finished.

    VLOG(2) << __FUNCTION__ << " found " << count << " XOR gates.";
    stat(Statistics::GATES_XOR) += count;
    return count;
}


/**
 * \brief Implements semantic gate detection.
 *
 * A semantic gate definition is a set of clauses which uniquely determines the
 * value of a output variable.
 *
 * For more information see the following two papers:
 * - Rabe, Seshia: <i>Incremental Determinization</i>,
 *   Proc. of SAT 2016, Lecture Notes in Computer Science, vol. 9710,
 *   pages 375--392, Springer. DOI: 10.1007/978-3-319-40970-2_23
 * - Balabanov, Jiang, Mishchenko, and Scholl:
 *   <i>Clauses Versus Gates in CEGAR-Based 2QBF Solving</i>,
 *   AAAI Workshop on Beyond NP, 2016,
 *   URL: http://www.aaai.org/ocs/index.php/WS/AAAIW16/paper/view/12660
 * \return the number of found gates.
 */
std::size_t Formula::findSemanticGates(std::vector<std::vector<Variable>>& gateDep)
{
    using namespace internal;

    val_assert(_unit_stack.empty());

    ScopeTimer sc(getTimer(WhichTimer::GATE_SEMANTIC_DETECTION));

    // If semantic checks are requested, use a SAT solver for the variables
    // which have not yet been identified as a gate output.
    std::vector<Variable> parents;

    std::size_t count = 0;
    Timer sem_time;
    sem_time.start();
    std::vector<ClauseID> used_clauses;

    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (sem_time.read() > 5.0) break;
        if (!isExistential(var) || _gates.isGateOutput(var)) continue;

        PicoSAT* sat = picosat_init();
        picosat_enable_trace_generation(sat);

        // Add all clauses containing var or ~var, in which var is the
        // right-most variable, but remove all occurrences of var before.
        const Literal out_pos = var2lit(var);
        const Literal out_neg = negate(out_pos);

        used_clauses.clear();
        used_clauses.reserve( _occ_list[out_pos].size() + _occ_list[out_neg].size() );

        for (const ClauseID c_nr: _occ_list[out_pos]) {
            val_assert(!clauseDeleted(c_nr));

            if (_gates.isGateClause(c_nr)) continue;
            bool add = true;
            for (Literal lit: _clauses[c_nr]) {
                if (!dependenciesSubset(lit2var(lit), var)) {
                    add = false;
                    break;
                }
            }
            if (add) {
                for (Literal lit: _clauses[c_nr]) {
                    if (lit2var(lit) != var) picosat_add(sat, lit2dimacs(lit));
                }
                picosat_add(sat, 0);
                used_clauses.push_back(c_nr);
            }
        }

        for (const ClauseID c_nr: _occ_list[out_neg]) {
            val_assert(!clauseDeleted(c_nr));

            if (_gates.isGateClause(c_nr)) continue;
            bool add = true;
            for (Literal lit: _clauses[c_nr]) {
                if (!dependenciesSubset(lit2var(lit), var)) {
                    add = false;
                    break;
                }
            }
            if (add) {
                for (Literal lit: _clauses[c_nr]) {
                    if (lit2var(lit) != var) picosat_add(sat, lit2dimacs(lit));
                }
                picosat_add(sat, 0);
                used_clauses.push_back(c_nr);
            }
        }

        // Check if the formula is satisfiable.
        const int ret = picosat_sat(sat, -1);

        if (ret != PICOSAT_UNSATISFIABLE) {
            picosat_reset(sat);
            continue;
        }

        parents.clear();
        for (Variable other_var = minVarIndex(); other_var <= maxVarIndex(); ++other_var) {
            if (varDeleted(other_var)) continue;
            if (other_var != var && picosat_corelit(sat, other_var)) parents.push_back(other_var);
        }

        if (_gates.isGateInput(var) > 0 && cyclicDep(gateDep, var, parents)) {
            picosat_reset(sat);
            continue;
        }

        Gate g(GateType::UNKNOWN, var2lit(var, false));
        g._input_literals.reserve(parents.size());
        for (Variable in: parents) {
            g._input_literals.push_back(var2lit(in));
        }
        for (std::size_t cl = 0; cl < used_clauses.size(); ++cl) {
            if (picosat_coreclause(sat, static_cast<int>(cl))) {
                g._encoding_clauses.push_back(used_clauses[cl]);
                _clauses[used_clauses[cl]].setStatus(ClauseStatus::MANDATORY);
            }
        }
        picosat_reset(sat);

        // Check if the gate is consistent
        Antom solver;
        solver.setMaxIndex(maxVarIndex() + static_cast<Variable>(g._encoding_clauses.size()));
        std::vector<Literal> pos;
        std::vector<Literal> neg;
        unsigned int current_index = maxVarIndex() + 1;

        std::vector<Literal> long_clause;
        std::vector<Literal> bin_clause(2,0);
        for (const ClauseID c_nr: g._encoding_clauses) {
            long_clause.clear();
            bin_clause[0] = var2lit(current_index, true);
            long_clause.push_back(var2lit(current_index, false));

            for (const Literal lit: _clauses[c_nr]) {
                if (lit == g._output_literal) {
                    pos.push_back(var2lit(current_index));
                } else if (lit == negate(g._output_literal)) {
                    neg.push_back(var2lit(current_index));
                } else {
                    bin_clause[1] = negate(lit);
                    long_clause.push_back(lit);
                    solver.addClause(bin_clause);
                }
            }
            solver.addClause(long_clause);
            ++current_index;
        }


        if (pos.empty()) {
            pushUnit(negate(g._output_literal), PureStatus::UNIT);
            unitPropagation();
        } else if (neg.empty()) {
            pushUnit(g._output_literal, PureStatus::UNIT);
            unitPropagation();
        } else {
            solver.addClause(pos);
            solver.addClause(neg);
            const auto result = solver.solve();

            if (result == TruthValue::FALSE) {
                addDep(gateDep, var, parents);
                _gates.addGate(std::move(g));
                ++count;
            }
        }
    }

    VLOG(2) << __FUNCTION__ << " found " << count << " semantic gates.";
    stat(Statistics::GATES_SEMANTIC) += count;

    return count;
}

bool Formula::isGateOutput(const Variable var) const noexcept
{
    return _gates.isGateOutput(var);
}

} // end namespace hqspre
