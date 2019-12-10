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

#include <cstdint>
#include <random>
#include <vector>

#include <easylogging++.hpp>

#include "aux.hpp"
#include "clause.hpp"
#include "formula.hpp"
#include "literal.hpp"
#include "prefix.hpp"
#include "sat_solver.hpp"
#include "timer.hpp"

/**
 * \file formula_sat.cpp
 * \brief Implementations of functions which try to identify unit
 *        literals, implications, and equivalences using SAT-solver calls.
 * \author Ralf Wimmer
 * \author Sven Reimer
 * \date 03-06/2016
 */

namespace hqspre {


#ifndef SatSolver
#define SatSolver Antom
#endif
#define macro_print(TXT) str(TXT)
#define str(s) #s


/**
 * \brief Tries to find backbones by SAT calls.
 *
 * Backbones are variables which have the same value in
 * all satisfying assignments of the matrix. They can be handled
 * like unit literals and replaced by a fixed value.
 * \return true if the formula was modified.
 */
bool Formula::findConstantsBySAT()
{
    VLOG(1) << __FUNCTION__;

    ScopeTimer unit_sat(getTimer(WhichTimer::UNIT_SAT));

    _process_limit.setLimit(PreproMethod::UNIT_BY_SAT);

    std::size_t count = 0;
    std::size_t num_sat_calls = 0;
    std::size_t number_of_literals = 0;
    std::size_t number_of_clauses = 0;

    val_assert(_seen.size() > maxVarIndex());
    std::vector<Variable> candidates(maxVarIndex() + 1);

    // Create a SAT solver and feed it with all clauses.
    Antom sat_solver;
    sat_solver.setMaxIndex(maxVarIndex());
    sat_solver.setTimeout(1.0);
    std::vector<Literal> antom_clause;

    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (clauseDeleted(c_nr)) continue;
        if (clauseOptional(c_nr)) continue;
        antom_clause.assign(_clauses[c_nr]._literals.cbegin(), _clauses[c_nr]._literals.cend());

        sat_solver.addClause(antom_clause);
        number_of_literals += _clauses[c_nr].size();
        number_of_clauses++;
    }

    if (number_of_clauses == 0) return false;
    std::size_t limit = number_of_literals * (number_of_literals / number_of_clauses);

    std::size_t number_of_candidates = 0;
    for (Variable v = minVarIndex(); v <= maxVarIndex(); ++v) {
        if (varAssignment(v) != 0 || varDeleted(v)) continue;
        candidates[number_of_candidates] = v;
        ++number_of_candidates;
    }

    candidates.resize(number_of_candidates);

    if (_interrupt) return false;

    VLOG(3) << __FUNCTION__ << " checking " << number_of_candidates << " candidates.";

    // Solve the formula without any assumptions and record
    // the assignment of the variables.
    TruthValue sat_result = sat_solver.solve();
    ++num_sat_calls;

    if (_interrupt) return false;
    _process_limit.decreaseLimitBy(limit);

    // run into timeout
    if (sat_result == TruthValue::UNKNOWN) {
        VLOG(2) << __FUNCTION__ << " initial SAT-instance runs into timeout -> abort.";
        clearSeen();
        return false;
    } else if (sat_result == TruthValue::FALSE) {
        throw UNSATException("The matrix is already unsatisfiable");
    }

    for (std::size_t i = 0; i != candidates.size(); ++i) {
        const TruthValue value = sat_solver.getValue(candidates[i]);
        val_assert( value != TruthValue::UNKNOWN );
        _seen[candidates[i]] = (value == TruthValue::TRUE);
        // force solver to take opposite polarity
        // strategy = 3 => set always true, strategy = 2 => set always false
        sat_solver.setPolarity(candidates[i], value == TruthValue::FALSE);
    }

    // For the assignments that haven't been seen so far,
    // call the SAT-solver incrementally with the assumption
    // forcing the current variable to the unseen assignment.
    // If formula becomes UNSAT -> variable is unit.
    // Record the satisfying assignment for all remaining variables.
    std::vector<Literal> assumptions(1,0);

    std::size_t size = candidates.size();
    for (std::size_t i = 0; i < size; ++i) {

        if (_process_limit.reachedLimit()) {
            VLOG(2) << "Terminate " << __FUNCTION__ << " due to process limit.";
            break;
        }
        if (_interrupt) break;

        const Variable nextvar = candidates[i];
        const bool polarity = _seen[nextvar];
        assumptions[0] = var2lit(nextvar,polarity);
        sat_result = sat_solver.solve(assumptions);
        ++num_sat_calls;

        _process_limit.decreaseLimitBy(limit);

        if (sat_result == TruthValue::TRUE) {
            // Find candidates which cannot be
            for (std::size_t j = i; j < size; ++j) {
                const Variable check_var = candidates[j];
                const TruthValue check_value = sat_solver.getValue(check_var);
                val_assert( check_value != TruthValue::UNKNOWN );
                // We have seen both polarities of "check_var" -> we can remove this candidate from the list
                if (_seen[check_var] == (check_value == TruthValue::FALSE) ) {
                    candidates[j--] = candidates[--size];
                    candidates.pop_back();
                }
            }
        } else if (sat_result == TruthValue::FALSE) {
            // no satisfying assignment where lit = true.
            // -> !lit is unit.
            pushUnit(var2lit(nextvar,!polarity), PureStatus::UNIT);
            VLOG(3) << __FUNCTION__ << "() " << nextvar << " is unit.";

            // Also add this information to the SAT solver
            sat_solver.addUnit(var2lit(nextvar,!polarity));
            ++count;
        }
    }

    VLOG(2) << __FUNCTION__ << " found " << count << " unit literals in " << num_sat_calls << " SAT solver calls.";
    VLOG_IF(_interrupt || _process_limit.reachedLimit(), 2) << __FUNCTION__ << ' ' << size << " candidates left";

    stat(Statistics::UNIT_SAT) += count;
    stat(Statistics::SAT_CALLS) += num_sat_calls;
    clearSeen();
    return count > 0;
}



/**
 * \brief Identifies all binary clauses that are consequences of the formula.
 *
 * For all binary clauses which are not derivable from implications by transitivity,
 * a SAT-solver is called to check if they are implied by the formula. If this is
 * the case, the binary clause is added.
 * \return true if the formula was modified.
 */
bool Formula::findImplicationsBySAT()
{
    VLOG(1) << __FUNCTION__;
    ScopeTimer impl_sat(getTimer(WhichTimer::IMPLICATIONS_SAT));

    _process_limit.setLimit(PreproMethod::IMPL_BY_SAT);

    std::size_t subsumptions = stat(Statistics::SUBSUMPTION);
    std::size_t count = 0;
    std::size_t unit_count = 0;
    std::size_t num_sat_calls = 0;

    const std::uint64_t max_entry = (static_cast<std::uint64_t>(maxLitIndex()) + 2ul) * static_cast<std::uint64_t>(maxLitIndex());

    // Exit if cache data structure would be too large
    if (max_entry > (1ul<<31)) {
        return false;
    }

    auto lit2index = [this](Literal lit1, Literal lit2) -> std::uint64_t
        {
            val_assert(lit1 < lit2);
            return (static_cast<std::uint64_t>(lit1) * static_cast<std::uint64_t>(maxLitIndex() + 1) + static_cast<std::uint64_t>(lit2));
        };

    std::vector<bool> cache(max_entry + 1, false);

    std::size_t number_of_literals = 0;
    std::size_t number_of_clauses = 0;

    // Create a SAT solver and feed it with all clauses
    Antom sat_solver;
    sat_solver.setMaxIndex(maxVarIndex());
    sat_solver.setTimeout(1.0);
    std::vector<Literal> antom_clause;

    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (clauseDeleted(c_nr)) continue;
        antom_clause.assign(_clauses[c_nr]._literals.cbegin(), _clauses[c_nr]._literals.cend());
        sat_solver.addClause(antom_clause);
        number_of_literals += _clauses[c_nr].size();
        number_of_clauses++;
    }

    if (number_of_clauses == 0 ) {
        throw SATException("No clauses left");
    }
    std::size_t limit = number_of_literals * (number_of_literals/number_of_clauses);

    // Create a first satisfying assignment of the matrix
    // and store the value of all pairs of variables in the
    // cache.
    const TruthValue result = sat_solver.solve();

    ++num_sat_calls;
    if (result == TruthValue::UNKNOWN ) {
        VLOG(2) << __FUNCTION__ << " initial SAT-instance runs into timeout -> abort";
        return false;
    } else if (result == TruthValue::FALSE) {
        throw UNSATException("Matrix is unsatisfiable.");
    }

    for (Variable var1 = minVarIndex(); var1 <= maxVarIndex(); ++var1)
    {
        if (varDeleted(var1) || varAssignment(var1) != 0) continue;
        const TruthValue value1 = sat_solver.getValue(var1);
        const Literal lit1 = var2lit(var1, value1 == TruthValue::FALSE);
        sat_solver.setPolarity(var1, value1 == TruthValue::FALSE);
        _process_limit.decreaseLimitBy(maxVarIndex());
        for (Variable var2 = var1 + 1; var2 <= maxVarIndex(); ++var2)
        {
            const TruthValue value2 = sat_solver.getValue(var2);
            const Literal lit2 = var2lit(var2, value2 == TruthValue::FALSE);
            if (varDeleted(var2) || varAssignment(var1) != 0) continue;
            if (value1 != TruthValue::UNKNOWN && value2 != TruthValue::UNKNOWN) {
                val_assert(lit2index(lit1, lit2) <= max_entry);
                cache[lit2index(lit1, lit2)] = true;
            }
        }
    }

    // For the missing assignments, enforce them by adding assumptions
    // and check if they exist. If not, we can add an implication.
    std::vector<Literal> assumptions(2, 0);
    for (Literal lit1 = minLitIndex(); lit1 < maxLitIndex(); ++lit1) {

        const Variable var1 = lit2var(lit1);
        if (varDeleted(var1) || varAssignment(var1) != 0) continue;

        for (Literal lit2 = lit1 + 1; lit2 <= maxLitIndex(); ++lit2) {

            if (_process_limit.reachedLimit()) {
                VLOG(2) << "Terminate " << __FUNCTION__ << " due to process limit.";
                break;
            }

            const Variable var2 = lit2var(lit2);
            if (varDeleted(var2) || var1 == var2 || varAssignment(var2) != 0) continue;
            if (varAssignment(var1) != 0) break;
            if (cache[lit2index(lit1, lit2)]) continue;
            _process_limit.decrementLimit();
            const auto has_impl = hasImplicationTransitive(lit1, lit2var(lit2));
            // ~lit1 implies lit2 and ~lit2 => lit1 is unit
            if (has_impl.first && has_impl.second) {
                VLOG(3) << __FUNCTION__ << ": " << lit2dimacs(lit1) << " implies " << lit2dimacs(lit2) << " and " << lit2dimacs(negate(lit2)) << " => add unit " << lit2dimacs(negate(lit1));
                pushUnit(negate(lit1), PureStatus::UNIT);
                if (!sat_solver.addUnit(negate(lit1)) ) {
                    throw UNSATException("Matrix is unsatisfiable.");
                }
                ++unit_count;
                continue;
            }

            if ((isNegative(lit2) && has_impl.first) || (isPositive(lit2) && has_impl.second)) continue;

            assumptions[0] = lit1;
            assumptions[1] = lit2;

            const TruthValue sat_result = sat_solver.solve(assumptions);
            ++num_sat_calls;
            _process_limit.decreaseLimitBy(limit);
            if (sat_result == TruthValue::FALSE) {
                ++count;
                val_assert(negate(lit1) < negate(lit2));
                assumptions[0] = negate(lit1);
                assumptions[1] = negate(lit2);

                VLOG(3) << __FUNCTION__ << "() add hidden implication " << lit2dimacs(assumptions[0]) << " " << lit2dimacs(assumptions[1]);

                // The function addClause marks the clause automatically
                // as mandatory if it subsumes an existing clause.
                addClause(assumptions, false, false, ClauseStatus::OPTIONAL);
            } else if (sat_result == TruthValue::TRUE) {
                // update the cache
                for (Variable curr_var1 = var1; curr_var1 <= maxVarIndex(); ++curr_var1) {
                    if (varDeleted(curr_var1) || varAssignment(curr_var1) != 0) continue;

                    const TruthValue curr_value1 = sat_solver.getValue(curr_var1);
                    const Literal curr_lit1 = var2lit(curr_var1, curr_value1 == TruthValue::FALSE);
                    sat_solver.setPolarity(curr_var1, curr_value1 == TruthValue::FALSE);

                    for (Variable curr_var2 = curr_var1 + 1; curr_var2 <= maxVarIndex(); ++curr_var2) {
                        if (varDeleted(curr_var2) || varAssignment(curr_var1) != 0) continue;
                        const TruthValue curr_value2 = sat_solver.getValue(curr_var2);
                        const Literal curr_lit2 = var2lit(curr_var2, curr_value2 == TruthValue::FALSE);
                        _process_limit.decrementLimit();
                        if (curr_value1 != TruthValue::UNKNOWN && curr_value2 != TruthValue::UNKNOWN) {
                            cache[lit2index(curr_lit1, curr_lit2)] = true;
                        }
                    }
                }
            }
        }

        if (_process_limit.reachedLimit()) {
            VLOG(2) << "Terminate " << __FUNCTION__ << " due to process limit";
            break;
        }
    }

    VLOG(2) << __FUNCTION__ << " added " << count << " implications and " << unit_count << " units after "
            << num_sat_calls << " SAT-solver calls and removed " << (stat(Statistics::SUBSUMPTION) - subsumptions)
            << " clauses by subsumption.";

    stat(Statistics::UNIT_SAT) += unit_count;
    stat(Statistics::IMPL_SAT) += count;
    stat(Statistics::SAT_CALLS) += num_sat_calls;

    return count > 0;
}


/**
 * \brief Performs an incomplete check for satisfiability.
 *
 * A SAT-solver is called for the formula obtained by removing all
 * universal literals. If the remaining purely existential formula
 * is satisfiable, so is the original (D)QBF.
 * \throw SATException
 */
void Formula::checkSAT()
{
    VLOG(1) << __FUNCTION__ << " using " << macro_print(SatSolver) << " as SAT solver";

    ScopeTimer sat_check(getTimer(WhichTimer::INCOMPLETE_SAT));

    std::vector<Literal> clause;
    clause.reserve(maxVarIndex());

    // Create a SAT solver and feed it with all clauses.
    SatSolver sat_solver;
    sat_solver.setMaxIndex(maxVarIndex());
    if (numUVars() > 0) {
        sat_solver.setTimeout(_settings.sat_timeout);
    } else {
        VLOG(1) << "Problem is a pure SAT problem!";
        sat_solver.setTimeout(_settings.pure_sat_timeout);
    }

    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (clauseDeleted(c_nr)) continue;
        if (clauseOptional(c_nr)) continue;
        for (Literal lit: _clauses[c_nr]) {
            if (isExistential(lit2var(lit))) clause.push_back(lit);
        }
        sat_solver.addClause(clause);
        clause.clear();
    }

    const TruthValue result = sat_solver.solve();
    ++stat(Statistics::SAT_CALLS);

    if (numUVars() == 0) {
        if (result == TruthValue::TRUE) {
            throw SATException("Formula is a SAT problem and satisfiable");
        } else if (result == TruthValue::FALSE) {
            throw UNSATException("Formula is a SAT problem and unsatisfiable");
        } else {
            VLOG(2) << __FUNCTION__ << ": SAT solver returned UNKNOWN.";
            _settings.sat_incomplete = false;
        }
        return;
    }

    if (result == TruthValue::TRUE) {
        throw SATException("Formula has constant Skolem functions.");
    }
}

/**
 * \brief Performs an incomplete check for unsatisfiability
 *
 * We take an assignment of the universal variables which satisfies the
 * least number of clauses and solve the resulting SAT instance. If it
 * is unsatisfiable, so is the (D)QBF. Otherwise we flip the value of
 * the universal variables and check again.
 */
void Formula::checkUNSAT(const std::size_t num_random_patterns)
{
    VLOG(1) << __FUNCTION__ << " using " << macro_print(SatSolver) << " as SAT solver.";

    ScopeTimer sat_check(getTimer(WhichTimer::INCOMPLETE_SAT));

    // Create a SAT solver and feed it with all clauses.
    SatSolver sat_solver;
    sat_solver.setMaxIndex(maxVarIndex());
    sat_solver.setTimeout(_settings.sat_timeout);
    std::vector<Literal> antom_clause;

    for (ClauseID c_nr = 0; c_nr < _clauses.size(); ++c_nr) {
        if (clauseDeleted(c_nr)) continue;
        if (clauseOptional(c_nr)) continue;
        antom_clause.assign(_clauses[c_nr]._literals.cbegin(), _clauses[c_nr]._literals.cend());
        sat_solver.addClause(antom_clause);
    }

    std::vector<Literal> assumptions;
    assumptions.reserve(_prefix->numUVars());

    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (!isUniversal(var)) continue;
        const Literal lit = var2lit(var, false);
        assumptions.push_back( _occ_list[lit].size() < _occ_list[negate(lit)].size() ? lit : negate(lit) );
    }

    TruthValue result = sat_solver.solve(assumptions);
    ++stat(Statistics::SAT_CALLS);

    if (numUVars() == 0) {
        if (result == TruthValue::TRUE) {
            throw SATException("Formula is a SAT problem and satisfiable.");
        } else if (result == TruthValue::FALSE) {
            throw UNSATException("Formula is a SAT problem and unsatisfiable.");
        }
        return;
    }

    if (result != TruthValue::FALSE) {
        for (Literal& lit: assumptions) lit = negate(lit);
        result = sat_solver.solve(assumptions);
        ++stat(Statistics::SAT_CALLS);
    }

    if (result == TruthValue::FALSE) {
        throw UNSATException("Matrix unsatisfiable for assignment of universals.");
    }

    if (num_random_patterns == 0) return;

    // Simulate a number of random patterns for the universal variables.
    std::mt19937 gen(394871); // Take an arbitrary, but fixed seed for reproducability
    std::bernoulli_distribution d(0.5);

    VLOG(2) << __FUNCTION__ << " checking " << num_random_patterns << " random patterns for the universal variables.";

    for (std::size_t current_random_pattern = 0; current_random_pattern < num_random_patterns; ++current_random_pattern) {
        for (std::size_t pos = 0; pos < assumptions.size(); ++pos) {
            const Variable current_var = lit2var(assumptions[pos]);

            assumptions[pos] = var2lit(current_var, d(gen));
        }
        result = sat_solver.solve(assumptions);
        ++stat(Statistics::SAT_CALLS);
        if (result == TruthValue::FALSE) {
            throw UNSATException("Matrix unsatisfiable for random assignment of universals.");
        } else if (result == TruthValue::UNKNOWN) {
            break;
        }
    }

}


} // end namespace hqspre

