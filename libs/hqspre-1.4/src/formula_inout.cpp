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
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <easylogging++.hpp>
#include "aux.hpp"
#include "formula.hpp"
#include "literal.hpp"
#include "prefix.hpp"

/**
 * \file formula_inout.cpp
 * \brief Implementation of methods and functions for reading and writing formulas.
 * \author Ralf Wimmer
 * \date 01/2016
 */


namespace hqspre {

/**
 * \brief Reads a formula in DQDIMACS format from an input stream
 *
 * If the current formula is not empty, it is cleared first by calling Formula::reset().
 */
void Formula::read(std::istream& stream)
{
    val_assert(_prefix);

    reset();

    std::string token;
    Variable num_vars = 0;
    std::size_t num_clauses = 1;
    std::size_t clauses_read = 0;
    std::set<Variable> universal_vars;
    std::vector<Variable> var_map;
    Clause::ClauseData clause;

    while (!stream.eof()) {
        stream >> token;

        if (token == "p") {
            // header
            std::string s;
            stream >> s; // cnf
            stream >> num_vars;
            stream >> num_clauses;
            var_map.resize(num_vars + 1, 0);
            _clauses.reserve(num_clauses);
            setMaxVarIndex(num_vars);
        } else if (token == "a") {
            // list of universally quantified variables
            while (stream) {
                Variable var = 0;
                stream >> var;
                if (var == 0) break;
                val_assert(minVarIndex() <= var && var <= maxVarIndex());
                var_map[var] = addUVar();
                universal_vars.insert(var_map[var]);
            }
        } else if (token == "e") {
            // existential variable depending on all universal vars 
            // declared so far.
            Variable exist_var = 0;
            while (stream) {
                stream >> exist_var;
                if (exist_var == 0) break;
                val_assert(exist_var <= maxVarIndex());
                if (_prefix->type() == PrefixType::DQBF) {
                    var_map[exist_var] = addEVar(universal_vars);
                } else {
                    var_map[exist_var] = addEVar();
                }
            }
        } else if (token == "d") {
            // If we have a QBF prefix and encounter a "d"-line,
            // we have to switch to DQBF.
            if (_prefix->type() == PrefixType::QBF) {
                VLOG(2) << "Switching to DQBF prefix representation.";

                _dqbf_prefix = _qbf_prefix->convertToDQBF();
                delete _prefix;
                _prefix = _dqbf_prefix;
                _qbf_prefix = nullptr;
            }
            // existential variable with dependencies
            Variable exist_var = 0;
            Variable all_var = 0;

            stream >> exist_var;
            val_assert(exist_var <= maxVarIndex());

            std::set<Variable> deps;
            while (stream) {
                stream >> all_var;
                if (all_var == 0) break;
                deps.insert(var_map[all_var]);
            }
            var_map[exist_var] = addEVar(std::move(deps));
        } else if (token == "c") {
            std::getline(stream, token); // consume the rest of the line
        } else {
            // clause
            clause.clear();
            bool eof = false;
            int lit = std::atoi(token.c_str());
            do {
                if (lit == 0) break;
                const Variable var = var_map[abs(lit)];
                clause.push_back(var2lit(var, lit < 0));
                if (stream.eof()) {
                    eof = true;
                    break;
                }
                stream >> lit;
            } while (stream);
            if (!eof) {
                addClause(std::move(clause));
                ++clauses_read;
            } else break;
            if (clauses_read == num_clauses) break;
        }
    }

    // Check if we can switch from DQBF to QBF
    if (!_enforce_dqbf && _prefix->type() == PrefixType::DQBF && isQBF()) {
        val_assert(!_qbf_prefix && _dqbf_prefix);

        VLOG(2) << "Switching to QBF prefix representation.";
        _qbf_prefix = _dqbf_prefix->convertToQBF();
        delete _dqbf_prefix;
        _dqbf_prefix = nullptr;
        _prefix = _qbf_prefix;
    }

    VLOG(1) << "Read formula with " << numEVars() << " existential vars, " << numUVars() << " univ. vars, and " << numClauses() << " clauses.";
}


/**
 * \brief Writes the current formula in DQDIMACS format to the given output stream.
 *
 * If the formula is a QBF, then the output file is actually in QDIMACS format.
 * \param stream output stream to which the formula should be written
 * \param compact rename the variables such that there are no deleted variables in between
 */
void Formula::write(std::ostream& stream, bool compact) const
{
    val_assert(_prefix);

    if (!stream) {
        LOG(ERROR) << "Could not write to output stream!\n";
        std::exit(-1);
    }

    const_cast<Formula*>(this)->updateVars();
    val_assert(checkConsistency());

    std::vector<Variable> translation_table(maxVarIndex() + 1, 0);
    Variable current = 0;
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (varDeleted(var)) continue;
        if (compact) { translation_table[var] = ++current; }
        else translation_table[var] = var;
    }

    // If necessary, print the translation table
    if (compact) {
        for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
            if (!varDeleted(var)) {
                stream << "c ren " << var << " -> " << translation_table[var] << '\n';
            }
        }
    }
    // print the header
    stream << "p cnf " << (compact ? current : maxVarIndex()) << ' ' << numClauses() << '\n';
    _prefix->write(stream, &translation_table);
    writeClauses(stream, &translation_table);
}

/**
 * \brief Print all clauses in DIMACS format to the given stream.
 */
void Formula::writeClauses(std::ostream& stream, std::vector<Variable>* translation_table) const
{
    auto trans = [translation_table](const Literal lit) -> Literal {
        if (!translation_table) return lit;
        else return var2lit((*translation_table)[lit2var(lit)], isNegative(lit));
    };

    // Print the clauses.
    for (ClauseID c_nr = 0; c_nr <= maxClauseIndex(); ++c_nr) {
        if (clauseDeleted(c_nr)) continue;
        for (const Literal lit: _clauses[c_nr]) {
            stream << lit2dimacs(trans(lit)) << ' ';
        }
        stream << "0\n";
    }

    for (const Literal lit: _unit_stack) {
        stream << lit2dimacs(lit) << " 0\n";
    }
}


std::ostream& operator<<(std::ostream& stream, const Formula& formula)
{
    formula.write(stream);
    return stream;
}


std::istream& operator>>(std::istream& stream, Formula& formula)
{
    formula.read(stream);
    return stream;
}


} // end namespace hqspre
