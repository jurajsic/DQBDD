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

#ifndef HQSPRE_GATE_HPP_
#define HQSPRE_GATE_HPP_

#include <initializer_list>
#include <iosfwd>
#include <vector>

#include "aux.hpp"
#include "literal.hpp"

/**
 * \file gate.hpp
 * \brief Contains the declaration of gate-related data structures.
 * \author Ralf Wimmer
 * \date 04/2017
 */

namespace hqspre {

/**
 * \brief Type of a gate.
 */
enum class GateType: unsigned int
{
    AND_GATE, ///< AND-gate
    XOR_GATE, ///< XOR-gate
    MUX_GATE, ///< Multiplexer
    UNKNOWN   ///< unknown gate, typically from semantic gate detection
};


/**
 * \brief Representation of a gate extracted from a formula in CNF.
 */
class Gate
{
public:
    /**
     * \brief Constructs a gate without inputs and unknown type.
     *
     * \note Gates with no encoding clauses are considered invalid. Gates with type
     *       UNKNOWN can be valid; they may result from semantic gate detection.
     */
    Gate():
        _type(GateType::UNKNOWN),
        _output_literal(0),
        _input_literals(),
        _encoding_clauses()
    { }

    explicit Gate(GateType type, Literal output = 0):
        _type(type),
        _output_literal(output),
        _input_literals(),
        _encoding_clauses()
    { }

    Gate(GateType type, Literal output, const std::vector<Literal>& inputs, const std::vector<ClauseID>& clauses):
        _type(type),
        _output_literal(output),
        _input_literals(inputs),
        _encoding_clauses(clauses)
    { }

    Gate(GateType type, Literal output, std::vector<Literal>&& inputs, std::vector<ClauseID>&& clauses):
        _type(type),
        _output_literal(output),
        _input_literals(std::move(inputs)),
        _encoding_clauses(std::move(clauses))
    { }

    Gate(GateType type, Literal output, const std::initializer_list<Literal>& inputs, const std::initializer_list<ClauseID>& clauses):
        _type(type),
        _output_literal(output),
        _input_literals(inputs),
        _encoding_clauses(std::move(clauses))
    { }

    GateType _type;                          ///< Type of the gate (AND/XOR)
    Literal _output_literal;                 ///< output literal (might be negated for NAND/EQUIV)
    std::vector<Literal> _input_literals;    ///< inputs of the gate (might be negated)
    std::vector<ClauseID> _encoding_clauses; ///< IDs of the clauses encoding the gate
};

std::ostream& operator<<(std::ostream& stream, const Gate& gate);


/**
 * \brief Class for managing gates found in a formula in CNF
 */
class GateInfo
{
public:
    //@{
    /**
     * \name Construction, assignment, and destruction
     */
    GateInfo() = default;
    GateInfo(const GateInfo& other) = default;
    GateInfo(GateInfo&& other) = default;
    ~GateInfo() noexcept = default;
    GateInfo& operator=(const GateInfo& other) = default;
    GateInfo& operator=(GateInfo&& other) = default;
    //@}

    //@{
    /**
     * \name Getting iterators for traversing the gates
     */
    using iterator = std::vector<Gate>::iterator;
    using const_iterator = std::vector<Gate>::const_iterator;
    using reverse_iterator = std::vector<Gate>::reverse_iterator;
    using const_reverse_iterator = std::vector<Gate>::const_reverse_iterator;

    const_iterator cbegin()          const noexcept { return _gates.cbegin();  }
    const_iterator cend()            const noexcept { return _gates.cend();    }
    const_iterator begin()           const noexcept { return _gates.begin();   }
    const_iterator end()             const noexcept { return _gates.end();     }
    iterator begin()                 noexcept       { return _gates.begin();   }
    iterator end()                   noexcept       { return _gates.end();     }

    const_reverse_iterator crbegin() const noexcept { return _gates.crbegin(); }
    const_reverse_iterator crend()   const noexcept { return _gates.crend();   }
    const_reverse_iterator rbegin()  const noexcept { return _gates.rbegin();  }
    const_reverse_iterator rend()    const noexcept { return _gates.rend();    }
    reverse_iterator rbegin()        noexcept       { return _gates.rbegin();  }
    reverse_iterator rend()          noexcept       { return _gates.rend();    }
    //@}

    bool gateValid(const Gate& g) const;
    void invalidateGate(Gate& g);

    /**
     * \brief Marks a clause as modified since the last gate detection.
     *
     * This invalidates all gates that are encoded my this clause.
     */
    void touchClause(const ClauseID c_nr) noexcept
    {
        val_assert(c_nr < _clause_modified.size());
        val_assert(_clause_modified.size() == _gate_clause.size());

        _clause_modified[c_nr] = true;
    }

    void update();
    void sort();
    void clear();

    bool checkConsistency() const;

    /**
     * \brief Returns the number of gates (including invalid ones) in the database
     */
    std::size_t numGates() const noexcept { return _gates.size(); }

    std::size_t numValidGates() const;
    void addGate(const Gate& gate);
    void addGate(Gate&& gate);

    /**
     * \brief Provides read access to the gate database.
     */
    const std::vector<Gate>& getGates() const noexcept { return _gates; }

    /**
     * \brief Adjusts the number of clauses in the system.
     *
     * This increases the size of the clause-related internal data structures in the
     * gate handler appropriately.
     */
    void resizeClauses(const std::size_t num_clauses)
    {
        if (num_clauses > _clause_modified.size()) {
            _clause_modified.resize(num_clauses, true);
            _gate_clause.resize(num_clauses, false);
        }
    }

    /**
     * \brief Adjusts the number of variables in the system.
     *
     * This increases the size of the variable-related internal data structures in the
     * gate handler appropriately.
     */
    void resizeVars(const Variable max_var_index) {
        if (max_var_index >= _gate_output.size()) {
            _gate_output.resize(max_var_index + 1, false);
            _gate_input.resize(max_var_index + 1, 0u);
        }
    }

    /**
     * \brief Returns true iff the given variable is defined by a gate in the database
     * \note A return value of true does not necessarily mean that the corresponding gate
     *       is still valid. An encoding clause might have been modified in the mean time.
     *       (If the gate was invalidated by a call to GateInfo::invalidateGate(Gate&), then
     *       this function takes that into account and returns false for the gate output).
     */
    bool isGateOutput(const Variable var) const noexcept {
        val_assert(var < _gate_output.size());
        return _gate_output[var];
    }


    /**
     * \brief Returns the number of gate inputs the variable is connected to
     * \note The return value does not take into account if the corresponding gate
     *       is still valid. An encoding clause might have been modified in the mean time.
     *       (If the gate was invalidated by a call to GateInfo::invalidateGate(Gate&), then
     *       this function takes that into account and returns the correspondingly decreased value).
     */
    unsigned int isGateInput(const Variable var) const noexcept {
        val_assert(var < _gate_input.size());
        return _gate_input[var];
    }

    /**
     * \brief Returns true iff the given clause is used in the definition of a gate in the database
     * \note A return value of true does not necessarily mean that the corresponding gate
     *       is still valid. An encoding clause might have been modified in the mean time.
     *       (If the gate was invalidated by a call to GateInfo::invalidateGate(Gate&), then
     *       this function takes that into account and returns false for the encoding clause).
     */
    bool isGateClause(const ClauseID c_nr) const noexcept {
        val_assert(c_nr < _gate_clause.size());
        return _gate_clause[c_nr];
    }

private:
    std::vector<Gate> _gates;              ///< \brief The list of gates in the database
    std::vector<bool> _clause_modified;    ///< \brief Which clauses have been modified since the last call to GateInfo::update()
    std::vector<bool> _gate_clause;        ///< \brief Is a clause used to encode a gate in the database
    std::vector<bool> _gate_output;        ///< \brief Is a variable the output of a gate in the database
    std::vector<unsigned int> _gate_input; ///< \brief For each variable the number of gate the variable is an input of
};


} // end namespace hqspre

#endif
