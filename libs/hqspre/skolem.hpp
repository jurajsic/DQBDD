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

#ifndef HQSPRE_SKOLEM_HPP_
#define HQSPRE_SKOLEM_HPP_

#include <map>
#include <vector>
#include "auxil.hpp"
#include "clause.hpp"
#include "gate.hpp"
#include "literal.hpp"

namespace hqspre {

/**
 * \brief Wrapper class for creating representations of boolean functions.
 *
 * If the computation of Skolem functions is requested, an instantiation of
 * BoolFunctionManager is needed which creates the Skolem functions in form
 * of BDDs/AIGs/clauses/...
 * BoolFunctionManager is also responsible for storing the table of Skolem
 * functions as well as a buffer which can be used to store intermediate
 * results.
 */
class BoolFunctionManager
{
   public:
    BoolFunctionManager()                           = default;
    BoolFunctionManager(const BoolFunctionManager&) = default;
    BoolFunctionManager(BoolFunctionManager&&)      = default;
    virtual ~BoolFunctionManager()                  = default;
    BoolFunctionManager& operator=(const BoolFunctionManager&) = default;
    BoolFunctionManager& operator=(BoolFunctionManager&&) = default;

    /**
     * \brief Writes the constant 0 function either into the Skolem table or the
     * buffer. \param var position into which the constant 0 function should be
     * written \param buffer if true, the position is taken in the buffer,
     * otherwise in the Skolem table
     */
    virtual void Const0(Variable var, bool buffer) = 0;

    /**
     * \brief Writes the constant 1 function either into the Skolem table or the
     * buffer. \param var position into which the constant 1 function should be
     * written \param buffer if true, the position is taken in the buffer,
     * otherwise in the Skolem table
     */
    virtual void Const1(Variable var, bool buffer) = 0;

    /**
     * \brief Copies a Boolean function from one entry in the buffer or Skolem
     * function table to another one.
     *
     * \param origin the source location of the Boolean function
     * \param buffer_origin if true the origin location is in the buffer,
     * otherwise in the Skolem function table \param target the target location
     * for the copy \param buffer_target if true the target location is in the
     * buffer, otherwise int he Skolem function table. \pre If the buffer is
     * accessed, it must be large enough such that the accessed position exists.
     * \pre If the source location is in the Skolem table, the entry must exist.
     */
    virtual void copy(Variable origin, bool buffer_origin, Variable target, bool buffer_target) = 0;
    virtual void Not(Variable origin, bool buffer_origin, Variable target, bool buffer_target)  = 0;
    virtual void And(Variable origin1, bool buffer_origin1, Variable origin2, bool buffer_origin2, Variable target,
                     bool buffer_target)
        = 0;
    virtual void Or(Variable origin1, bool buffer_origin1, Variable origin2, bool buffer_origin2, Variable target,
                    bool buffer_target)
        = 0;
    virtual void Xor(Variable origin1, bool buffer_origin1, Variable origin2, bool buffer_origin2, Variable target,
                     bool buffer_target)
        = 0;
    virtual void Ite(Variable sel, bool buffer_sel, Variable input1, bool buffer_input1, Variable input0,
                     bool buffer_input0, Variable target, bool buffer_target)
        = 0;

    /**
     * \brief Removes all entries from the buffer.
     *
     * Afterwards the buffer has size 0, i.e., before the buffer can be used
     * again, its size must be increased.
     * \sa BoolFunctionManager::resizeBuffer(unsigned int)
     */
    virtual void clearBuffer()                   = 0;
    virtual void resizeBuffer(unsigned int size) = 0;

    /**
     * \brief Removes an entry from the Skolem function table.
     *
     * The entry corresponding to the given variable `var` is removed.
     */
    virtual void removeEntry(Variable var) = 0;

    /**
     * \brief Checks if a given entry exists in the Skolem table.
     * \return true if `var` is contained in the Skolem table.
     */
    virtual bool existsEntry(Variable var) const = 0;
};

#ifdef SKOLEM_CUDD

#    include <cuddObj.hh>

class CuddManager : public BoolFunctionManager
{
   public:
    explicit CuddManager(Cudd* manager) : _manager(manager) { val_assert(_manager != nullptr); }

    CuddManager(const CuddManager&) = default;
    CuddManager(CuddManager&&)      = default;
    ~CuddManager()                  = default;
    CuddManager& operator=(const CuddManager&) = default;
    CuddManager& operator=(CuddManager&&) = default;

    virtual void variable(Variable var, const BDD& var_bdd) { skolem_table[var] = var_bdd; }

    virtual void Const0(Variable var, bool buffer) override
    {
        val_assert(!buffer || var < buffer.size());

        if (buffer)
            buffer[var] = _manager->bddZero();
        else
            skolem_table[var] = _manager->bddZero();
    }

    virtual void Const1(Variable var, bool buffer) override
    {
        val_assert(!buffer || var < buffer.size());

        if (buffer)
            buffer[var] = _manager->bddOne();
        else
            skolem_table[var] = _manager->bddOne();
    }

    virtual void Not(Variable origin, bool buffer_origin, Variable target, bool buffer_target) override
    {
        val_assert(!buffer_origin || origin < buffer.size());
        val_assert(buffer_origin || skolem_table.find(origin) != skolem_table.end());
        val_assert(!buffer_target || target < buffer.size());

        (buffer_target ? buffer[target] : skolem_table[target])
            = (buffer_origin ? (!buffer[origin]) : (!skolem_table[origin]));
    }

    virtual void And(Variable origin1, bool buffer_origin1, Variable origin2, bool buffer_origin2, Variable target,
                     bool buffer_target) override
    {
        val_assert(!buffer_origin1 || origin1 < buffer.size());
        val_assert(buffer_origin1 || skolem_table.find(origin1) != skolem_table.end());
        val_assert(!buffer_origin2 || origin2 < buffer.size());
        val_assert(buffer_origin2 || skolem_table.find(origin2) != skolem_table.end());
        val_assert(!buffer_target || target < buffer.size());

        (buffer_target ? buffer[target] : skolem_table[target])
            = (buffer_origin1 ? buffer[origin1] : skolem_table[origin1])
              * (buffer_origin2 ? buffer[origin2] : skolem_table[origin2]);
    }

    virtual void Or(Variable origin1, bool buffer_origin1, Variable origin2, bool buffer_origin2, Variable target,
                    bool buffer_target) override
    {
        val_assert(!buffer_origin1 || origin1 < buffer.size());
        val_assert(buffer_origin1 || skolem_table.find(origin1) != skolem_table.end());
        val_assert(!buffer_origin2 || origin2 < buffer.size());
        val_assert(buffer_origin2 || skolem_table.find(origin2) != skolem_table.end());
        val_assert(!buffer_target || target < buffer.size());

        (buffer_target ? buffer[target] : skolem_table[target])
            = (buffer_origin1 ? buffer[origin1] : skolem_table[origin1])
              + (buffer_origin2 ? buffer[origin2] : skolem_table[origin2]);
    }

    virtual void Xor(Variable origin1, bool buffer_origin1, Variable origin2, bool buffer_origin2, Variable target,
                     bool buffer_target) override
    {
        val_assert(!buffer_origin1 || origin1 < buffer.size());
        val_assert(buffer_origin1 || skolem_table.find(origin1) != skolem_table.end());
        val_assert(!buffer_origin2 || origin2 < buffer.size());
        val_assert(buffer_origin2 || skolem_table.find(origin2) != skolem_table.end());
        val_assert(!buffer_target || target < buffer.size());

        (buffer_target ? buffer[target] : skolem_table[target])
            = (buffer_origin1 ? buffer[origin1] : skolem_table[origin1])
              ^ (buffer_origin2 ? buffer[origin2] : skolem_table[origin2]);
    }

    virtual void Ite(Variable sel, bool buffer_sel, Variable input1, bool buffer_input1, Variable input0,
                     bool buffer_input0, Variable target, bool buffer_target) override
    {
        val_assert(!buffer_sel || sel < buffer.size());
        val_assert(!buffer_input1 || input1 < buffer.size());
        val_assert(!buffer_input0 || input0 < buffer.size());
        val_assert(!buffer_target || target < buffer.size());
        val_assert(buffer_sel || skolem_table.find(sel) != skolem_table.end());
        val_assert(buffer_input1 || skolem_table.find(input1) != skolem_table.end());
        val_assert(buffer_input0 || skolem_table.find(input0) != skolem_table.end());

        (buffer_target ? buffer[target] : skolem_table[target])
            = (buffer_sel ? buffer[sel] : skolem_table[sel])
                  .Ite((buffer_input1 ? buffer[input1] : skolem_table[input1]),
                       (buffer_input0 ? buffer[input0] : skolem_table[input0]));
    }

    virtual void clearBuffer() override { buffer.clear(); }

    virtual void resizeBuffer(unsigned int size) override { buffer.resize(size); }

    virtual void removeEntry(Variable var) override { skolem_table.erase(var); }

    virtual bool existsEntry(Variable var) const override { return skolem_table.find(var) != skolem_table.cend(); }

    virtual void copy(Variable origin, bool buffer_origin, Variable target, bool buffer_target) override
    {
        val_assert(!buffer_origin || origin < buffer.size());
        val_assert(!buffer_target || target < buffer.size());
        val_assert(buffer_origin || skolem_table.find(origin) != skolem_table.end());

        (buffer_target ? buffer[target] : skolem_table[target])
            = (buffer_origin ? buffer[origin] : skolem_table[origin]);
    }

    Cudd*                   getManager() { return _manager; }
    std::map<Variable, BDD> getSkolemTable() const { return skolem_table; }

   private:
    Cudd*                   _manager;      ///< Pointer to the actual Cudd manager
    std::map<Variable, BDD> skolem_table;  ///< Table with all Skolem functions
    std::vector<BDD>        buffer;        ///< Table with temporary BDDs
};

#endif

class SkolemEntry
{
   public:
    SkolemEntry()          = default;
    virtual ~SkolemEntry() = default;

    virtual void computeSkolem(BoolFunctionManager* manager) = 0;
};

class SkolemUnitPure : public SkolemEntry
{
   public:
    SkolemUnitPure(Variable _var, bool _isNegated) : SkolemUnitPure(var2lit(_var, _isNegated)) {}

    explicit SkolemUnitPure(Literal _literal) : SkolemEntry(), literal(_literal) {}

    virtual ~SkolemUnitPure() = default;

    virtual void computeSkolem(BoolFunctionManager* manager) override
    {
        const Variable var = lit2var(literal);
        if (isNegative(literal))
            manager->Const0(var, false);
        else
            manager->Const1(var, false);
    }

   private:
    Literal literal;
};

class SkolemEquiv : public SkolemEntry
{
   public:
    SkolemEquiv(Variable _variable, Literal _replaced_by) :
        SkolemEntry(),
        variable(_variable),
        replaced_by(_replaced_by)
    {}

    virtual ~SkolemEquiv() = default;

    virtual void computeSkolem(BoolFunctionManager* manager) override
    {
        const Variable replaced_by_var = lit2var(replaced_by);

        if (isNegative(replaced_by)) {
            manager->Not(replaced_by_var, false, variable, false);
        } else {
            manager->copy(replaced_by_var, false, variable, false);
        }
    }

   private:
    Variable variable;
    Literal  replaced_by;
};

class SkolemDepElim : public SkolemEntry
{
   public:
    SkolemDepElim(Literal _uLit, Literal _eLit_old, Literal _eLit0, Literal _eLit1) :
        SkolemEntry(),
        uLit(_uLit),
        eLit_old(_eLit_old),
        eLit0(_eLit0),
        eLit1(_eLit1)
    {
        val_assert(isPositive(uLit));
        val_assert(isPositive(eLit_old));
        val_assert(eLit0 <= 1 || isPositive(eLit0));
        val_assert(eLit1 <= 1 || isPositive(eLit1));
    }

    virtual ~SkolemDepElim() = default;

    virtual void computeSkolem(BoolFunctionManager* manager) override
    {
        // TODO: Take into account that eLit0 and elit1 might be constants
        //        manager->Ite(uLit, false, eLit1, false, eLit0, false, eLit_old,
        //        false);
        if (eLit0 != eLit_old) {
            manager->removeEntry(lit2var(eLit0));
        }
        if (eLit1 != eLit_old) {
            manager->removeEntry(lit2var(eLit1));
        }
    }

   private:
    Literal uLit;      ///< Universal variable that was deleted from the dependency set
                       ///< of eVar_old.
    Literal eLit_old;  ///< Existential variable that contains uVar in its
                       ///< dependency set.
    Literal eLit0;     ///< Copy of the existential variable eVar_old in the 0-cofactor.
    Literal eLit1;     ///< Copy of the existential variable eVar_old in the 1-cofactor.
};

class SkolemUniversalExpansion : public SkolemEntry
{
   public:
    SkolemUniversalExpansion(Variable _uVar, const std::vector<Variable>& _vars_old,
                             const std::vector<Variable>& _vars0, const std::vector<Variable>& _vars1) :
        SkolemEntry(),
        uVar(_uVar),
        vars_old(_vars_old),
        vars0(_vars0),
        vars1(_vars1)
    {}

    virtual ~SkolemUniversalExpansion() = default;

    virtual void computeSkolem(BoolFunctionManager* manager) override
    {
        for (unsigned int i = 0; i < vars_old.size(); ++i) {
            val_assert(vars0[i] != vars1[i]);
            manager->Ite(uVar, false, vars1[i], false, vars0[i], false, vars_old[i], false);
            if (vars_old[i] != vars0[i]) {
                manager->removeEntry(vars0[i]);
            }
            if (vars_old[i] != vars1[i]) {
                manager->removeEntry(vars1[i]);
            }
        }
    }

   private:
    Variable              uVar;      ///< Eliminated universal variable;
    std::vector<Variable> vars_old;  ///< original existential variables depending on uVar
    std::vector<Variable> vars0;     ///< copies of the dependent variables in the 0-cofactor
    std::vector<Variable> vars1;     ///< copies of the dependent variables in the 1-cofactor
};

class SkolemGate : public SkolemEntry
{
   public:
    explicit SkolemGate(const Gate& _gate) : SkolemEntry(), gate(_gate) {}

    virtual ~SkolemGate() {}

    virtual void computeSkolem(BoolFunctionManager* manager) override
    {
        const unsigned int num_inputs = gate.input_literals.size();
        manager->resizeBuffer(num_inputs);

        Variable index = 0;
        for (const auto input : gate.input_literals) {
            const Variable input_var = lit2var(input);

            if (isNegative(input)) {
                manager->Not(input_var, false, index, true);
            } else {
                manager->copy(input_var, false, index, true);
            }
            ++index;
        }

        switch (gate.type) {
            case GateType::AND_GATE:
                // AND gate with an arbitrary number of inputs
                for (Variable pos = 1; pos < num_inputs; ++pos) {
                    manager->And(0, true, pos, true, 0, true);
                }
                if (isNegative(gate.output_literal)) {
                    manager->Not(0, true, lit2var(gate.output_literal), false);
                } else {
                    manager->copy(0, true, lit2var(gate.output_literal), false);
                }
                break;

            case GateType::XOR_GATE:
                // XOR gate with two inputs
                val_assert_msg(num_inputs == 2, "Only XOR gates with 2 inputs are supported!");
                if (isNegative(gate.output_literal)) {
                    manager->Xor(0, true, 1, true, 0, true);
                    manager->Not(0, true, lit2var(gate.output_literal), false);
                } else {
                    manager->Xor(0, true, 1, true, lit2var(gate.output_literal), false);
                }
                break;

            case GateType::MUX_GATE:
                // Multiplexer: first input is select, second input the value if select=1,
                // third input the value if select=0.
                val_assert_msg(num_inputs == 3, "A multiplexer needs exactly 3 inputs!");
                if (isNegative(gate.output_literal)) {
                    manager->Ite(0, true, 1, true, 2, true, 0, true);
                    manager->Not(0, true, lit2var(gate.output_literal), false);
                } else {
                    manager->Ite(0, true, 1, true, 2, true, lit2var(gate.output_literal), false);
                }
                break;

            default:
                val_assert_msg(false, "Invalid gate type!");
                break;
        }
        manager->clearBuffer();
    }

   private:
    Gate gate;  ///< The detected gate.
};

class SkolemResolution : public SkolemEntry
{
   public:
    SkolemResolution(Variable _eVar, const std::vector<Clause>& _clauses_pos, const std::vector<Clause>& _clauses_neg) :
        eVar(_eVar),
        clauses_pos(_clauses_pos),
        clauses_neg(_clauses_neg)
    {}

    SkolemResolution(Variable _eVar, std::vector<Clause>&& _clauses_pos, std::vector<Clause>&& _clauses_neg) :
        eVar(_eVar),
        clauses_pos(std::move(_clauses_pos)),
        clauses_neg(std::move(_clauses_neg))
    {}

    virtual void computeSkolem(BoolFunctionManager* /* manager */) override {}

   private:
    Variable            eVar;  ///< The eliminated variable
    std::vector<Clause> clauses_pos;
    std::vector<Clause> clauses_neg;
};

}  // end namespace hqspre

#endif /* HQSPRE_SKOLEM_HPP_ */
