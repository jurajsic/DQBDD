/*
 * This file is part of HQSpre.
 *
 * Copyright 2016-2018 Ralf Wimmer, Sven Reimer, Paolo Marin, Bernd Becker
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

#ifndef HQSPRE_PREFIX_IPP_
#define HQSPRE_PREFIX_IPP_

/**
 * \file prefix.ipp
 * \brief Inline functions for manipulating a formula's quantifier prefix
 * \author Ralf Wimmer
 * \date 06/2016, modified 11/2018
 */

namespace hqspre {

/**
 * \brief Returns true if the variable is universal and has not been deleted.
 */
inline bool
Prefix::isUniversal(Variable var) const noexcept
{
    val_assert(minVarIndex() <= var);

    if (var > maxVarIndex())
        return false;
    else
        return _var_status[var] == VariableStatus::UNIVERSAL;
}

/**
 * \brief Returns true if the variable is existential and has not been deleted.
 */
inline bool
Prefix::isExistential(Variable var) const noexcept
{
    val_assert(minVarIndex() <= var);

    if (var > maxVarIndex())
        return false;
    else
        return _var_status[var] == VariableStatus::EXISTENTIAL;
}

/**
 * \brief Returns true if the variable has been deleted.
 */
inline bool
Prefix::varDeleted(Variable var) const noexcept
{
    val_assert(minVarIndex() <= var);

    if (var > maxVarIndex())
        return true;
    else
        return _var_status[var] == VariableStatus::DELETED;
}

inline VariableStatus
Prefix::getStatus(const Variable var) const noexcept
{
    val_assert(minVarIndex() <= var);

    if (var > maxVarIndex())
        return VariableStatus::DELETED;
    else
        return _var_status[var];
}

//--------------------------------------------------------

/**
 * \brief Creates a copy of the given variable.
 *
 * The copy inherits the dependencies of the given variable.
 * \return the index of the newly created variable.
 */
inline void
Prefix::makeCopy(const Variable var, const Variable original)
{
    val_assert(varDeleted(var));
    val_assert(!varDeleted(original));

    if (var > maxVarIndex()) setMaxVarIndex(var);
    if (isExistential(original))
        Prefix::addEVar(var);
    else if (isUniversal(original))
        Prefix::addUVar(var);
}

inline void
DQBFPrefix::makeCopy(const Variable var, const Variable original)
{
    val_assert(!varDeleted(original));
    Prefix::makeCopy(var, original);
    if (_in_rmb[original])
        moveToRMB(var);
    else
        setDependencies(var, _dependencies[original]);
}

inline void
QBFPrefix::makeCopy(const Variable var, const Variable original)
{
    val_assert(!varDeleted(original));
    Prefix::makeCopy(var, original);
    _depth[var] = _depth[original];
    _blocks[_depth[original]].insert(var);
}

inline void
SATPrefix::makeCopy(const Variable var, const Variable original)
{
    val_assert(!varDeleted(original));
    Prefix::makeCopy(var, original);
    _vars.insert(var);
}

//--------------------------------------------------------

/**
 * \brief Clears the prefix and removes all variables.
 */
inline void
Prefix::clear()
{
    _var_status.resize(1);
    _num_e_vars = 0;
    _num_u_vars = 0;
}

/**
 * \brief Clears the prefix and removes all variables.
 */
inline void
QBFPrefix::clear()
{
    Prefix::clear();
    _depth.resize(1);
    _blocks.clear();
}

/**
 * \brief Clears the prefix and removes all variables.
 */
inline void
DQBFPrefix::clear()
{
    Prefix::clear();
    _dependencies.resize(1);
    _rmb.clear();
    _in_rmb.resize(1);
    _univ_vars.clear();
}

/**
 * \brief Clears the prefix and removes all variables.
 */
inline void
SATPrefix::clear()
{
    Prefix::clear();
    _vars.clear();
}

//--------------------------------------------------------

/**
 * \brief Creates a new existential variable
 * \param var the index of the variable that should be created
 */
inline void
Prefix::addEVar(const Variable var)
{
    val_assert(varDeleted(var));

    _var_status[var] = VariableStatus::EXISTENTIAL;
    ++_num_e_vars;
}

inline void
DQBFPrefix::addEVar(const Variable var)
{
    val_assert(varDeleted(var));
    Prefix::addEVar(var);
    moveToRMB(var);
}

inline void
SATPrefix::addEVar(const Variable var)
{
    val_assert(varDeleted(var));
    Prefix::addEVar(var);
    _vars.insert(var);
}

inline void
DQBFPrefix::addEVar(const Variable var, const std::set<Variable>& dependencies)
{
    Prefix::addEVar(var);
    if (!dependencies.empty()) setDependencies(var, dependencies);
}

inline void
DQBFPrefix::addEVar(const Variable var, std::set<Variable>&& dependencies)
{
    Prefix::addEVar(var);
    if (!dependencies.empty()) setDependencies(var, std::move(dependencies));
}

inline void
QBFPrefix::addEVar(const Variable var)
{
    val_assert(varDeleted(var));
    Prefix::addEVar(var);

    std::size_t var_depth = 0;
    if (_blocks.empty()) {
        var_depth = 0;
        _blocks.push_back({var});
    } else if (_blocks.back().empty() || isExistential(*(_blocks.back().begin()))) {
        var_depth = _blocks.size() - 1;
        _blocks[var_depth].insert(var);
    } else {
        var_depth = _blocks.size();
        _blocks.push_back({var});
    }
    if (_depth.size() <= var) _depth.resize(var + 1);
    _depth[var] = var_depth;
}

//--------------------------------------------------------

/**
 * \brief Creates a new universal variable.
 * \param var the index of the variable that should be created.
 */
inline void
Prefix::addUVar(const Variable var)
{
    val_assert(varDeleted(var));

    _var_status[var] = VariableStatus::UNIVERSAL;
    ++_num_u_vars;
}

/**
 * \brief Does nothing as SAT formulas do not have universal variables.
 */
inline void
SATPrefix::addUVar(const Variable)
{
    return;
}

inline void
DQBFPrefix::addUVar(const Variable var)
{
    val_assert(varDeleted(var));

    Prefix::addUVar(var);

    // All variables which were in the RMB so far,
    // become ordinary variables.
    for (Variable evar : _rmb) {
        _in_rmb[evar]       = false;
        _dependencies[evar] = _univ_vars;
    }
    _rmb.clear();

    _univ_vars.insert(var);
}

inline void
QBFPrefix::addUVar(const Variable var)
{
    val_assert(varDeleted(var));

    Prefix::addUVar(var);

    std::size_t var_depth = 0;
    if (_blocks.empty()) {
        var_depth = _blocks.size();
        _blocks.push_back({var});
    } else if (_blocks.back().empty() || isUniversal(*(_blocks.back().begin()))) {
        var_depth = _blocks.size() - 1;
        _blocks[var_depth].insert(var);
    } else {
        var_depth = _blocks.size();
        _blocks.push_back({var});
    }
    if (_depth.size() <= var) _depth.resize(var + 1);
    _depth[var] = var_depth;
}

//--------------------------------------------------------

/**
 * \brief Sets the maximal index of the variables that will be created.
 * \param index the maximal index
 *
 * The function is used to resize the data structures such that all
 * occurring variables can be stored without further size adaptations.
 * Creating more variables than specified is allowed (it is just less
 * efficient).
 */
inline void
Prefix::setMaxVarIndex(const Variable var)
{
    if (var > maxVarIndex()) {
        _var_status.resize(var + 1, VariableStatus::DELETED);
    }
}

inline void
DQBFPrefix::setMaxVarIndex(const Variable var)
{
    if (var > maxVarIndex()) {
        Prefix::setMaxVarIndex(var);
        _dependencies.resize(var + 1);
        _in_rmb.resize(var + 1, false);
    }
}

inline void
QBFPrefix::setMaxVarIndex(const Variable var)
{
    if (var > maxVarIndex()) {
        Prefix::setMaxVarIndex(var);
        _depth.resize(var + 1, 0);
    }
}

//--------------------------------------------------------

/**
 * \brief Removes a given variable from the prefix.
 * \param var the index of the variable that should be removed
 */
inline void
Prefix::removeVar(const Variable var)
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    if (varDeleted(var)) return;

    if (isExistential(var))
        --_num_e_vars;
    else if (isUniversal(var))
        --_num_u_vars;
    _var_status[var] = VariableStatus::DELETED;
}

inline void
SATPrefix::removeVar(const Variable var)
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    if (varDeleted(var)) return;
    _vars.erase(var);
    Prefix::removeVar(var);
}

/**
 * \brief Removes a given variable from the prefix.
 * \param var the index of the the variable that should be removed
 */
inline void
DQBFPrefix::removeVar(const Variable var)
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    if (varDeleted(var)) return;

    if (isUniversal(var)) _univ_vars.erase(var);

    if (_in_rmb[var]) {
        _in_rmb[var] = false;
        _rmb.erase(var);
        for (Variable univ : _univ_vars) _dependencies[univ].erase(var);
    } else {
        for (Variable dep : _dependencies[var]) {
            _dependencies[dep].erase(var);
        }
        _dependencies[var].clear();
    }

    Prefix::removeVar(var);
}

/**
 * \brief Removes a given variable from the prefix.
 * \param var the index of the variable that should be removed
 */
inline void
QBFPrefix::removeVar(const Variable var)
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    if (varDeleted(var)) return;

    const std::size_t var_depth = _depth[var];
    _blocks[var_depth].erase(var);
    _depth[var] = 0;
    Prefix::removeVar(var);
}

//--------------------------------------------------------

/**
 * \brief Returns the type of the prefix.
 */
inline PrefixType
DQBFPrefix::type() const noexcept
{
    return PrefixType::DQBF;
}

/**
 * \brief Returns the type of the prefix.
 */
inline PrefixType
QBFPrefix::type() const noexcept
{
    return PrefixType::QBF;
}

/**
 * \brief Returns the type of the prefix.
 */
inline PrefixType
SATPrefix::type() const noexcept
{
    return PrefixType::SAT;
}

//--------------------------------------------------------

/**
 * \brief Checks if one variable depends on another one.
 * \pre One variable must be universal, the other one existential
 * \param var1 one of the variables
 * \param var2 the other variable
 * \return true if the existential variable depends on the universal one
 * \note The order in which the variables are given does not matter.
 */
inline bool
DQBFPrefix::depends(const Variable var1, const Variable var2) const
{
    val_assert(minVarIndex() <= var1 && var1 <= maxVarIndex());
    val_assert(minVarIndex() <= var2 && var2 <= maxVarIndex());
    val_assert(!varDeleted(var1) && !varDeleted(var2));
    val_assert(isUniversal(var1) != isUniversal(var2));

    if (_in_rmb[var1] || _in_rmb[var2])
        return true;
    else
        return _dependencies[var1].find(var2) != _dependencies[var1].cend();
}

/**
 * \brief Checks if one variable depends on another one.
 * \pre One variable must be universal, the other one existential
 * \param var1 one of the variables
 * \param var2 the other variable
 * \return true if the existential variable depends on the universal one,
 *           i.e., if the existential is right of the universal variable.
 * \note The order in which the variables are given does not matter.
 */
inline bool
QBFPrefix::depends(const Variable var1, const Variable var2) const
{
    val_assert(minVarIndex() <= var1 && var1 <= maxVarIndex());
    val_assert(minVarIndex() <= var2 && var2 <= maxVarIndex());
    val_assert(!varDeleted(var1) && !varDeleted(var2));
    val_assert(isUniversal(var1) != isUniversal(var2));

    if (isUniversal(var1) && getLevel(var1) < getLevel(var2)) return true;
    if (isUniversal(var2) && getLevel(var2) < getLevel(var1)) return true;
    return false;
}

inline bool
SATPrefix::depends(const Variable, const Variable) const
{
    return false;
}

//--------------------------------------------------------

inline std::size_t
QBFPrefix::getLevel(const Variable var) const noexcept
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());

    return _depth[var];
}

//--------------------------------------------------------

/**
 * \brief Removes a dependency between two variables.
 *
 * If the variables do not depend on each other, this operation
 * does nothing.
 * \pre One variable must be universal, the other one existential.
 */
inline void
DQBFPrefix::removeDependency(const Variable var1, const Variable var2)
{
    val_assert(minVarIndex() <= var1 && var1 <= maxVarIndex());
    val_assert(isUniversal(var1) != isUniversal(var2));
    val_assert(!varDeleted(var1) && !varDeleted(var2));

    if (isExistential(var1) && _in_rmb[var1]) {
        _dependencies[var1] = _univ_vars;
        _rmb.erase(var1);
        _in_rmb[var1] = false;
    } else if (isExistential(var2) && _in_rmb[var2]) {
        _dependencies[var2] = _univ_vars;
        _rmb.erase(var2);
        _in_rmb[var2] = false;
    }
    _dependencies[var1].erase(var2);
    _dependencies[var2].erase(var1);
}

/**
 * \brief Adds a dependency between two variables.
 *
 * If the variables already depend on each other, this operation
 * does nothing.
 * \pre One variable must be universal, the other one existential.
 */
inline void
DQBFPrefix::addDependency(Variable var1, Variable var2)
{
    val_assert(minVarIndex() <= var1 && var1 <= maxVarIndex());
    val_assert(minVarIndex() <= var2 && var2 <= maxVarIndex());
    val_assert(!varDeleted(var1) && !varDeleted(var2));
    val_assert(isUniversal(var1) != isUniversal(var2));

    _dependencies[var1].insert(var2);
    _dependencies[var2].insert(var1);

    if (isExistential(var1) && _dependencies[var1].size() == numUVars()) {
        _in_rmb[var1] = true;
        _rmb.insert(var1);
        _dependencies[var1].clear();
    } else if (isExistential(var2) && _dependencies[var2].size() == numUVars()) {
        _in_rmb[var2] = true;
        _rmb.insert(var2);
        _dependencies[var2].clear();
    }
}

inline const std::set<Variable>&
DQBFPrefix::getDependencies(const Variable var) const noexcept
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());

    if (_in_rmb[var])
        return _univ_vars;
    else
        return _dependencies[var];
}

inline const std::set<Variable>&
QBFPrefix::getVarBlock(const std::size_t level) const noexcept
{
    val_assert(level < _blocks.size());
    return _blocks[level];
}

}  // end namespace hqspre

#endif
