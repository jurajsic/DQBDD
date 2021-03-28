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

#include "prefix.hpp"

#include <algorithm>
#include <iterator>
#include <map>
#include <ostream>
#include <utility>

#include <easylogging++.hpp>
#include "auxil.hpp"

/**
 * \file prefix.cpp
 * \brief Functions for manipulating the quantifier prefix of a (D)QBF
 * \author Ralf Wimmer
 * \date 06/2016
 */

namespace hqspre {

/**
 * \brief Returns the existential variables in increasing order.
 */
std::vector<Variable>
Prefix::getExistVars() const noexcept
{
    std::vector<Variable> result;
    result.reserve(numEVars());
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (isExistential(var)) result.push_back(var);
    }

    return result;
}

/**
 * \brief Returns the universal variables in increasing order.
 */
std::vector<Variable>
Prefix::getUnivVars() const noexcept
{
    std::vector<Variable> result;
    result.reserve(numUVars());
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (isUniversal(var)) result.push_back(var);
    }

    return result;
}

/**
 * \brief Replaces the dependencies of `var1` by the intersection of
 *        `var1` and `var2`'s dependencies.
 */
void
DQBFPrefix::intersectDependencies(const Variable var1, const Variable var2)
{
    val_assert(minVarIndex() <= var1 && var1 <= maxVarIndex());
    val_assert(minVarIndex() <= var2 && var2 <= maxVarIndex());
    val_assert(isExistential(var1) && isExistential(var2));

    if (var1 == var2 || _in_rmb[var2]) return;

    std::set<Variable> new_deps;
    std::set_intersection(_dependencies[var1].cbegin(), _dependencies[var1].cend(), _dependencies[var2].cbegin(),
                          _dependencies[var2].cend(), std::inserter(new_deps, new_deps.end()));
    setDependencies(var1, std::move(new_deps));
}

/**
 * \brief Replaces the dependencies of `var1` by the intersection of
 *        `var1` and `var2`'s dependencies.
 */
void
QBFPrefix::intersectDependencies(const Variable var1, const Variable var2)
{
    val_assert(minVarIndex() <= var1 && var1 <= maxVarIndex());
    val_assert(minVarIndex() <= var2 && var2 <= maxVarIndex());
    val_assert(isExistential(var1) && isExistential(var2));

    if (_depth[var2] < _depth[var1]) setLevel(var1, _depth[var2]);
}

/*
 * \brief Does nothing as SAT formulas do not have dependencies.
 */
void
SATPrefix::intersectDependencies(const Variable, const Variable)
{
    return;
}

//--------------------------------------------------------

/**
 * \brief Replaces the dependencies of `var1` by the union of
 *        `var1` and `var2`'s dependencies.
 */
void
DQBFPrefix::uniteDependencies(const Variable var1, const Variable var2)
{
    val_assert(minVarIndex() <= var1 && var1 <= maxVarIndex());
    val_assert(minVarIndex() <= var2 && var2 <= maxVarIndex());
    val_assert(isExistential(var1) && isExistential(var2));

    if (var1 == var2 || _in_rmb[var1]) return;

    std::set<Variable> new_deps;
    std::set_union(_dependencies[var1].cbegin(), _dependencies[var1].cend(), _dependencies[var2].cbegin(),
                   _dependencies[var2].cend(), std::inserter(new_deps, new_deps.end()));
    setDependencies(var1, std::move(new_deps));
}

/**
 * \brief Replaces the dependencies of `var1` by the union of
 *        `var1` and `var2`'s dependencies.
 */
void
QBFPrefix::uniteDependencies(const Variable var1, const Variable var2)
{
    val_assert(minVarIndex() <= var1 && var1 <= maxVarIndex());
    val_assert(minVarIndex() <= var2 && var2 <= maxVarIndex());
    val_assert(isExistential(var1) && isExistential(var2));

    setLevel(var1, std::max(_depth[var1], _depth[var2]));
}

/*
 * \brief Does nothing as SAT formulas do not have dependencies.
 */
void
SATPrefix::uniteDependencies(const Variable, const Variable)
{
    return;
}
//--------------------------------------------------------

/**
 * \brief Checks if `var1`'s dependencies are a subset of `var2`'s dependencies.
 *
 * If `var1` is universal, this function returns whether `var2` depends on
 * `var1`. If both variables are existential, this function checks whether the
 * dependency set of `var1` is a subset of `var2`'s dependencies.
 *
 * \note `var2` has to be existential.
 */
bool
DQBFPrefix::dependenciesSubset(const Variable var1, const Variable var2) const
{
    val_assert(minVarIndex() <= var1 && var1 <= maxVarIndex());
    val_assert(minVarIndex() <= var2 && var2 <= maxVarIndex());
    val_assert(!varDeleted(var1) && !varDeleted(var2));
    val_assert(isExistential(var2));

    if (var1 == var2) return true;
    if (_in_rmb[var2]) return true;
    if (_in_rmb[var1] && !_in_rmb[var2]) return false;
    if (isUniversal(var1)) return depends(var2, var1);
    if (_dependencies[var1].size() > _dependencies[var2].size()) return false;

    return std::includes(_dependencies[var2].cbegin(), _dependencies[var2].cend(), _dependencies[var1].cbegin(),
                         _dependencies[var1].cend());
}

/**
 * \brief Checks if `var1`'s dependencies are a subset of `var2`'s dependencies.
 *
 * If `var1` is universal, this function returns whether `var2` depends on
 * `var1`. If both variables are existential, this function checks whether the
 * dependency set of `var1` is a subset of `var2`'s dependencies. If both
 * variables are universal, the function returns whether 'var1' depends on more
 * variables than 'var2' \note `var2` has to be existential.
 */
bool
QBFPrefix::dependenciesSubset(const Variable var1, const Variable var2) const
{
    val_assert(minVarIndex() <= var1 && var1 <= maxVarIndex());
    val_assert(minVarIndex() <= var2 && var2 <= maxVarIndex());
    val_assert(!varDeleted(var1) && !varDeleted(var2));

    return _depth[var1] <= _depth[var2];
}

bool
SATPrefix::dependenciesSubset(const Variable, const Variable) const
{
    return true;
}

//--------------------------------------------------------

bool
Prefix::checkConsistency() const
{
    // Correctness of num_u_vars and num_e_vars
    std::size_t e_vars = 0;
    std::size_t u_vars = 0;
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (isUniversal(var)) ++u_vars;
        if (isExistential(var)) ++e_vars;
    }

    if (numEVars() != e_vars) {
        LOG(ERROR) << "num_e_vars is not correct (actual: " << numEVars() << ", correct: " << e_vars;
        return false;
    }
    if (numUVars() != u_vars) {
        LOG(ERROR) << "num_u_vars is not correct (actual: " << numUVars() << ", correct: " << u_vars;
        return false;
    }

    return true;
}

bool
DQBFPrefix::checkConsistency() const
{
    if (!Prefix::checkConsistency()) return false;

    if (_dependencies.size() != maxVarIndex() + 1 || _in_rmb.size() != maxVarIndex() + 1) {
        LOG(ERROR) << "Data structures have wrong size.";
        return false;
    }

    if (_univ_vars.size() != numUVars()) {
        LOG(ERROR) << "numUVars() = " << numUVars() << " != univ_vars.size() = " << _univ_vars.size();
        return false;
    }
    for (Variable uvar : _univ_vars) {
        if (!isUniversal(uvar)) {
            LOG(ERROR) << "Variable " << uvar << " is not universal, but in univ_vars.";
            return false;
        }
    }

    for (Variable var = 1; var <= maxVarIndex(); ++var) {
        if (_in_rmb[var] && _rmb.find(var) == _rmb.cend()) {
            LOG(ERROR) << "RMB-flag of " << var << " set, but variable is not in rmb.";
            return false;
        } else if (!_in_rmb[var] && _rmb.find(var) != _rmb.cend()) {
            LOG(ERROR) << "RMB-flag of " << var << " not set, but variable is in rmb.";
            return false;
        }
    }

    // Consistency of dependency lists
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (varDeleted(var)) continue;
        if (_in_rmb[var]) {
            if (!isExistential(var)) {
                LOG(ERROR) << "Universal variable " << var << " in DQBF-rmb.";
                return false;
            }
            if (!_dependencies[var].empty()) {
                LOG(ERROR) << "Variable " << var << " is in RMB, but dependency data structure is not empty.";
                return false;
            }
            for (Variable univ : _univ_vars) {
                if (_dependencies[univ].find(var) == _dependencies[univ].cend()) {
                    LOG(ERROR) << "RMB-Variable " << var << " does not appear in the dependency list of " << univ
                               << '.';
                    return false;
                }
            }
            continue;
        }

        // only for non-rmb variables
        for (Variable dep : _dependencies[var]) {
            if (varDeleted(dep)) {
                LOG(ERROR) << "Variable " << var << " depends on " << dep << ", but " << dep << " is deleted.";
                return false;
            }
            if (!_in_rmb[dep] && _dependencies[dep].find(var) == _dependencies[dep].cend()) {
                LOG(ERROR) << "Variable " << var << " depends on " << dep << ", but not vice versa.";
                return false;
            }
        }
    }

    return true;
}

bool
QBFPrefix::checkConsistency() const
{
    if (!Prefix::checkConsistency()) return false;

    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (varDeleted(var)) continue;
        if (_depth[var] >= _blocks.size()) {
            LOG(ERROR) << "Variable " << var << " in a block (" << _depth[var] << ") which does not exist.";
            return false;
        }

        if (_blocks[_depth[var]].find(var) == _blocks[_depth[var]].cend()) {
            LOG(ERROR) << "Variable " << var << " is not in the block in which it should be.";
            return false;
        }
    }

    for (const auto& block : _blocks) {
        if (block.empty()) continue;
        const bool block_exist = isExistential(*(block.cbegin()));
        for (Variable var : block) {
            if (varDeleted(var)) {
                LOG(ERROR) << "Variable " << var << " is deleted but still contained in a quantifier block.";
                return false;
            }

            if (isExistential(var) != block_exist) {
                LOG(ERROR) << "Variable " << var << " is in a block with opposite quantifier.";
                return false;
            }
        }
    }

    return true;
}

bool
SATPrefix::checkConsistency() const
{
    if (!Prefix::checkConsistency()) return false;
    for (Variable var : _vars) {
        if (varDeleted(var)) return false;
        if (!isExistential(var)) return false;
    }
    return true;
}

//-----------------------------------------------------

/**
 * \brief Writes the DQBF prefix to an output stream.
 * \param stream the stream to write to
 * \param translation_table Used to rename the variables in case deleted
 * variables should be omitted
 */
void
DQBFPrefix::write(std::ostream& stream, std::vector<Variable>* translation_table) const
{
    auto trans = [translation_table](Variable var) -> Variable {
        return (translation_table ? (*translation_table)[var] : var);
    };

    // print universal variables
    stream << "a ";
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (isUniversal(var)) {
            stream << trans(var) << ' ';
        }
    }
    stream << "0\n";

    // print existential variables which do not depend on all universals
    for (Variable var = 1; var <= maxVarIndex(); ++var) {
        if (!isExistential(var)) continue;
        if (_in_rmb[var]) continue;

        stream << "d " << trans(var) << ' ';
        for (const Variable dep : _dependencies[var]) {
            stream << trans(dep) << ' ';
        }
        stream << "0\n";
    }

    // print existential variable which depend on all universals.
    if (!_rmb.empty()) {
        stream << "e ";
        for (const Variable var : _rmb) {
            stream << trans(var) << ' ';
        }
        stream << "0\n";
    }
}

/**
 * \brief Writes the QBF prefix to an output stream.
 * \param stream the stream to write to
 * \param translation_table Used to rename the variables in case deleted
 * variables should be omitted
 */
void
QBFPrefix::write(std::ostream& stream, std::vector<Variable>* translation_table) const
{
    auto trans = [translation_table](Variable var) -> Variable {
        return (translation_table ? (*translation_table)[var] : var);
    };

    for (const auto& block : _blocks) {
        if (block.empty()) continue;
        if (isExistential(*(block.begin())))
            stream << "e ";
        else
            stream << "a ";
        for (Variable var : block) stream << trans(var) << ' ';
        stream << "0\n";
    }
}

void
SATPrefix::write(std::ostream&, std::vector<Variable>*) const
{
    return;
}

//---------------------------------------------------------

void
Prefix::updateVars()
{}

void
SATPrefix::updateVars()
{
    return;
}

/**
 * \brief Update the variable data structures.
 *
 * Empty blocks are removed, adjacent blocks with the same
 * quantifier are merged. This changes the depth of the variables.
 */
void
QBFPrefix::updateVars()
{
    Prefix::updateVars();

    if (!_blocks.empty()) {
        std::size_t pos_done = 0;
        std::size_t pos_curr = 1;
        while (pos_curr < _blocks.size()) {
            if (_blocks[pos_curr].empty()) {
                ++pos_curr;
                continue;
            } else if (_blocks[pos_done].empty()) {
                _blocks[pos_done] = std::move(_blocks[pos_curr]);
                _blocks[pos_curr].clear();
                ++pos_curr;
            } else if (getLevelQuantifier(pos_done) != getLevelQuantifier(pos_curr)) {
                ++pos_done;
            } else if (getLevelQuantifier(pos_done) == getLevelQuantifier(pos_curr) && pos_curr != pos_done) {
                // Merge the two blocks
                _blocks[pos_done].insert(_blocks[pos_curr].cbegin(), _blocks[pos_curr].cend());
                _blocks[pos_curr].clear();
                ++pos_curr;
            } else {
                ++pos_curr;
            }
        }

        if (_blocks.size() > pos_done + 1) {
            _blocks.resize(pos_done + 1);
            for (std::size_t k = 0; k < _blocks.size(); ++k) {
                for (Variable var : _blocks[k]) {
                    _depth[var] = k;
                }
            }
            val_assert((!_blocks.back().empty()) || (numVars() == 0));
        }
    }
}

//---------------------------------------------------------

/**
 * \brief Returns the cumulated number of dependencies of all existential
 * variables.
 */
std::size_t
DQBFPrefix::numDependencies() const noexcept
{
    std::size_t num = _rmb.size() * numUVars();
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (isExistential(var)) num += numDependencies(var);
    }

    return num;
}

/**
 * \brief Returns the cumulated number of dependencies of all existential
 * variables.
 */
std::size_t
QBFPrefix::numDependencies() const noexcept
{
    std::size_t result = 0;
    std::size_t univs  = 0;
    for (const auto& block : _blocks) {
        if (block.empty()) continue;
        if (isExistential(*(block.begin())))
            result += univs * block.size();
        else
            univs += block.size();
    }

    return result;
}

std::size_t
SATPrefix::numDependencies() const noexcept
{
    return 0;
}

//---------------------------------------------------------

/**
 * \brief Returns the number of dependencies of a given existential variable.
 */
std::size_t
DQBFPrefix::numDependencies(const Variable var) const noexcept
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    val_assert(isExistential(var));

    if (_in_rmb[var])
        return _univ_vars.size();
    else
        return _dependencies[var].size();
}

/**
 * \brief Returns the number of dependencies of a given existential variable.
 */
std::size_t
QBFPrefix::numDependencies(const Variable var) const noexcept
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    if (varDeleted(var)) return 0;

    const std::size_t block_nr = getLevel(var);
    std::size_t       result   = 0;

    for (std::size_t b = 0; b < block_nr; ++b) {
        if (_blocks[b].empty()) continue;
        if (isUniversal(*(_blocks[b].begin()))) {
            result += _blocks[b].size();
        }
    }

    return result;
}

std::size_t
SATPrefix::numDependencies(const Variable) const noexcept
{
    return 0;
}

//---------------------------------------------------------

/**
 * \brief Removes all dependencies of a variable.
 */
void
DQBFPrefix::clearDependencies(Variable var)
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());

    if (_in_rmb[var]) {
        _in_rmb[var] = false;
        _rmb.erase(var);
        for (Variable univ : _univ_vars) {
            _dependencies[univ].erase(var);
        }
    } else {
        for (Variable dep : _dependencies[var]) {
            _dependencies[dep].erase(var);
        }
        _dependencies[var].clear();
    }
}

/**
 * Checks if a DQBF is actually a QBF.
 */
bool
DQBFPrefix::isQBF() const
{
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (!isExistential(var)) continue;
        if (_in_rmb[var]) continue;  // variable depends on all universals

        for (Variable other_var = var + 1; other_var <= maxVarIndex(); ++other_var) {
            if (!isExistential(other_var)) continue;
            if (_in_rmb[other_var]) continue;  // variable depends on all universals

            // remove from the first dep-set each element from the second dep-set
            // remove from the second dep-set each element from the first dep-set
            const auto result
                = two_sided_difference_empty(_dependencies[var].cbegin(), _dependencies[var].cend(),
                                             _dependencies[other_var].cbegin(), _dependencies[other_var].cend());
            // if both differences are non-empty, we have a proper DQBF
            if (!result.first && !result.second) return false;
        }
    }
    return true;
}

/**
 * \brief Converts the DQBF prefix into an equivalent QBF prefix (if possible).
 *
 * The QBF prefix is allocated on the heap. The caller is
 * responsible for deleting the object.

 * \pre The DQBF prefix must have an equivalent QBF prefix.
 * \sa DQBFPrefix::isQBF()
 */
QBFPrefix*
DQBFPrefix::convertToQBF() const
{
    if (!isQBF()) return nullptr;

    auto result = new QBFPrefix();
    result->setMaxVarIndex(maxVarIndex());

    // Get all existential variables
    std::vector<Variable> exist_vars;
    exist_vars.reserve(numEVars());
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (isExistential(var) && !_in_rmb[var]) exist_vars.push_back(var);
    }

    // Sort existential variables according to their number of dependencies.
    // This brings them into the right order for the QBF prefix.
    std::sort(exist_vars.begin(), exist_vars.end(),
              [this](Variable a, Variable b) { return _dependencies[a].size() < _dependencies[b].size(); });

    // Insert the universal variables into the order.
    std::vector<bool> finished_universal(maxVarIndex() + 1, false);
    for (const auto exist : exist_vars) {
        for (const auto dep : _dependencies[exist]) {
            if (!finished_universal[dep]) {
                finished_universal[dep] = true;
                result->addUVar(dep);
            }
        }
        result->addEVar(exist);
    }

    // Append the remaining universal variables which do not appear in
    // any dependency sets.
    for (Variable univ = minVarIndex(); univ <= maxVarIndex(); ++univ) {
        if (isUniversal(univ) && !finished_universal[univ]) {
            finished_universal[univ] = true;
            result->addUVar(univ);
        }
    }

    // Append the rmb-variables
    for (Variable var : _rmb) result->addEVar(var);

    return result;
}

void
QBFPrefix::setLevel(const Variable var, const std::size_t level)
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    val_assert(level <= _blocks.size());
    val_assert(!varDeleted(var));

    if (level == _depth[var]) return;

    val_assert(level == _blocks.size() || _blocks[level].empty()
               || isExistential(var) == isExistential(*(_blocks[level].begin())));
    val_assert(level < _blocks.size() || _blocks.back().empty()
               || isExistential(*(_blocks.back().begin())) != isExistential(var));

    _blocks[_depth[var]].erase(var);
    _depth[var] = level;
    if (level >= _blocks.size()) _blocks.resize(level + 1);
    _blocks[level].insert(var);
}

/**
 * \brief Turns a QBF prefix into an equivalent DQBF prefix.
 *
 * The DQBF prefix is allocated on the heap. The caller is
 * responsible for deleting the object.
 */
DQBFPrefix*
QBFPrefix::convertToDQBF() const
{
    VLOG(2) << "Switching to the DQBF prefix representation.";

    auto result = new DQBFPrefix();
    result->setMaxVarIndex(this->maxVarIndex());
    std::set<Variable> deps;
    for (const auto& block : _blocks) {
        for (Variable var : block) {
            if (varDeleted(var)) continue;
            if (isUniversal(var)) {
                result->addUVar(var);
                deps.insert(var);
            } else {
                result->addEVar(var);
                result->setDependencies(var, deps);
            }
        }
    }
    return result;
}

template<typename Container>
void
DQBFPrefix::setDependencies(Variable var, const Container& deps)
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    val_assert(isExistential(var));
    val_assert(_dependencies.size() > var);
#ifndef NDEBUG
    for (auto dep : deps) {
        val_assert(isUniversal(dep));
    }
#endif

    if (deps.size() == numUVars()) {
        moveToRMB(var);
        return;
    }

    if (_in_rmb[var]) {
        for (Variable dep : _univ_vars) _dependencies[dep].erase(var);
        _dependencies[var].clear();
        _in_rmb[var] = false;
        _rmb.erase(var);
    } else {
        // remove old dependencies
        for (Variable dep : _dependencies[var]) _dependencies[dep].erase(var);
        _dependencies[var].clear();
    }

    // add new dependencies
    auto pos = _dependencies[var].begin();
    for (Variable dep : deps) {
        _dependencies[dep].insert(var);
        pos = _dependencies[var].insert(pos, dep);
    }
}

void
DQBFPrefix::setDependencies(const Variable var, std::set<Variable>&& deps)
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    val_assert(isExistential(var));
#ifndef NDEBUG
    for (auto dep : deps) {
        val_assert(isUniversal(dep));
    }
#endif

    if (deps.size() == numUVars()) {
        moveToRMB(var);
        return;
    }

    if (_in_rmb[var]) {
        for (Variable dep : _univ_vars) _dependencies[dep].erase(var);
        _dependencies[var].clear();
        _in_rmb[var] = false;
        _rmb.erase(var);
    } else {
        // remove old dependencies
        for (Variable dep : _dependencies[var]) _dependencies[dep].erase(var);
        _dependencies[var].clear();
    }
    val_assert(checkConsistency());

    // add new dependencies
    for (Variable dep : deps) {
        _dependencies[dep].insert(var);
    }
    _dependencies[var] = std::move(deps);
}

/**
 * \brief Returns the set of variables in the right-most quantifier block.
 */
const std::set<Variable>&
QBFPrefix::getRMB() const noexcept
{
    int max_num = static_cast<int>(_blocks.size()) - 1;
    while (max_num >= 0 && _blocks[max_num].empty()) --max_num;
    if (max_num < 0 || getLevelQuantifier(max_num) != VariableStatus::EXISTENTIAL) {
        const_cast<QBFPrefix*>(this)->_blocks.resize(max_num + 2);
    } else {
        const_cast<QBFPrefix*>(this)->_blocks.resize(max_num + 1);
    }

    return _blocks.back();
}

/**
 * \brief Returns the set of existential variables in the right-most block.
 *
 * The right-most block consists of all those existential variables which
 * depend on <i>all</i> universal variables. These often play a special role
 * because they can be eliminated by resolution, all gate outputs after
 * Tseitin transformation can be moved to this set etc.
 */
const std::set<Variable>&
DQBFPrefix::getRMB() const noexcept
{
    return _rmb;
}

const std::set<Variable>&
SATPrefix::getRMB() const noexcept
{
    return _vars;
}

//-----------------------------------------------

/**
 * \brief Moves the given variable to the right-most block.
 *
 * Afterwards, the variable depends on all universal variables.
 * \pre The variable must be existential.
 */
void
QBFPrefix::moveToRMB(const Variable var)
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    val_assert(!varDeleted(var));
    val_assert(isExistential(var));

    if (getLevel(var) == _blocks.size() - 1) return;
    _blocks[_depth[var]].erase(var);

    int max_num = static_cast<int>(_blocks.size()) - 1;
    while (max_num >= 0 && _blocks[max_num].empty()) --max_num;
    if (max_num < 0 || getLevelQuantifier(max_num) != VariableStatus::EXISTENTIAL) {
        _blocks.resize(max_num + 2);
    } else {
        _blocks.resize(max_num + 1);
    }

    _blocks.back().insert(var);
    _depth[var] = _blocks.size() - 1;
}

/**
 * \brief Moves the given variable to the right-most block.
 *
 * Afterwards, the variable depends on all universal variables.
 * \pre The variable must be existential.
 */
void
DQBFPrefix::moveToRMB(const Variable var)
{
    val_assert(minVarIndex() <= var && var <= maxVarIndex());
    val_assert(!varDeleted(var));
    val_assert(isExistential(var));

    if (_dependencies[var].size() == numUVars()) return;

    _in_rmb[var] = true;
    _dependencies[var].clear();
    _rmb.insert(var);
    for (Variable univ : _univ_vars) _dependencies[univ].insert(var);
}

void
SATPrefix::moveToRMB(const Variable)
{
    return;
}

//-------------------------------------------------------------------------------
/**
 * \brief Returns the quantifier of a given QBF quantifier level.
 * \return VariableStatus::EXISTENTIAL or VariableStatus::UNIVERSAL
 * \pre the given level has to exist.
 */
VariableStatus
QBFPrefix::getLevelQuantifier(const std::size_t level) const noexcept
{
    if (level >= _blocks.size() || _blocks[level].empty())
        return VariableStatus::DELETED;
    else if (isExistential(*(_blocks[level].begin())))
        return VariableStatus::EXISTENTIAL;
    else
        return VariableStatus::UNIVERSAL;
}

SATPrefix*
Prefix::convertToSAT() const
{
    val_assert(this->isSAT());

    SATPrefix* result = new SATPrefix();
    result->setMaxVarIndex(this->maxVarIndex());

    for (Variable i = 0; i <= this->maxVarIndex(); ++i) {
        if (!varDeleted(i)) result->addEVar(i);
    }

    return result;
}

}  // end namespace hqspre
