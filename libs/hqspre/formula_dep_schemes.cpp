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
#include <set>
#include <stack>
#include <utility>
#include <vector>

#include <easylogging++.hpp>

#include "auxil.hpp"
#include "formula.hpp"
#include "gate.hpp"
#include "literal.hpp"
#include "prefix.hpp"
#include "timer.hpp"

/**
 * \file formula_dep_schemes.cpp
 * \brief Implementation of dependency schemes for DQBF
 * \author Ralf Wimmer
 * \date 12/2015--03/2016
 */

namespace hqspre {

/**
 * \brief Tries to find paths in the clause-literal incidence graph.
 *
 * \param start_lits the literals the paths may start at
 * \param forbidden a functor telling over which variable nodes the paths may
 * not run \param[out] seen a vector which stores the variables that can be
 * reached by paths. \sa stdTriDep(bool, bool) \sa invStdTriDep(bool, bool) \sa
 * sstdQuadDep(bool, bool) \sa invSstdQuadDep(bool, bool)
 */
template <typename Function>
void
Formula::searchPath(const std::vector<Literal>& start_lits, Function forbidden, std::vector<unsigned char>& seen) const
{
    val_assert(seen.size() > maxLitIndex());
    std::stack<Literal> open_lits;

    // Invariant:
    // - only literals are added to open_vars for which seen[lit] = false;
    // - all literals in open_vars are marked with seen[lit] = true
    // This ensures that no literal is added to open_vars twice.
    for (Literal lit : start_lits) {
        if (seen[lit] == 0) {
            open_lits.push(lit);
        }
        seen[lit] = 1u;
    }

    while (!open_lits.empty()) {
        const Literal current_lit = open_lits.top();
        open_lits.pop();
        for (ClauseID c_nr : _occ_list[current_lit]) {
            val_assert(!clauseDeleted(c_nr));
            if (clauseOptional(c_nr)) {
                continue;
            }  ///< \todo May we skip optional clauses?

            for (Literal next_lit : _clauses[c_nr]) {
                if (seen[next_lit] == 0u && !forbidden(next_lit)) {
                    open_lits.push(next_lit);
                }
                if (seen[negate(next_lit)] == 0u && !forbidden(negate(next_lit))) {
                    open_lits.push(negate(next_lit));
                }

                seen[next_lit]         = 1u;
                seen[negate(next_lit)] = 1u;
            }
        }
    }
}

/**
 * \brief Tries to find resolution paths in the clause-literal incidence graph.
 *
 * A resolution path connects clauses in which the connecting variable
 * appears in different polarities. Additionally, the same variable may not
 * be used twice in sequence along the path.
 * \param start_lits the literals the resolution paths may start at
 * \param forbidden a functor telling over which variable nodes the paths may
 * not run \param[out] seen a vector which stores the variables that can be
 * reached by paths. \sa stdTriDep(bool, bool) \sa invStdTriDep(bool, bool) \sa
 * sstdQuadDep(bool, bool) \sa invSstdQuadDep(bool, bool)
 */
template <typename Function>
void
Formula::searchResolutionPath(const std::vector<Literal>& start_lits, Function forbidden,
                              std::vector<unsigned char>& seen) const
{
    val_assert(seen.size() > maxLitIndex());
    std::stack<Literal> open_lits;

    // Invariant:
    // - only literals are added to open_vars for which seen[lit] = false;
    // - all literals in open_vars are marked with seen[lit] = true
    // This ensures that no literal is added to open_vars twice.
    for (const Literal lit : start_lits) {
        for (const ClauseID c_nr : _occ_list[lit]) {
            val_assert(!clauseDeleted(c_nr));
            if (clauseOptional(c_nr)) {
                continue;
            }  ///< \todo May we skip optional clauses?

            for (Literal next_lit : _clauses[c_nr]) {
                if (!forbidden(next_lit) && !seen[next_lit]) {
                    open_lits.push(next_lit);
                }
                seen[next_lit] = 1u;
            }
        }
    }

    while (!open_lits.empty()) {
        const Literal current_lit = negate(open_lits.top());
        open_lits.pop();
        for (const ClauseID c_nr : _occ_list[current_lit]) {
            val_assert(!clauseDeleted(c_nr));
            if (clauseOptional(c_nr)) {
                continue;
            }  ///< \todo May we skip optional clauses?

            for (Literal next_lit : _clauses[c_nr]) {
                if (seen[next_lit] == 0u && !forbidden(next_lit)) {
                    open_lits.push(next_lit);
                }
                seen[next_lit] = 1u;
            }
        }
    }
}

bool
Formula::invStdTriDep(bool resolution_paths, bool triangle)
{
    bool modified = false;
    for (Variable u_var = minVarIndex(); u_var <= maxVarIndex(); ++u_var) {
        if (!isUniversal(u_var)) {
            continue;
        }
        if (invStdTriDep(u_var, resolution_paths, triangle)) {
            modified = true;
        }
    }

    return modified;
}

/**
 * \brief Implementation of the standard and reflexive triangle dependency
 * scheme for addinging pseudo-dependencies
 *
 * The dependency scheme can either use normal
 * or resolution paths. Identified pseudo-dependencies are
 * added to the dependency sets.
 * \note For a QBF prefix, the function may only be called if the
 *       parameter `pseudo_deps` is not nullptr.
 * \param u_var the universal variable whose proper dependencies should be
 * computed \param resolution_paths if true, use resolution paths \param
 * triangle if true, use the reflexive triangle dependency scheme \param
 * pseudo_deps if not nullptr, the addable dependencies are stored in
 * pseudo_deps; otherwise they are added. \return true if the formula was
 * modified \sa stdTriDep(bool, bool) \sa sstdQuadDep(bool, bool) \sa
 * invSstdQuadDep(bool, bool)
 */
bool
Formula::invStdTriDep(Variable u_var, bool resolution_paths, bool triangle, std::set<Variable>* pseudo_deps)
{
    val_assert(minVarIndex() <= u_var && u_var <= maxVarIndex());
    val_assert(isUniversal(u_var));
    val_assert(_prefix);
    val_assert(_dqbf_prefix || pseudo_deps);

    std::size_t count = 0;

    const auto forbidden = [this, u_var](Literal lit) -> bool {
        return this->isUniversal(lit2var(lit)) || !this->depends(lit2var(lit), u_var);
    };

    if (resolution_paths) {
        searchResolutionPath({var2lit(u_var, true), var2lit(u_var, false)}, forbidden, _seen);
    } else {
        searchPath({var2lit(u_var, true), var2lit(u_var, false)}, forbidden, _seen);
    }

    for (Variable e_var = 1; e_var <= maxVarIndex(); ++e_var) {
        if (!isExistential(e_var) || depends(e_var, u_var)) {
            continue;
        }

        const Literal ep = var2lit(e_var, false);
        const Literal en = var2lit(e_var, true);

        if ((triangle && ((_seen[ep] == 0u) || (_seen[en] == 0u)))
            || (!triangle && (_seen[ep] == 0u) && (_seen[en] == 0u))) {
            if (pseudo_deps == nullptr) {
                addDependency(u_var, e_var);
            } else {
                pseudo_deps->insert(e_var);
            }
            ++count;
        }
    }

    clearSeen();
    stat(Statistics::ADD_DEPENDENCY_SCHEMES) += count;
    return count > 0;
}

bool
Formula::invSstdQuadDep(bool resolution_paths, bool quadrangle)
{
    bool modified = false;
    for (Variable u_var = minVarIndex(); u_var <= maxVarIndex(); ++u_var) {
        if (!isUniversal(u_var)) {
            continue;
        }
        if (invSstdQuadDep(u_var, resolution_paths, quadrangle)) {
            modified = true;
        }
    }

    return modified;
}

/**
 * \brief Implementation of the strict standard and reflexive quadrangle
 * dependency scheme for addinging pseudo-dependencies
 *
 * The dependency scheme can either use normal
 * or resolution paths. Identified pseudo-dependencies are
 * added to the dependency sets.
 * \note For a QBF prefix, the function may only be called if the
 *       parameter `pseudo_deps` is not nullptr.
 * \param u_var the universal variable whose proper dependencies should be
 * computed \param resolution_paths if true, use resolution paths \param
 * quadrangle if true, use the reflexive quadrangle dependency scheme \param
 * pseudo_deps if not nullptr, the addable dependencies are stored in
 * pseudo_deps; otherwise they are added. \return true if the formula was
 * modified \sa stdTriDep(bool, bool) \sa invStdTriDep(bool, bool) \sa
 * sstdQuadDep(bool, bool)
 */
bool
Formula::invSstdQuadDep(Variable u_var, bool resolution_paths, bool quadrangle, std::set<Variable>* pseudo_deps)
{
    val_assert(_prefix);
    val_assert(_dqbf_prefix || pseudo_deps);

    std::size_t                 count    = 0;
    std::vector<unsigned char>& seen_pos = _seen;
    std::vector<unsigned char>  seen_neg(maxLitIndex() + 1, 0u);

    const auto forbidden = [this, u_var](Literal lit) -> bool {
        return this->isUniversal(lit2var(lit)) || !this->depends(lit2var(lit), u_var);
    };

    if (resolution_paths) {
        searchResolutionPath({var2lit(u_var, false)}, forbidden, seen_pos);
        searchResolutionPath({var2lit(u_var, true)}, forbidden, seen_neg);
    } else {
        searchPath({var2lit(u_var, false)}, forbidden, seen_pos);
        searchPath({var2lit(u_var, true)}, forbidden, seen_neg);
    }

    for (Variable e_var = minVarIndex(); e_var <= maxVarIndex(); ++e_var) {
        if (!isExistential(e_var) || depends(e_var, u_var)) {
            continue;
        }

        const Literal ep = var2lit(e_var, false);
        const Literal en = var2lit(e_var, true);

        if (quadrangle) {
            if (((seen_pos[ep] == 0u) || (seen_neg[ep] == 0u)) && ((seen_neg[ep] == 0u) || (seen_pos[en] == 0u))) {
                if (pseudo_deps == nullptr) {
                    addDependency(u_var, e_var);
                } else {
                    pseudo_deps->insert(e_var);
                }
                ++count;
            }
        } else {
            if (((seen_pos[ep] == 0u) && (seen_pos[en] == 0u)) || ((seen_neg[ep] == 0u) && (seen_neg[en] == 0u))) {
                if (pseudo_deps == nullptr) {
                    addDependency(u_var, e_var);
                } else {
                    pseudo_deps->insert(e_var);
                }
                ++count;
            }
        }
    }

    clearSeen();
    stat(Statistics::ADD_DEPENDENCY_SCHEMES) += count;
    return count > 0;
}

bool
Formula::stdTriDep(bool resolution_paths, bool triangle)
{
    bool modified = false;
    for (Variable u_var = minVarIndex(); u_var <= maxVarIndex(); ++u_var) {
        if (!isUniversal(u_var)) {
            continue;
        }
        if (stdTriDep(u_var, resolution_paths, triangle)) {
            modified = true;
        }
    }

    return modified;
}

/**
 * \brief Implementation of the standard and reflexive triangle dependency
 * scheme for removing pseudo-dependencies
 *
 * The dependency scheme can either use normal
 * or resolution paths. Identified pseudo-dependencies are
 * removed from the dependency sets.
 * \note For a QBF prefix, the function may only be called if the
 *       parameter `pseudo_deps` is not nullptr.
 * \param u_var the universal variables whose dependencies are to be analyzed
 * \param resolution_paths if true, use resolution paths
 * \param triangle if true use the reflexive triangle instead of the standard
 * dependency scheme \param pseudo_deps if not nullptr, the removable
 * dependencies are stored in pseudo_deps; otherwise they are removed. \return
 * true if the formula was modified \sa invStdTriDep(bool, bool) \sa
 * sstdQuadDep(bool, bool) \sa invSstdQuadDep(bool, bool)
 */
bool
Formula::stdTriDep(Variable u_var, bool resolution_paths, bool triangle, std::set<Variable>* pseudo_deps)
{
    val_assert(minVarIndex() <= u_var && u_var <= maxVarIndex());
    val_assert(isUniversal(u_var));
    val_assert(_prefix);
    val_assert(_dqbf_prefix || pseudo_deps);

    std::size_t count = 0;
    const auto  forbidden
        = [this, u_var](Literal lit) -> bool { return isUniversal(lit2var(lit)) || !depends(lit2var(lit), u_var); };

    // Peforming depth-first search starting at 'u_var'.
    if (resolution_paths) {
        searchResolutionPath({var2lit(u_var, true), var2lit(u_var, false)}, forbidden, _seen);
    } else {
        searchPath({var2lit(u_var, true), var2lit(u_var, false)}, forbidden, _seen);
    }

    for (Variable e_var = 1; e_var <= maxVarIndex(); ++e_var) {
        if (!isExistential(e_var) || !depends(e_var, u_var)) {
            continue;
        }

        const Literal ep = var2lit(e_var, false);
        const Literal en = var2lit(e_var, true);

        if ((triangle && ((_seen[en] == 0u) || (_seen[ep] == 0u)))
            || (!triangle && (_seen[en] == 0u) && (_seen[ep] == 0u))) {
            if (pseudo_deps == nullptr) {
                removeDependency(e_var, u_var);
            } else {
                pseudo_deps->insert(e_var);
            }
            ++count;
        }
    }

    clearSeen();
    stat(Statistics::REM_DEPENDENCY_SCHEMES) += count;
    return count > 0;
}

bool
Formula::sstdQuadDep(bool resolution_paths, bool quadrangle)
{
    bool modified = false;
    for (Variable u_var = minVarIndex(); u_var <= maxVarIndex(); ++u_var) {
        if (!isUniversal(u_var)) {
            continue;
        }
        if (sstdQuadDep(u_var, resolution_paths, quadrangle)) {
            modified = true;
        }
    }

    return modified;
}

/**
 * \brief Implementation of the strict standard and reflexive quadrangle
 * dependency scheme for removing pseudo-dependencies
 *
 * The dependency scheme can either use normal
 * or resolution paths. Identified pseudo-dependencies are
 * removed from the dependency sets.
 * \note For a QBF prefix, the function may only be called if the
 *       parameter `pseudo_deps` is not nullptr.
 * \param u_var the universal variable whose proper dependencies should be
 * computed \param resolution_paths if true, use resolution paths \param
 * quadrangle if true, reflexive quadrangle dependencies are computed \param
 * pseudo_deps if not nullptr, the removable dependencies are stored in
 * pseudo_deps; otherwise they are removed. \return true if the formula was
 * modified \sa stdTriDep(bool, bool) \sa invStdTriDep(bool, bool) \sa
 * invSstdQuadDep(bool, bool)
 */
bool
Formula::sstdQuadDep(Variable u_var, bool resolution_paths, bool quadrangle, std::set<Variable>* pseudo_deps)
{
    val_assert(minVarIndex() <= u_var && u_var <= maxVarIndex());
    val_assert(isUniversal(u_var));
    val_assert(_prefix);
    val_assert(_dqbf_prefix || pseudo_deps);

    std::size_t                 count    = 0;
    std::vector<unsigned char>& seen_pos = _seen;
    std::vector<unsigned char>  seen_neg(maxLitIndex() + 1, 0u);

    const auto forbidden
        = [this, u_var](Literal lit) -> bool { return isUniversal(lit2var(lit)) || !depends(lit2var(lit), u_var); };

    seen_pos.assign(seen_pos.size(), 0u);
    seen_neg.assign(seen_neg.size(), 0u);

    if (resolution_paths) {
        // Peforming depth-first search starting at 'u_var'.
        searchResolutionPath({var2lit(u_var, false)}, forbidden, seen_pos);

        // Peforming depth-first search starting at '!u_var'.
        searchResolutionPath({var2lit(u_var, true)}, forbidden, seen_neg);
    } else {
        // Peforming depth-first search starting at 'u_var'.
        searchPath({var2lit(u_var, false)}, forbidden, seen_pos);

        // Peforming depth-first search starting at '!u_var'.
        searchPath({var2lit(u_var, true)}, forbidden, seen_neg);
    }

    for (Variable e_var = 1; e_var <= maxVarIndex(); ++e_var) {
        if (!isExistential(e_var) || !depends(e_var, u_var)) {
            continue;
        }

        const Literal ep = var2lit(e_var, false);
        const Literal en = var2lit(e_var, true);

        if (quadrangle) {
            if (((seen_pos[ep] == 0u) || (seen_neg[ep] == 0u)) || ((seen_neg[ep] == 0u) && (seen_pos[en] == 0u))) {
                if (pseudo_deps == nullptr) {
                    removeDependency(e_var, u_var);
                } else {
                    pseudo_deps->insert(e_var);
                }
                ++count;
            }
        } else {
            if (((seen_pos[en] == 0u) && (seen_pos[ep] == 0u)) || ((seen_neg[en] == 0u) && (seen_neg[ep] == 0u))) {
                if (pseudo_deps == nullptr) {
                    removeDependency(e_var, u_var);
                } else {
                    pseudo_deps->insert(e_var);
                }
                ++count;
            }
        }
    }

    clearSeen();
    stat(Statistics::REM_DEPENDENCY_SCHEMES) += count;
    return count > 0;
}

/**
 * \brief Uses gate detection to identify pseudo-dependencies.
 * \param operation specifies whether pseudo-dependencies are added, removed, or
 * only identified. \return true if the formula was modified.
 */
bool
Formula::gateDependencies(const DependencyOperation operation)
{
    val_assert(checkConsistency());
    val_assert(_prefix);
    val_assert(operation != DependencyOperation::REMOVE || _prefix->type() == PrefixType::DQBF);

    std::size_t count = 0;

    _gates.update();

    if (operation == DependencyOperation::ADD) {
        // We can make all gate outputs depend on all universal variables
        // as they are implied by the gate inputs
        for (const auto& g : _gates) {
            if (!_gates.gateValid(g)) {
                continue;
            }

            const Variable output_var = lit2var(g._output_literal);
            count += (_prefix->numUVars() - numDependencies(output_var));
            _prefix->moveToRMB(output_var);
        }
        stat(Statistics::ADD_DEPENDENCY_SCHEMES) += count;
    } else if (operation == DependencyOperation::REMOVE) {
        for (const auto& g : _gates) {
            if (!_gates.gateValid(g)) {
                continue;
            }

            const Variable     outp_var = lit2var(g._output_literal);
            std::set<Variable> new_deps;
            for (const Literal inp_lit : g._input_literals) {
                const Variable inp_var = lit2var(inp_lit);
                if (isUniversal(inp_var)) {
                    new_deps.insert(inp_var);
                } else {
                    const auto& inp_deps = _dqbf_prefix->getDependencies(inp_var);
                    new_deps.insert(inp_deps.cbegin(), inp_deps.cend());
                }
            }

            count += _dqbf_prefix->numDependencies(outp_var) - new_deps.size();
            if (new_deps.size() == numUVars()) {
                _dqbf_prefix->moveToRMB(outp_var);
            } else {
                _dqbf_prefix->setDependencies(outp_var, std::move(new_deps));
            }
        }
        stat(Statistics::REM_DEPENDENCY_SCHEMES) += count;
    }

    return count != 0;
}

/**
 * \brief Applies a dependency scheme to the formula to manipulate the
 * dependency sets. \param scheme specifies which dependency scheme is applied
 * \param operation determines whether dependencies should be added or removed
 * \return true if the formula was modified.
 */
bool
Formula::applyDependencyScheme(DependencyScheme scheme, DependencyOperation operation)
{
    val_assert(_prefix);
    val_assert_msg(_prefix->type() == PrefixType::DQBF,
                   "Dependency schemes are currently only supported for DQBF prefixes");
    val_assert(_dqbf_prefix);

    ScopeTimer dependency(getTimer(WhichTimer::DEPENDENCY_SCHEMES));

#ifdef SKOLEM
    if (operation == DependencyOperation::ADD && scheme != DependencyScheme::TRIVIAL
        && scheme != DependencyScheme::GATE) {
        LOG(WARNING) << "Adding dependencies using dependency schemes is not "
                        "supported when Skolem functions are computed.\n";
        return false;
    }
#endif

    switch (scheme) {
        case DependencyScheme::TRIVIAL:
            return false;
        case DependencyScheme::STANDARD:
            if (operation == DependencyOperation::REMOVE) {
                return stdTriDep(false, false);
            } else if (operation == DependencyOperation::ADD) {
                return invStdTriDep(false, false);
            }
            break;
        case DependencyScheme::STRICT_STANDARD:
            if (operation == DependencyOperation::REMOVE) {
                return sstdQuadDep(false, false);
            } else if (operation == DependencyOperation::ADD) {
                return invSstdQuadDep(false, false);
            }
            break;
        case DependencyScheme::REF_TRIANGLE:
            if (operation == DependencyOperation::REMOVE) {
                return stdTriDep(false, true);
            } else if (operation == DependencyOperation::ADD) {
                return invStdTriDep(false, true);
            }
            break;
        case DependencyScheme::REF_QUADRANGLE:
            if (operation == DependencyOperation::REMOVE) {
                return sstdQuadDep(false, true);
            } else if (operation == DependencyOperation::ADD) {
                return invSstdQuadDep(false, true);
            }
            break;
        case DependencyScheme::RP_STANDARD:
            if (operation == DependencyOperation::REMOVE) {
                return stdTriDep(true, false);
            } else if (operation == DependencyOperation::ADD) {
                return invStdTriDep(true, false);
            }
            break;
        case DependencyScheme::RP_STRICT_STANDARD:
            if (operation == DependencyOperation::REMOVE) {
                return sstdQuadDep(true, false);
            } else if (operation == DependencyOperation::ADD) {
                return invSstdQuadDep(true, false);
            }
            break;
        case DependencyScheme::RP_REF_TRIANGLE:
            if (operation == DependencyOperation::REMOVE) {
                return stdTriDep(true, true);
            } else if (operation == DependencyOperation::ADD) {
                return invStdTriDep(true, true);
            }
            break;
        case DependencyScheme::RP_REF_QUADRANGLE:
            if (operation == DependencyOperation::REMOVE) {
                return sstdQuadDep(true, true);
            } else if (operation == DependencyOperation::ADD) {
                return invSstdQuadDep(true, true);
            }
            break;
        case DependencyScheme::GATE:
            return gateDependencies(operation);
            break;
        default:
            LOG(ERROR) << "Invalid dependency scheme.";
            return false;
            break;
    }  // end switch

    return false;
}

std::vector<std::vector<Variable>>
Formula::identifyDontCares()
{
    VLOG(1) << __FUNCTION__;

    if (_prefix->type() != PrefixType::DQBF) {
        return std::vector<std::vector<Variable>>();
    }

    val_assert(_dqbf_prefix && !_qbf_prefix);
    std::vector<std::vector<Literal>> result(maxVarIndex() + 1);

    // First maximize the dependencies

#ifndef SKOLEM
    // Adding dependencies is only possible when we do not compute Skolem
    // functions.
    applyDependencyScheme(DependencyScheme::RP_REF_QUADRANGLE, DependencyOperation::ADD);
#endif
    applyDependencyScheme(DependencyScheme::GATE, DependencyOperation::ADD);

    DQBFPrefix max_prefix = (*_dqbf_prefix);

    applyDependencyScheme(DependencyScheme::GATE, DependencyOperation::REMOVE);
    applyDependencyScheme(DependencyScheme::RP_REF_QUADRANGLE, DependencyOperation::REMOVE);

    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (!isExistential(var)) {
            continue;
        }
        if (max_prefix.inRMB(var)) {
            _prefix->moveToRMB(var);
            continue;
        }

        const auto& max_dep = max_prefix.getDependencies(var);
        const auto& min_dep = this->getDependencies(var);
        std::set_difference(max_dep.cbegin(), max_dep.cend(), min_dep.cbegin(), min_dep.cend(),
                            std::back_inserter(result[var]));
    }

    std::size_t num_dc = 0;
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (isExistential(var)) {
            if (max_prefix.inRMB(var)) {
                _prefix->moveToRMB(var);
                result[var].clear();
            } else {
                num_dc += result[var].size();
            }
        }
    }

    VLOG(2) << __FUNCTION__ << " found " << num_dc << " don't-care dependencies.";

    return result;
}

}  // end namespace hqspre
