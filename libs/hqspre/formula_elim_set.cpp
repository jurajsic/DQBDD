// $Id: formula_elim_set.cpp 2644 2019-09-07 20:46:54Z wimmer $

/*
 * This file is part of HQSpre.
 *
 * Copyright 2016-18 Ralf Wimmer, Sven Reimer, Paolo Marin, Bernd Becker
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

#define ELPP_STL_LOGGING

#include "formula.hpp"

#include <cmath>
#include <deque>
#include <functional>
#include <map>
#include <unordered_set>
#include <vector>

#include <easylogging++.hpp>

#include <glpk.h>

#include "auxil.hpp"
#include "literal.hpp"
#include "milp_solver.hpp"
#include "prefix.hpp"
#include "timer.hpp"

namespace std {
template <typename T1, typename T2>
class hash<std::pair<T1, T2>>
{
   public:
    std::size_t operator()(const std::pair<T1, T2>& s) const
    {
        const std::size_t h1 = std::hash<T1>()(s.first);
        const std::size_t h2 = std::hash<T2>()(s.second);
        return h1 ^ (h2 << 1u);
    }
};
}  // end namespace std

namespace hqspre {

namespace {

struct GlpkData
{
    GlpkData(const std::vector<std::vector<int>>& graph, const std::vector<bool>& classuniv,
             const std::vector<MilpSolver::VarType>& fy, const std::vector<MilpSolver::VarType>& zy,
             const std::map<std::pair<int, int>, MilpSolver::VarType>& dxy) noexcept :
        _graph(graph),
        _classuniv(classuniv),
        _fy(fy),
        _zy(zy),
        _dxy(dxy),
        _num_secant_cuts(0),
        _num_loop_cuts(0)
    {}

    const std::vector<std::vector<int>>&                      _graph;
    const std::vector<bool>&                                  _classuniv;
    const std::vector<MilpSolver::VarType>&                   _fy;
    const std::vector<MilpSolver::VarType>&                   _zy;
    const std::map<std::pair<int, int>, MilpSolver::VarType>& _dxy;
    unsigned int                                              _num_secant_cuts;
    unsigned int                                              _num_loop_cuts;
};

void
glpk_callback(glp_tree* T, void* info)
{
    // We need the callback only when requested to add new rows.
    if (glp_ios_reason(T) != GLP_IROWGEN) {
        return;
    }

    auto* data   = static_cast<GlpkData*>(info);
    auto  solver = glp_ios_get_prob(T);

    // Only add a cut if an integer solution has been found.
    for (int y = 0; y < static_cast<int>(data->_graph.size()); ++y) {
        if (data->_classuniv[y]) {
            continue;
        }
        const double actual_val_zy = glp_get_col_prim(solver, data->_zy[y]);
        if (std::floor(actual_val_zy) != actual_val_zy) {
            return;
        }
        const double actual_val_fy = glp_get_col_prim(solver, data->_fy[y]);
        if (std::floor(actual_val_fy) != actual_val_fy) {
            return;
        }
    }

    for (const auto& d : data->_dxy) {
        const double val = glp_get_col_prim(solver, d.second);
        if (std::floor(val) != val) {
            return;
        }
    }

    // Now check if we can add lazy constraints.
    std::vector<MilpSolver::VarType> vars;
    std::vector<double>              coeffs;
    bool                             added = false;

    // First consider the exponential constraint zy = 2^{z_y}
    for (int y = 0; y < static_cast<int>(data->_graph.size()); ++y) {
        if (data->_classuniv[y]) {
            continue;
        }

        const double actual_val_zy = glp_get_col_prim(solver, data->_zy[y]);
        const double actual_val_fy = glp_get_col_prim(solver, data->_fy[y]);
        const double wanted_val_zy = std::pow(2.0, actual_val_fy);

        if (actual_val_zy < wanted_val_zy) {
            vars.clear();
            coeffs.clear();

            vars.push_back(0);
            coeffs.push_back(0.0);  // Glpk ignores position 0 in the vectors!
            vars.push_back(data->_fy[y]);
            coeffs.push_back(wanted_val_zy);
            vars.push_back(data->_zy[y]);
            coeffs.push_back(-1.0);
            const auto rownr = glp_add_rows(solver, 2);

            glp_set_mat_row(solver, rownr, static_cast<int>(vars.size()) - 1, &vars[0], &coeffs[0]);
            glp_set_row_bnds(solver, rownr, GLP_UP, 0.0, wanted_val_zy * (actual_val_fy - 1.0));

            coeffs[1] *= 0.5;
            glp_set_mat_row(solver, rownr + 1, static_cast<int>(vars.size()) - 1, &vars[0], &coeffs[0]);
            glp_set_row_bnds(solver, rownr + 1, GLP_UP, 0.0, 0.5 * wanted_val_zy * (actual_val_fy - 2.0));

            data->_num_secant_cuts += 2;
            added = true;
        }
    }

    if (added) {
        return;
    }

    // Now check if the graph has become acyclic
    const std::function<bool(int, int)> edge_available = [solver, data](int x, int y) -> bool {
        const auto var = data->_dxy.find(std::make_pair(x, y));
        if (var == data->_dxy.cend()) {
            return true;
        } else {
            return glp_get_col_prim(solver, var->second) < 0.1;
        }
    };

    std::vector<bool> seen(data->_graph.size(), false);
    std::vector<int>  onPath(data->_graph.size(), -1);
    std::vector<int>  path;

    // Performs a depth-first search starting from the given node to find a cycle
    // in the graph. The function returns true iff a cycle has been found. In this
    // case, 'path' stores the lasso-shaped path, ending with the first node on
    // the cycle. The vector onPath contains the positions of the nodes on the
    // path, where onPath[path.back()] is the position of the first cycle node on
    // the lasso.
    const std::function<bool(int)> findCycle
        = [data, &seen, &onPath, &path, &edge_available, &findCycle](int node) -> bool {
        if (seen[node]) {
            return false;
        }

        if (onPath[node] > 0) {
            path.push_back(node);
            return true;
        }
        onPath[node] = static_cast<int>(path.size());
        path.push_back(node);

        for (int succ : data->_graph[node]) {
            if (!data->_classuniv[node] || edge_available(node, succ)) {
                bool result = findCycle(succ);
                if (result) {
                    return result;
                }
            }
        }
        seen[node]          = true;
        onPath[path.back()] = -1;
        path.pop_back();

        return false;
    };

    // Now perform a depth-first search from every node. If a cycle is found, add
    // the corresponding cycle exclusion constraint.
    for (int x = 0; x < static_cast<int>(data->_graph.size()); ++x) {
        if (!findCycle(x)) {
            continue;
        }

        vars.resize(1);
        coeffs.resize(1);

        for (int node_pos = onPath[path.back()]; node_pos < static_cast<int>(path.size()); ++node_pos) {
            const int x2 = path[node_pos];
            if (data->_classuniv[x2] && node_pos < static_cast<int>(path.size()) - 1) {
                const int  y     = path[node_pos + 1];
                const auto found = data->_dxy.find(std::make_pair(x2, y));
                vars.push_back(found->second);
                coeffs.push_back(1.0);
            }
        }

        const auto rownr = glp_add_rows(solver, 1);
        glp_set_mat_row(solver, rownr, static_cast<int>(vars.size()) - 1, &vars[0], &coeffs[0]);
        glp_set_row_bnds(solver, rownr, GLP_LO, 1.0, 0.0);
        ++data->_num_loop_cuts;
        return;
    }
}

}  // end anonymous namespace

static std::unordered_set<std::pair<int, int>>
solveMilp(const std::vector<std::vector<int>>& graph, const std::vector<std::vector<Variable>>& class2vars,
          const std::vector<bool>& classuniv)
{
    std::unique_ptr<MilpSolver> solver = nullptr;

    solver = std::make_unique<GlpkSolver>();

    val_assert(solver);

    solver->setVerbosity(false);

    // Vectors for storing linear expressions
    std::vector<MilpSolver::VarType> vars;
    vars.reserve(class2vars.size());
    std::vector<double> coeffs;
    coeffs.reserve(class2vars.size());

    // Create all variables
    std::vector<MilpSolver::VarType>                   zy(class2vars.size());
    std::vector<MilpSolver::VarType>                   fy(class2vars.size());
    std::map<std::pair<int, int>, MilpSolver::VarType> dxy;

    for (int y = 0; y < static_cast<int>(class2vars.size()); ++y) {
        if (classuniv[y]) {
            continue;  // skip universal variables
        }

        zy[y] = solver->addVariable(MilpSolver::VarSort::INT, true, 1.0, false, 0.0);
        fy[y] = solver->addVariable(MilpSolver::VarSort::INT, true, 0.0, false, 0.0);

        // create the objective function
        vars.push_back(zy[y]);
        coeffs.push_back(static_cast<double>(class2vars[y].size()));
    }

    for (int x = 0; x < static_cast<int>(graph.size()); ++x) {
        if (!classuniv[x]) {
            continue;
        }
        for (const int y : graph[x]) {
            dxy[std::make_pair(x, y)] = solver->addVariable(MilpSolver::VarSort::INT, true, 0.0, true, 1.0);
        }
    }

    // Set the objective function
    solver->setObjectiveDirection(MilpSolver::ObjectiveType::MINIMIZE);
    solver->setObjective(vars, coeffs);

    // For defining the f_y constraints, we need the reversed graph
    std::vector<std::vector<int>> reverse_graph(graph.size());
    for (int u = 0; u < static_cast<int>(graph.size()); ++u) {
        for (const int v : graph[u]) {
            reverse_graph[v].push_back(u);
        }
    }

    for (int y = 0; y < static_cast<int>(class2vars.size()); ++y) {
        if (classuniv[y]) {
            continue;  // skip universal variables
        }
        vars.clear();
        coeffs.clear();

        vars.push_back(fy[y]);
        coeffs.push_back(-1.0);
        for (const int x : reverse_graph[y]) {
            vars.push_back(dxy[std::make_pair(x, y)]);
            coeffs.push_back(static_cast<double>(class2vars[x].size()));
        }

        solver->addConstraint(vars, coeffs, true, 0.0, true, 0.0);
    }

    // Find all cycles of length 4
    vars.resize(2);
    coeffs.resize(2);
    for (std::size_t x1 = 0; x1 < graph.size(); ++x1) {
        if (!classuniv[x1]) {
            continue;
        }
        for (std::size_t y1 : graph[x1]) {
            for (std::size_t x2 : graph[y1]) {
                if (x2 <= x1) {
                    continue;
                }
                for (std::size_t y2 : graph[x2]) {
                    if (std::find(graph[y2].cbegin(), graph[y2].cend(), x1) != graph[y2].cend()) {
                        // 4-cycle x1 -> y1 -> x2 -> y2 -> x1 found.
                        VLOG(3) << "cycle: " << x1 << " --> " << y1 << " --> " << x2 << " --> " << y2 << " --> " << x1;
                        val_assert(dxy.find(std::make_pair(x1, y1)) != dxy.cend());
                        vars[0]   = dxy[std::make_pair(x1, y1)];
                        coeffs[0] = 1.0;

                        val_assert(dxy.find(std::make_pair(x2, y2)) != dxy.cend());
                        vars[1]   = dxy[std::make_pair(x2, y2)];
                        coeffs[1] = 1.0;
                        solver->addConstraint(vars, coeffs, true, 1.0, true, 2.0);
                    }
                }
            }
        }
    }

    TruthValue result = TruthValue::UNKNOWN;

    GlpkData callback_data(graph, classuniv, fy, zy, dxy);
    result = dynamic_cast<GlpkSolver*>(solver.get())->solve(glpk_callback, &callback_data);
    VLOG(2) << "Added " << callback_data._num_secant_cuts << " secant constraints.";
    VLOG(2) << "Added " << callback_data._num_loop_cuts << " loop exclusion constraints.";

    std::unordered_set<std::pair<int, int>> to_eliminate;

    if (result == TruthValue::TRUE) {
        for (const auto d : dxy) {
            if (solver->getIntValue(d.second) == 1) {
                to_eliminate.insert(std::make_pair(d.first.first, d.first.second));
            }
        }
    } else {
        throw ElimSetException("MILP-solver did not find a solution.");
    }

    return to_eliminate;
}

/**
 * \brief Kahn's algorithm for topogical sorting
 *
 * For details see:<br>
 * Arthur B. Kahn: <i>Topological sorting of large networks</i>,
 * Communications of the ACM, 5 (11): 558–562, 1962.
 * DOI: 10.1145/368996.369025
 */
static std::vector<unsigned int>
topoSort(const std::vector<std::vector<int>>& graph, const std::function<bool(int, int)>& edge_forbidden)
{
    std::vector<unsigned int> indegree(graph.size(), 0);

    for (int x = 0; x < static_cast<int>(graph.size()); ++x) {
        for (const int y : graph[x]) {
            if (!edge_forbidden(x, y)) {
                ++indegree[y];
            }
        }
    }

    std::deque<int>           queue;
    std::vector<unsigned int> result(graph.size(), 0);
    unsigned int              current_number = 0;

    for (int x = 0; x < static_cast<int>(graph.size()); ++x) {
        if (indegree[x] == 0) {
            queue.push_back(x);
        }
    }

    while (!queue.empty()) {
        const unsigned int current_node = queue.front();
        queue.pop_front();
        result[current_node] = current_number;
        ++current_number;

        for (int succ : graph[current_node]) {
            if (!edge_forbidden(current_node, succ)) {
                --indegree[succ];
                if (indegree[succ] == 0) {
                    queue.push_back(succ);
                }
            }
        }
    }

    return result;
}

/**
 * \brief Computes a cost-minimal set of dependencies whose elimination turns
 * the DQBF into a QBF.
 */
std::vector<std::vector<Variable>>
Formula::computeDepElimSet()
{
    ScopeTimer et(getTimer(WhichTimer::DEP_ELIM_SET));

    val_assert(_dqbf_prefix != nullptr);

    int                                class_index = 0;
    std::map<std::set<Variable>, int>  dep2class;
    std::vector<int>                   var2class(maxVarIndex() + 1, -1);
    std::vector<std::vector<Variable>> class2vars;
    std::vector<bool>                  classuniv;

    // Perform symmetry reduction - all variables with the same dependencies
    // belong to the same class.
    for (Variable var = minVarIndex(); var <= maxVarIndex(); ++var) {
        if (varDeleted(var) || _prefix->inRMB(var)) {
            continue;
        }
        const auto& deps = _dqbf_prefix->getDependencies(var);

        const auto found = dep2class.find(deps);
        if (found == dep2class.end()) {
            dep2class[deps] = class_index;
            var2class[var]  = class_index;
            class2vars.emplace_back(1, var);
            classuniv.push_back(isUniversal(var));
            ++class_index;
        } else {
            const int index = found->second;
            var2class[var]  = index;
            class2vars[index].push_back(var);
            val_assert(isUniversal(var) == classuniv[index]);
        }
    }

    // Create the symmetry-reduced dependency graph
    std::vector<std::vector<int>> graph(class_index);
    for (int x = 0; x < class_index; ++x) {
        const Variable var_x = class2vars[x].front();
        for (int y = x; y < class_index; ++y) {
            const Variable var_y = class2vars[y].front();

            if (isExistential(var_x) == isExistential(var_y)) {
                continue;
            } else if (isExistential(var_y)) {
                // x universal, y existential
                if (_prefix->depends(var_y, var_x)) {
                    graph[x].push_back(y);
                } else {
                    graph[y].push_back(x);
                }
            } else {
                // x existential, y universal
                if (_prefix->depends(var_x, var_y)) {
                    graph[y].push_back(x);
                } else {
                    graph[x].push_back(y);
                }
            }
        }
    }

    // We determine an elimination set of the symmetry-reduced graph by solving an
    // MILP
    const auto to_elim = solveMilp(graph, class2vars, classuniv);

    // Compute a topological ordering of the dependency graph without the selected
    // edges. Insert the selected edges such that they point into the right
    // direction w.r.t. the topological ordering.
    const std::vector<unsigned int> order = topoSort(
        graph, [&to_elim](int x, int y) -> bool { return to_elim.find(std::make_pair(x, y)) != to_elim.cend(); });

    // Now undo the symmetry reduction and return the elimination set in terms of
    // DQBF variables. Only those edges of the graph are relevant which point into
    // the "wrong" direction w.r.t. the topological order.
    std::vector<std::vector<Variable>> result(maxVarIndex() + 1);
    for (const auto& edge : to_elim) {
        if (order[edge.first] > order[edge.second]) {
            for (const Variable exist : class2vars[edge.second]) {
                result[exist].insert(result[exist].end(), class2vars[edge.first].cbegin(),
                                     class2vars[edge.first].cend());
            }
        }
    }

    return result;
}

}  // end namespace hqspre
