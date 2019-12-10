// $Id: milp_solver.cpp 1965 2018-01-27 16:11:37Z wimmer $

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


#include <cmath>
#include <cstdlib>
#include <vector>

#if defined(HAVE_GLPK)
#include <glpk.h>
#endif

#include <easylogging++.hpp>
#include "milp_solver.hpp"


namespace hqspre {


#if defined(HAVE_GLPK)

GlpkSolver::GlpkSolver():
    _solver(nullptr),
    _max_constraint(1),
    _curr_constraint(1)
{
    _solver = glp_create_prob();
    setObjectiveDirection(ObjectiveType::MINIMIZE);
}


GlpkSolver::~GlpkSolver()
{
    glp_delete_prob(_solver);
}


GlpkSolver::VarType GlpkSolver::addVariable(MilpSolver::VarSort sort,
                                            bool lower_bounded,
                                            double lower_bound,
                                            bool upper_bounded,
                                            double upper_bound)
{
    // Create a new variable
    GlpkSolver::VarType result = glp_add_cols(_solver, 1);

    // Set the variable's type
    if (sort == MilpSolver::VarSort::INT) {
        glp_set_col_kind(_solver, result, GLP_IV);
    } else if (sort == MilpSolver::VarSort::REAL) {
        glp_set_col_kind(_solver, result, GLP_CV);
    }

    // Set the variable's bounds
    if (lower_bounded && upper_bounded) {
        // bounded on both sides
        if (lower_bound == upper_bound) {
            // constant variable
            glp_set_col_bnds(_solver, result, GLP_FX, lower_bound, upper_bound);
        } else {
            glp_set_col_bnds(_solver, result, GLP_DB, lower_bound, upper_bound);
        }
    } else if (lower_bounded && !upper_bounded) {
        // only lower bounded
        glp_set_col_bnds(_solver, result, GLP_LO, lower_bound, upper_bound);
    } else if (!lower_bounded && upper_bounded) {
        // only upper bounded
        glp_set_col_bnds(_solver, result, GLP_UP, lower_bound, upper_bound);
    } else {
        // free variable (unbounded on both sides)
        glp_set_col_bnds(_solver, result, GLP_FR, lower_bound, upper_bound);
    }

    return result;
}


void GlpkSolver::setObjectiveDirection(MilpSolver::ObjectiveType dir)
{
    if (dir == ObjectiveType::MAXIMIZE) {
        glp_set_obj_dir(_solver, GLP_MAX);
    } else {
        glp_set_obj_dir(_solver, GLP_MIN);
    }
}

void GlpkSolver::setObjective(const std::vector<MilpSolver::VarType>& variables,
                              const std::vector<double>& coeffs)
{
    for (std::size_t i = 0; i < variables.size(); ++i) {
        glp_set_obj_coef(_solver, variables[i], coeffs[i]);
    }
}

double GlpkSolver::getObjectiveValue()
{
    return glp_get_obj_val(_solver);
}


void GlpkSolver::addConstraint(const std::vector<MilpSolver::VarType>& variables,
                               const std::vector<double>& coeffs,
                               bool lower_bounded,
                               double lower_bound,
                               bool upper_bounded,
                               double upper_bound)
{
    if (_max_constraint >= _curr_constraint) {
        glp_add_rows(_solver, 5);
        _max_constraint += 5;
    }

    std::vector<MilpSolver::VarType> vars(variables.size() + 1);
    std::vector<double> cs(coeffs.size() + 1);

    for (std::size_t i = 0; i < variables.size(); ++i) {
        vars[i + 1] = variables[i];
        cs[i + 1] = coeffs[i];
    }

    glp_set_mat_row(_solver, _curr_constraint, static_cast<int>(variables.size()), &vars[0], &cs[0]);

    if (lower_bounded && upper_bounded) {
        if (lower_bound == upper_bound) {
            glp_set_row_bnds(_solver, _curr_constraint, GLP_FX, lower_bound, upper_bound);
        } else {
            glp_set_row_bnds(_solver, _curr_constraint, GLP_DB, lower_bound, upper_bound);
        }
    } else if (lower_bounded && !upper_bounded) {
        glp_set_row_bnds(_solver, _curr_constraint, GLP_LO, lower_bound, upper_bound);
    } else if (!lower_bounded && upper_bounded) {
        glp_set_row_bnds(_solver, _curr_constraint, GLP_UP, lower_bound, upper_bound);
    } else {
        LOG(ERROR) << "Trying to add unbounded constraint!";
        std::exit(-1);
    }

    ++_curr_constraint;
}


TruthValue GlpkSolver::solve()
{
    glp_iocp opt;
    glp_init_iocp(&opt);
    opt.presolve = GLP_ON;
    opt.mir_cuts = GLP_ON;
    opt.fp_heur = GLP_ON;
    if (!_verbose) {
        opt.msg_lev = GLP_MSG_OFF;
    }
    glp_intopt(_solver, &opt);

    const int status = glp_mip_status(_solver);

    if (status == GLP_FEAS || status == GLP_OPT) return TruthValue::TRUE;
    else if (status == GLP_NOFEAS) return TruthValue::FALSE;
    else return TruthValue::UNKNOWN;
}


TruthValue GlpkSolver::solve(void(*callback)(glp_tree*, void*), void* info)
{
    glp_smcp lp_opt;
    glp_init_smcp(&lp_opt);
    lp_opt.msg_lev = GLP_MSG_OFF;
    glp_simplex(_solver, &lp_opt);
    val_assert(glp_get_status(_solver) == GLP_OPT);

    glp_iocp opt;
    glp_init_iocp(&opt);
    opt.mir_cuts = GLP_ON;
    opt.fp_heur = GLP_ON;
    opt.cb_func = callback;
    opt.cb_info = info;
    if (!_verbose) opt.msg_lev = GLP_MSG_OFF;

    glp_intopt(_solver, &opt);

    const int status = glp_mip_status(_solver);

    if (status == GLP_FEAS || status == GLP_OPT) return TruthValue::TRUE;
    else if (status == GLP_NOFEAS) return TruthValue::FALSE;
    else return TruthValue::UNKNOWN;
}


int GlpkSolver::getIntValue(MilpSolver::VarType var)
{
    return static_cast<int>(std::round(glp_mip_col_val(_solver, var)));
}


double GlpkSolver::getDoubleValue(MilpSolver::VarType var)
{
    return glp_mip_col_val(_solver, var);
}


#endif

} // end namespace hqspre
