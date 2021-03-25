// $Id: milp_solver.hpp 2644 2019-09-07 20:46:54Z wimmer $

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

#ifndef HQSPRE_MILP_SOLVER_HPP_
#define HQSPRE_MILP_SOLVER_HPP_

#include <vector>

#include "auxil.hpp"

struct glp_tree;
struct glp_prob;

namespace hqspre {

class MilpSolver
{
   public:
    using VarType = int;
    enum class VarSort
    {
        INT,
        REAL
    };
    enum class ObjectiveType
    {
        MINIMIZE,
        MAXIMIZE
    };

    MilpSolver()          = default;
    virtual ~MilpSolver() = default;

    virtual VarType addVariable(VarSort sort, bool lower_bounded, double lower_bound, bool upper_bounded,
                                double upper_bound)
        = 0;

    virtual void addConstraint(const std::vector<VarType>& variables, const std::vector<double>& coeffs,
                               bool lower_bounded, double lower_bound, bool upper_bounded, double upper_bound)
        = 0;

    virtual void       setObjectiveDirection(ObjectiveType dir)                                               = 0;
    virtual void       setObjective(const std::vector<VarType>& variables, const std::vector<double>& coeffs) = 0;
    virtual double     getObjectiveValue()                                                                    = 0;
    virtual TruthValue solve()                                                                                = 0;

    virtual int    getIntValue(VarType var)    = 0;
    virtual double getDoubleValue(VarType var) = 0;

    void setVerbosity(bool verbose) { _verbose = verbose; }

   protected:
    bool _verbose = true;
};

class GlpkSolver : public MilpSolver
{
   public:
    GlpkSolver();
    virtual ~GlpkSolver();
    GlpkSolver(const GlpkSolver&) = delete;
    GlpkSolver(GlpkSolver&&)      = delete;
    GlpkSolver& operator=(const GlpkSolver&) = delete;
    GlpkSolver& operator=(GlpkSolver&&) = delete;

    VarType addVariable(VarSort sort, bool lower_bounded, double lower_bound, bool upper_bounded,
                        double upper_bound) override;

    void addConstraint(const std::vector<VarType>& variables, const std::vector<double>& coeffs, bool lower_bounded,
                       double lower_bound, bool upper_bounded, double upper_bound) override;

    void   setObjectiveDirection(ObjectiveType dir) override;
    void   setObjective(const std::vector<VarType>& variables, const std::vector<double>& coeffs) override;
    double getObjectiveValue() override;

    TruthValue solve() override;
    TruthValue solve(void (*callback)(glp_tree*, void*), void* info);
    int        getIntValue(VarType var) override;
    double     getDoubleValue(VarType var) override;

   private:
    glp_prob*    _solver;
    unsigned int _max_constraint;
    unsigned int _curr_constraint;
};

}  // end namespace hqspre

#endif
