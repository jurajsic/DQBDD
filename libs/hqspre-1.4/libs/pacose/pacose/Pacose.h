#ifndef PACOSE_H
#define PACOSE_H

#include "Encodings.h"
#include "Settings.h"
#include <iostream>

namespace Glucose{

struct SoftClause;

// Pacose is an extension to the glucose solver
// containing
/**
 * @brief The Pacose class
 *          Is an extension to the glucose syrup solver containing:
 *      1.  (Partial) (Weighted) MaxSAT Solver
 *          Derived from QMaxSAT (2017 Competition Version)
 *      2.  Solving Functions for Applications
 *          Derived from Antom SAT Solver (University of Freiburg 2011-2018✝, IRA)
 *          Needed for Phaeton as formal method backend solver.
 */
class Pacose : public Solver
{
    //    friend class Settings;

public:
    Pacose();

    ~Pacose();



    lbool MaxSolve(IntOption comp, int card, IntOption pmodel);


    // Functions needed for external programms including Pacose as libary!

    /**
     * @brief AddSoftClause
     * @param Clause
     * @param weight
     */
    void AddSoftClause(vec<Lit>& Clause, long long int weight);

    /**
     * @brief DeduceAssumptions
     *          Try to find trivial contradictions for given assumptions.
     * @return l_True       no trivial contradictions found
     *         l_False      trivial contradictions found
     *         l_Unknown    else...
     */
    lbool DeduceAssumptions(const vec<Lit>& externalAssumptions);

    // Resets the SAT solving core. The following status flags/variables remain untouched:
    // -- The SAT solving threads unique ID number: "_id".
    // -- The pointer to the "Control" object: "_control".
    // -- The number of variables for which memory has been reserved: "_variables".
    void Reset(void);

    void InstanceReset(void);

    void qMaxSatHeuristic();

    void DGPW18Heuristic();

    Settings _settings;
    Encodings _encodings;

    vec< SoftClause* > _softClauses;

    inline int _numberOfSoftVars()
    {
        return _softClauses.size();
    }

private:

    /**
     * @brief wbSort    Creates a copy of weights and blockings and saves them into sweights and sblockings.
     * @param weights
     * @param blockings
     * @param sweights
     * @param sblockings
     */
    void wbSort(vec<long long> &weights, vec<Lit> &blockings, vec<long long> &sweights, vec<Lit> &sblockings);

    /**
     * @brief wbFilter  Filters out all non possible weights smaller than o-Value.
     * @param UB
     * @param S
     * @param lits
     * @param weights
     * @param blockings
     * @param sweights
     * @param sblockings
     */
    void wbFilter(long long UB, Solver &S, vec<Lit> &lits, vec<long long> &weights, vec<Lit> &blockings, vec<long long> &sweights, vec<Lit> &sblockings);

    /**
     * @brief genCardinals
     * @param card
     * @param comp
     * @param weights
     * @param blockings
     * @param max
     * @param k
     * @param divisor
     * @param S
     * @param lits
     * @param linkingVar
     * @param linkingWeight
     * @param divisors
     * @param linkingVars
     * @param linkingWeights
     */
    void genCardinals(int &card, int comp, vec<long long> &weights, vec<Lit> &blockings, long long max, long long k, long long &divisor, vec<Lit> &lits, vec<Lit> &linkingVar, vec<long long> &linkingWeight, vec<long long> &divisors, vec<vec<Lit> > &linkingVars, vec<vec<long long> > &linkingWeights);
    long long sumWeight(const vec<long long> &weights) const;

    // Vars
    int _nbOfOrigVars;

    enum class SATSolverResult {
            UNKNOWN                 =  0,
            SAT                     = 10,
            UNSAT                   = 20,
            UNSAT_WITH_ASSUMPTION   = 30,
    } _extendedResult;

};
}

#endif // PACOSE_H
