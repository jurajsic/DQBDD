#ifndef ENCODINGS_H
#define ENCODINGS_H

#include "core/Solver.h"
//#include "core/SolverTypes.h"
#include "Settings.h"

namespace Glucose{

class Encodings
{

  public:
    Encodings():
      _settings(nullptr)
      {}


    bool generateEncoding();


//private:
    Settings* _settings;

    // Warners adder encoding
    void genWarnersHalf(Lit &a, Lit &b, Lit &carry, Lit &sum, int comp, Solver &S, vec<Lit> &lits);
    void genWarnersFull(Lit &a, Lit &b, Lit &c, Lit &carry, Lit &sum, int comp, Solver &S, vec<Lit> &lits);
    void genWarners(vec<long long> &weights, vec<Lit> &blockings, long long max, int k, int comp, Solver &S, const Lit zero, vec<Lit> &lits, vec<Lit> &linkingVar);
    void genWarners0(vec<long long> &weights, vec<Lit> &blockings, long long max, long long k, int comp, Solver &S, vec<Lit> &lits, vec<Lit> &linkingVar);

    // Bailleux totalizer encoding
    void genBailleux(vec<long long> &weights, vec<Lit> &blockings, long long total, Lit zero, Lit one, int comp, Solver &S, vec<Lit> &lits, vec<Lit> &linkingVar, long long UB);
    void genBailleux0(vec<long long> &weights, vec<Lit> &blockings, long long max, long long k, int comp, Solver &S, vec<Lit> &lits, vec<Lit> &linkingVar);

    // [Asin et. al 2011] encoding
    void sComparator(Lit &a, Lit &b, Lit &c1, Lit &c2, int comp, Solver &S, vec<Lit> &lits);
    void genSMerge(vec<Lit> &linkA, vec<Lit> &linkB, Lit zero, Lit one, int comp, Solver &S, vec<Lit> &lits, vec<Lit> &linkingVar, long long UB);
    void genKCard(vec<long long> &weights, vec<Lit> &blockings, long long total, Lit zero, Lit one, int comp, Solver &S, vec<Lit> &lits, vec<Lit> &linkingVar, long long UB);
    void genAsin(vec<long long> &weights, vec<Lit> &blockings, long long max, long long k, int comp, Solver &S, vec<Lit> &lits, vec<Lit> &linkingVar);

    // Ogawa Modulo Totalizer
    void genOgawa(long long weightX, vec<Lit> &linkingX, long long weightY, vec<Lit> &linkingY, long long &total, long long divisor, Lit zero, Lit one, int comp, Solver &S, vec<Lit> &lits, vec<Lit> &linkingVar, long long UB);
    void genOgawa(vec<long long> &weights, vec<Lit> &blockings, long long &total, long long divisor, Lit zero, Lit one, int comp, Solver &S, vec<Lit> &lits, vec<Lit> &linkingVar, long long UB);
    void genOgawa0(int &card, vec<long long> &weights, vec<Lit> &blockings, long long max, long long k, long long &divisor, int comp, Solver &S, vec<Lit> &lits, vec<Lit> &linkingVar);

    // BailleuxW2
    void genBailleuxW2(vec<long long> &weights, vec<Lit> &blockings, long long total, Lit zero, Lit one, int comp, Solver &S, vec<Lit> &lits, vec<Lit> &linkingVar, vec<long long> &linkingW, long long UB);
    void genBailleuxW20(vec<long long> &weights, vec<Lit> &blockings, long long max, long long k, int comp, Solver &S, vec<Lit> &lits, vec<Lit> &linkingVar, vec<long long> &linkingWeight);
    void genCCl(Lit a, Solver &S, vec<Lit> &lits, Var varZero);
    void genCCl(Lit a, Lit b, Solver &S, vec<Lit> &lits, Var varZero);
    void genCCl(Lit a, Lit b, Lit c, Solver &S, vec<Lit> &lits, Var varZero);
    void genCCl1(Lit a, Lit b, Lit c, Solver &S, vec<Lit> &lits, Var varZero);
    void genCCl(Lit a, Lit b, Lit c, Lit d, Solver &S, vec<Lit> &lits, Var varZero);
    void genCCl(Lit a, Lit b, Lit c, Lit d, Lit e, Solver &S, vec<Lit> &lits, Var varZero);

    // Weighted MaxSAT Totalizer
    void genKWMTO(vec<long long> &weights, vec<Lit> &blockings, vec<long long> &weightsTable, int from, int to, int div, Lit zero, vec<Lit> &lower, vec<long long> &lowerW, vec<Lit> &upper, vec<long long> &upperW, Solver &S, long long ub, vec<Lit> &lits, Var varZero);
    void genKWMTO0(int &card, vec<long long> &weights, vec<Lit> &blockings, long long max, long long k, vec<long long> &divisors, Solver &S, vec<Lit> &lits, vec<vec<Lit> > &linkingVars, vec<vec<long long> > &linkingWeights);

    // Mixed Radix Weighted Totalizer
    void genMRWTO(vec<long long> &weights, vec<Lit> &blockings, vec<long long> &weightsTable, int from, int to, vec<long long> &divisors, Lit zero, vec<vec<Lit> > &linkingVars, vec<vec<long long> > &linkingWeights, Solver &S, long long ub, vec<Lit> &lits, Var varZero);
    void genMRWTO0(int &card, vec<long long> &weights, vec<Lit> &blockings, long long max, long long k, vec<long long> &divisors, Solver &S, vec<Lit> &lits, vec<vec<Lit> > &linkingVars, vec<vec<long long> > &linkingWeights);
};
}

#endif // ENCODINGS_H
