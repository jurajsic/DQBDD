#ifndef SOFTCLAUSE_H
#define SOFTCLAUSE_H

#include <vector>
#include <core/SolverTypes.h>
#include <mtl/Vec.h>

namespace Glucose{

// For MAX-SAT we store the soft clauses separately
struct SoftClause
{

  // Constructor
  SoftClause(Lit relaxlit, const vec< Lit >& givenClause, long long int givenWeight = 1) :
    relaxationLit(relaxlit),
    weight(givenWeight)
  {
     givenClause.memCopyTo(clause);
  }

  // Destructor
  ~SoftClause(void) {}

  Lit relaxationLit;
  vec< Lit > clause;
  long long int weight;

  static bool bigger(const SoftClause* s1, const SoftClause* s2) { return s1->weight > s2->weight; }

};
}

#endif // SOFTCLAUSE_H
