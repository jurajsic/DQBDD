/********************************************************************************************
softclause.h -- Copyright (c) 2016-2017, Sven Reimer, Tobias Paxian

Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software 
without restriction, including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons 
to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

********************************************************************************************/

#ifndef ANTOM_SOFTCLAUSE_H_
#define ANTOM_SOFTCLAUSE_H_

#include <vector>

namespace antom
{
    // For MAX-SAT we store the soft clauses separately
	struct SoftClause 
	{

	  // Constructor
	  SoftClause(uint32_t relaxlit, const std::vector< uint32_t >& clause, uint64_t weight = 1) :
		relaxationLit(relaxlit),
		clause(clause),
		lastassignment(0),
		contra(0),
		weight(weight)
	  {};

	  // Destructor
	  ~SoftClause(void) {}

	  uint32_t relaxationLit;
	  std::vector<uint32_t> clause;

	  // Stores last found assignment for the relaxation lit
	  uint32_t lastassignment;

	  uint32_t contra;
	  uint64_t weight;

	  static bool bigger(const SoftClause* s1, const SoftClause* s2) { return s1->weight > s2->weight; }

	};

  struct SoftClauseContraSorter {
	bool operator()(SoftClause* s1, SoftClause* s2) const;
  };
  
  inline bool SoftClauseContraSorter::operator()(SoftClause* s1, SoftClause* s2) const
  { return s1->contra > s2->contra; }
  
    /*
  struct SoftClausePopCountSorter {
    bool operator()(SoftClause* s1, SoftClause* s2) const;
  };

  inline bool SoftClausePopCountSorter::operator()(SoftClause* s1, SoftClause* s2) const
  { return __builtin_popcountll(s1->weight) > __builtin_popcountll(s2->weight); }
*/
}

#endif
