
/********************************************************************************************
bce.h -- Copyright (c) 2016, Sven Reimer

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

#ifndef ANTOM_BCE_H_
#define ANTOM_BCE_H_

#include "antombase.h"
#include "watcher.h"
#include "statistics.h"
#include "clause.h"
#include "varheap.h"

namespace antom {

  class ModelRebuilder;

  class BCE {
  public:

	explicit BCE(Preprocessor* prepro);
	// Performs blocked clause elimination,
	bool DoBlockedClauseElimination(void);

	// Try to add new binary implications which were originally blocked 
	void AddBlockedImplications(void);

  private:

	// Copy constructor.
    BCE (const BCE&) = default;

    // Assignment operator.
    BCE& operator = (const BCE&) = default;

	// Check whether watchlist of "lit" contains blocked clauses
	bool CheckBlockedLit(uint32_t lit);
	// Check whether the binary clause (lit, secondlit) is blocked by "lit"
	bool CheckBlocked(uint32_t lit, uint32_t secondlit) const;
	// Check whether the clause "clause" is blocked by "lit" 
	bool CheckBlocked(uint32_t lit, const Clause& clause) const;

	// Remove the blocked clauses from database and store them in "b_blockedClauses"
	void RemoveBlocked(uint32_t bLit, const Watcher& binClause);
	void RemoveBlocked(uint32_t bLit, CRef clause); 

	Preprocessor* _prepro;
	Settings* _setting;
	ModelRebuilder* _rebuilder;
	ClauseAllocator& _ca;
	std::vector< uint32_t >& _model;
	Statistics& _statistics;

	VarHeap<helper::AscendingOrder, size_t> _candidates;
  };

}

#endif
