
/********************************************************************************************
modelrebuilder.h -- Copyright (c) 2016, Sven Reimer

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

#ifndef ANTOM_MODELREBUILDER_H_
#define ANTOM_MODELREBUILDER_H_

#include "antombase.h"
#include "clause.h"

enum PreproOperation {
  NO_OPERATION,
  VAR_ELIMINATION,
  VAR_REPLACEMENT,
  BLOCKED_CLAUSE,
};

namespace antom {

  class Control;
  class Preprocessor;

  class ModelRebuilder {

  public:

	explicit ModelRebuilder(Preprocessor* prepro);

	uint32_t IsReplaced(uint32_t lit) const;

	// Restore model for prepro operations which does not preseve logical equivalence
	// This method should be called before model is returned
	// Assumes that model of the solver core is already set to "p_model"
	void ExtendModel(void);

	void AddVarEquivalence(uint32_t replacedVar, uint32_t replacingLit);
	void AddVarElimination(uint32_t eliminatedVar);
	void AddBlockedClause(uint32_t blockingLit, const std::vector< uint32_t >& blockingClause);

	friend class Preprocessor;
	friend class UPLA;
	
  protected:

	struct BlockedClause 
	{
	BlockedClause(void) :
	  bLit(0),
		bClause()
	  {}

	BlockedClause(uint32_t blockingLiteral, std::vector< uint32_t > blockingClause) :
	  bLit(blockingLiteral),
		bClause(blockingClause)
	  {}

	  uint32_t bLit;
	  std::vector< uint32_t > bClause;
	};

	struct PreproMethod 
	{
	PreproMethod(void) :
	  preproOp(NO_OPERATION),
		data(0),
		endpos(0)
	  {}

	PreproMethod(PreproOperation op, size_t input) :
	  preproOp(op),
		data(input),
		endpos(0)
	  {}
		
	  PreproOperation preproOp;

	  // represents :
	  // beginning position in _elimclauses for VAR_ELIMINATION,
	  // replaced variable for VAR_REPLACEMENT,
	  // position in _blockedClauses for BLOCKED_CLAUSE
	  size_t data;
	  // endposition in _elimclauses for VAR_ELIMINATION
	  size_t endpos;
	};


	void ReconstructBCE(size_t pos);
	void RestoreVariableElimination(size_t finalpos);

	void ClearRestoreData(uint32_t begin, uint32_t end);

	Preprocessor* _prepro;
	Control* _control;
	ClauseAllocator& _ca;

	std::vector< uint32_t >& _model;

	// Stores equivalences
	std::vector< uint32_t > _replacedBy;
	// Stores Variable elimination
	std::vector< uint32_t > _elimClauses;
	int32_t _currentPos;
	// Stores Blocked clauses
	std::vector< BlockedClause > _blockedClauses;
	// Stores all operations which does not preseve logical equivalence in order
	std::vector< PreproMethod > _eliminationSteps;

  private :
	// Copy constructor.
    ModelRebuilder (const ModelRebuilder&) = default;

    // Assignment operator.
    ModelRebuilder& operator = (const ModelRebuilder&) = default;
  };
}

#endif
