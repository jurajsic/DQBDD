/********************************************************************************************
modelrebuilder.cpp -- Copyright (c) 2016, Sven Reimer

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

#include "modelrebuilder.h"
#include "preprocessor.h"

namespace antom 
{

  ModelRebuilder::ModelRebuilder(Preprocessor* prepro):
	_prepro(prepro),
	_control(prepro->_control),
	_ca(prepro->_ca),
	_model(prepro->_model),
	_replacedBy(),
	_elimClauses(),
	_currentPos(0),
	_blockedClauses(),
	_eliminationSteps()
  {
  } 

  uint32_t ModelRebuilder::IsReplaced(uint32_t lit) const 
  {
	uint32_t var = lit>>1;
	bool polarity = lit&1;
	assert( var < _replacedBy.size() );
	uint32_t rlit(_replacedBy[var]);

	while ( _replacedBy[rlit>>1] != 0 )
	  {
		rlit = _replacedBy[rlit>>1]^(rlit&1);
	  }

	return rlit^polarity;
  } 

  void ModelRebuilder::ExtendModel(void)
  {
	_currentPos = static_cast<int32_t>(_elimClauses.size())-1;

	// Preceed methods in reverse order
	for ( int32_t i = static_cast<int32_t>(_eliminationSteps.size())-1; i >= 0; --i )
	  {
		switch (_eliminationSteps[i].preproOp)
		  {
		  case VAR_ELIMINATION : 
			{
			  RestoreVariableElimination ( _eliminationSteps[i].data );
			  break;
			}
		  case VAR_REPLACEMENT :
			{
			  uint32_t replacedVar = static_cast<uint32_t>(_eliminationSteps[i].data);
			  uint32_t replacingLit = _replacedBy[replacedVar];
			  uint32_t replacingModel = _model[replacingLit>>1];
			  assert ( replacingModel != 0 );

			  bool polarity = (replacingLit^replacingModel)&1;
			  _model[replacedVar] = (replacedVar<<1)^polarity;
			  break;
			}
		  case BLOCKED_CLAUSE :
			{
			  ReconstructBCE ( _eliminationSteps[i].data );			  
			  break;
			}
		  default :
			{
			  assert(false);
			  exit(0);
			}
		  }
	  }

	assert( _currentPos == -1 );
  }

  void ModelRebuilder::AddVarEquivalence(uint32_t replacedLit, uint32_t replacingLit)
  {
	_eliminationSteps.push_back( PreproMethod( VAR_REPLACEMENT, replacedLit>>1 ) );
	_replacedBy[replacedLit>>1] = replacingLit^(replacedLit&1);
  }


  void ModelRebuilder::AddVarElimination(uint32_t eliminatedVar)
  {
	PreproMethod currentVarElim = PreproMethod(VAR_ELIMINATION, _elimClauses.size());

	uint32_t poslit(eliminatedVar<<1);
	uint32_t neglit(poslit^1);

	size_t possize(_prepro->OccurenceCount(poslit));
	size_t negsize(_prepro->OccurenceCount(neglit));

	assert(!_prepro->_assignment[poslit] && !_prepro->_assignment[neglit]);

	// Preserve model
	uint32_t firstlit( poslit );
	uint32_t secondlit( neglit );

	// Determine the polarity with less occurences
	if (negsize < possize)
	  {
		firstlit = neglit;
		secondlit = poslit;
	  }
  
	// Store binary clauses of one polarity (the one with less occurences)
	for (size_t n = 0; n != _prepro->_binaries[firstlit].size(); ++n)
	  {
		if ( _prepro->_binaries[firstlit][n].IsLearnedBinary() )
		  { continue; }
		_elimClauses.push_back(firstlit);
		_elimClauses.push_back(_prepro->_binaries[firstlit][n].GetSecondLit());
		// Size of the clause
		_elimClauses.push_back(2);
	  }

	// Store n-nary clauses of one polarity (the one with less occurences)
	size_t osize(_prepro->_occur[firstlit].size());
	for (size_t n = 0; n != osize; ++n)
	  {
		Clause& clause = _ca[_prepro->_occur[firstlit][n]];

		if (clause.IsLearned())
		  { continue; }

		uint32_t first = static_cast<uint32_t>(_elimClauses.size());
		uint32_t csize = static_cast<uint32_t>(clause.size());

		for (uint32_t pos = 0; pos != csize; ++pos)
		  {
			_elimClauses.push_back(clause[pos]);
			// Store "firstlit" on first position;
			if (clause[pos] == firstlit)
			  { 
				_elimClauses[first+pos] = _elimClauses[first];
				_elimClauses[first] = firstlit;
			  }
		  }
		// Size of the clause 
		_elimClauses.push_back(csize);
	  }

	// Store opposite polarity of the literal
	_elimClauses.push_back(secondlit);
	_elimClauses.push_back(1);

	currentVarElim.endpos = _elimClauses.size();
	_eliminationSteps.push_back(currentVarElim);

  }
  void ModelRebuilder::AddBlockedClause(uint32_t blockingLit, const std::vector< uint32_t >& blockingClause)
  {
	_eliminationSteps.push_back(PreproMethod(BLOCKED_CLAUSE, _blockedClauses.size() ));

	_blockedClauses.push_back( BlockedClause( blockingLit, blockingClause ) );
  }

  // Adapted (and translated in readable code) from minisat/SATelite
  // Proceed deleted clauses during resolution and reset model assignments if necessary
  // Clauses are stored in "_elimClauses" in following manner:
  // 1. If "var" is eliminated, the complete occurence list of one polarity is stored
  // 2. The occurence of "var" is stored in front of the clause
  // 3. After each clause the size of the clause is stored
  // Ex.: clause "( v_1 + v_2 + var )" => "var v_2 v_1 3"
  // The vector is proceeded backwards, checking if all clauses are satisfied without "var"
  // If one clause is not satisfied, the assignment of var is flipped to the opposite polarity
  void ModelRebuilder::RestoreVariableElimination(size_t finalpos)
  {
	uint32_t clausesize( 0 );

	assert( static_cast<int32_t>(finalpos) <= _currentPos );

	while ( static_cast<int32_t>(finalpos) <= _currentPos )
	  {
		// Did we have a clause size argument?
		if ( clausesize == 0 )
		  { 
			clausesize = _elimClauses[_currentPos]; --_currentPos;
		  }

		uint32_t pos = _elimClauses[_currentPos];

		// Model of literal correctly setted? Skip clause
		if ( _model[pos>>1] == pos )
		  {
			// Goto next clause
			_currentPos -= (clausesize-1);
			clausesize = 0;
		  }
		// Model of literal incorrectly setted? Go to next literal
		else if ( _model[pos>>1] != pos )
		  { 
			// Reached last element in clause and clause is not satisfied? Change the assignment!
			if ( --clausesize == 0 ) 
			  { 
				_model[pos>>1] = pos;
				//std::cout << "change model of " << (pos>>1) << " to " << Lit(pos) << std::endl;
			  }
		  }
		// variable was eliminated before setting
		// Set model for variable to true and skip clause
		else if ( _model[pos>>1] == 0 )
		  {
			_model[pos>>1] = pos;
			//std::cout << "change model of " << (pos>>1) << " to " << Lit(pos) << std::endl;
		  }
		--_currentPos;
	  }
  }

  void ModelRebuilder::ReconstructBCE(size_t pos)
  {
	const std::vector< uint32_t >& clause = _blockedClauses[pos].bClause;

	uint32_t i = 0;
	for ( ; i != clause.size(); ++i )
	  {
		if ( _model[clause[i]>>1] == clause[i] )
		  { break; }
	  }

	if ( i == clause.size() )
	  {
		uint32_t blockinglit( _blockedClauses[pos].bLit );
		_model[blockinglit>>1] = blockinglit;
	  }
  }

  void ModelRebuilder::ClearRestoreData(uint32_t begin, uint32_t end)
  {
	// First clear replacelist
	for ( uint32_t v = 0; v <= end; ++v )
	  {
		if ( (v >= begin) && (v <= end) )
		  {
			_replacedBy[v] = 0;
		  }
#ifndef NDEBUG
		else 
		  {
			// Consistency checks
			assert( _replacedBy[v] == 0 || (_replacedBy[v]>>1) < begin || (_replacedBy[v]>>1) > end );
		  }
#endif
	  }

	// Now clear "_elimClauses"
	uint32_t clausesize(0);
	uint32_t clausebegin(0);
	bool unit(false);
	bool deleteclauses(false);
	for ( int32_t i = (int32_t)_elimClauses.size()-1; i >= 0; --i )
	  {
		// Did we have a clause size argument?
		if ( clausesize == 0 )
		  { 
			clausesize = _elimClauses[i];
			clausebegin = i;
			--i; 

			unit = (clausesize==1);

			if ( deleteclauses )
			  { 
				// reached a new unit, clear deleteclauses status
				if ( unit )
				  { deleteclauses = false; }
				// Otherwise delete clausesize entry
				else
				  { _elimClauses[i+1] = 0; }
			  }
			
		  }

		if ( deleteclauses )
		  { 
			_elimClauses[i] = 0; --clausesize; continue;
		  }

		uint32_t var( _elimClauses[i] >> 1);
		// reached deleted var?
		if ( (var>=begin) && (var<=end) )
		  {
			if ( unit )
			  {
				// We have reached a clause from a deleted var
				// Delete all these clauses
				deleteclauses = true;
				_elimClauses[i] = 0;
				_elimClauses[i+1] = 0;
			  }
			else
			  {
				// We found a deleted var in list of a clause from another deleted var
				// Just delete this clause
				for ( int32_t j = clausebegin; j > i; --j )
				  { 
					_elimClauses[j] = 0; 
				  }
				for ( ; clausesize == 0 ; --clausesize )
				  { 
					_elimClauses[i] = 0; --i;
				  }
			  }
		  }

		--clausesize;
	  }

	// Finally update the elimantionstep positions for blocking clauses and variable eliminations
	size_t elimOld = 0;
	size_t elimNew = 0;
	uint32_t blockOld = 0;
	uint32_t blockNew = 0;
	uint32_t stepOld = 0;
	uint32_t stepNew = 0;

	for (; stepOld != _eliminationSteps.size(); ++stepOld)
	  {
		_eliminationSteps[stepNew] = _eliminationSteps[stepOld];

		if (_eliminationSteps[stepNew].preproOp == VAR_ELIMINATION)
		  {
			_eliminationSteps[stepNew].data = elimNew;
			for ( ; elimOld != _eliminationSteps[stepNew].endpos; ++elimOld )
			  {
				if ( _elimClauses[elimOld] == 0 )
				  { continue;}

				_elimClauses[elimNew] = _elimClauses[elimOld];
				++elimNew;
			  }
			
			// If we have deleted all clauses of the resolution step, delete the complete elimination step
			// Otherwise update the endposition
			if (_eliminationSteps[stepNew].data != elimNew)
			  {
				_eliminationSteps[stepNew].endpos = elimNew;
				++stepNew;
			  }
		  }
		else if (_eliminationSteps[stepNew].preproOp == BLOCKED_CLAUSE)
		  {
			const std::vector< uint32_t >& clause = _blockedClauses[blockOld].bClause;
			// Should we delete this clause?
			uint32_t i = 0;
			for (; i != clause.size(); ++i )
			  {
				if (((clause[i]>>1)>=begin) && ((clause[i]>>1)<=end))
				  {
					break;
				  }
			  }

			// Keep clause? -> Update elimination step and blocking clause positions
			if( i == clause.size() )
			  {
				_blockedClauses[blockNew] = _blockedClauses[blockOld];
				_eliminationSteps[stepNew].data = blockNew;
				++blockNew;
				++stepNew;
			  }
			++blockOld;
		  }
	  }

	_eliminationSteps.resize(stepNew);
	_elimClauses.resize(elimNew);
	_blockedClauses.resize(blockNew);
  }

}
