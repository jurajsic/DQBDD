/********************************************************************************************
debug.cpp -- Copyright (c) 2014-2016, Sven Reimer

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

#include <algorithm>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

#include "antom.h"
#include "sorter.h"
#include "core.h"
#include "preprocessor.h"

namespace antom {
  
  // prints literals of clause
  void Clause::PrintLits(void) const
  {
	const uint32_t* literals = GetLiterals();
	assert( literals != NULL );
	
	for( uint32_t pos = 0; pos != _header.length; ++pos)
	  {
		std::cout << helper::Lit(literals[pos]) << " "; 
	  }
	
	std::cout << std::endl;
  }

  void Clause::Print(void) const
  {
	PrintLits();
	std::cout << "length: " << size() << " lbd: " << Lbd() << " activity: " << Activity() << " reloced: " << Reloced() << std::endl;
  }

  void Core::PrintClause(Clause* clause, bool assignment) const
  {
	assert( clause != NULL );
	uint32_t size = clause->size();
	for( uint32_t pos = 0; pos != size; ++pos)
	  {
		std::cout << helper::Lit( (*clause)[pos]) << " "; 
		if( assignment )
		  { std::cout << "[" << IsAssigned((*clause)[pos]) << "] "; }
	  }
	std::cout << std::endl;
  }

  bool Core::CheckClauses() const
  {
	bool okay = true;
	for( uint32_t i = 0; i != _clauseDatabase.size(); ++i )
	  {
		uint32_t pos(0);
		const Clause& lits = _ca[_clauseDatabase[i]];
		for( ; pos != lits.size(); ++pos )
		  {
			if (lits[pos] == 0)
			  {
				std::cout << "clause " << i << " has corrupted content at position " <<  pos << " : ";
				lits.Print();
				okay = false;
				break;
			  }
		  }
	  }
	return okay;
  }
  

  // Returns:
  // "-1" if literal is assigned to false
  // " 1" if literal is assigned to true
  // " 0" if literal is unassigned
  int32_t Core::IsAssigned(uint32_t lit) const 
  {
  if( _assignment[lit] )
	  {
		assert( !_assignment[lit^1] );
		return 1;
	  }
	else if ( _assignment[lit^1] )
	  {
		return -1;
	  }
	return 0;
  }

  int32_t Preprocessor::IsAssigned(uint32_t lit) const
  {
	if( _assignment[lit] )
	  {
		assert( !_assignment[lit^1] );
		return 1;
	  }
	else if ( _assignment[lit^1] )
	  { 
		return -1;
	  }

	return 0;
  }

  void Preprocessor::PrintOccurenceList(uint32_t lit) const
  {
	for( uint32_t o = 0; o != _occur[lit].size(); ++o )
	  {
		std::cout << o << ": ";
		_ca[_occur[lit][o]].Print();
	  }
	std::cout << std::endl;
  }

  void Preprocessor::PrintWatchList(uint32_t lit) const
  {
	std::cout << "Watch list of " << helper::Lit(lit) << ":" << std::endl;
	for( uint32_t o = 0; o != _binaries[lit].size(); ++o )
	  {
		std::cout << o << ": ";
		_binaries[lit][o].Print(_ca);
	  }
	std::cout << std::endl;
  }

  void Preprocessor::PrintCompleteLists(uint32_t var) const
  {
	uint32_t poslit( var <<1 );
	uint32_t neglit( poslit^1);

	std::cout << helper::Lit(poslit) << ":" << std::endl;
	for( uint32_t o = 0; o != _binaries[poslit].size(); ++o )
	  {
		if( _binaries[poslit][o].IsBinary() )
		  {
			std::cout << o << ": " << helper::Lit(poslit) << " [" << IsAssigned(poslit) << "] " << helper::Lit( _binaries[poslit][o].GetSecondLit()) << " [" << IsAssigned(_binaries[poslit][o].GetSecondLit()) << "] " << _binaries[poslit][o].IsLearnedBinary() << std::endl ;
		  }
	  }
	std::cout << std::endl;
	for( uint32_t o = 0; o != _occur[poslit].size(); ++o )
	  {
		std::cout << o << ": ";
		const Clause& clause( _ca[_occur[poslit][o]] );

		for( uint32_t k = 0; k != clause.size(); ++k )
		  {
			std::cout << helper::Lit(clause[k]) << " [" << IsAssigned(clause[k]) << "] ";
		  }
		std::cout << "lbd: " << clause.Lbd() << " length: " << clause.size() << std::endl;
	  }

	std::cout << std::endl << helper::Lit(neglit) << ":" << std::endl;
	for( uint32_t o = 0; o != _binaries[neglit].size(); ++o )
	  {
		if( _binaries[neglit][o].IsBinary() )
		  {
			std::cout << o << ": " << helper::Lit(neglit) << " [" << IsAssigned(neglit) << "] " << helper::Lit( _binaries[neglit][o].GetSecondLit()) << " [" << IsAssigned(_binaries[neglit][o].GetSecondLit()) << "] "  << _binaries[neglit][o].IsLearnedBinary() << std::endl;
		  }
	  }
	std::cout << std::endl;
	for( uint32_t o = 0; o != _occur[neglit].size(); ++o )
	  {
		std::cout << o << ": ";
		const Clause& clause( _ca[_occur[neglit][o]] );

		for( uint32_t k = 0; k != clause.size(); ++k )
		  {
			std::cout << helper::Lit(clause[k]) << " [" << IsAssigned(clause[k]) << "] ";
		  }
		std::cout << "lbd: " << clause.Lbd() << " length: " << clause.size() << std::endl;
	  }
	std::cout << std::endl;
  }
  
  void Core::PrintBinaryList(uint32_t lit) const
  {
	std::cout << "Binary clauses of " << helper::Lit(lit) << ":" << std::endl;
	for( uint32_t o = 0; o != _watches[lit].size(); ++o )
	  {
		if(_watches[lit][o].IsBinary())
		  {
			std::cout << o << ": " << helper::Lit(lit) << " " << helper::Lit(_watches[lit][o].GetSecondLit()) << std::endl;
		  }
	  }
	std::cout << std::endl;
  }
  
  void Core::PrintWatchList(uint32_t lit) const
  {
	std::cout << "Watch list of " << helper::Lit(lit) << ":" << std::endl;
	for( uint32_t o = 0; o != _watches[lit].size(); ++o )
	  {
		std::cout << o << ": ";
		_watches[lit][o].Print(_ca);
	  }
	std::cout << std::endl;
  }

  // Watch lists in preprocessor should only contain binary clauses
  void Preprocessor::CheckWatchLists(void) const
  {
	for( uint32_t v = 1; v <= _variables; ++v )
	  {
		uint32_t lit = (v<<1)^1;
		for( uint32_t j = 0; j != 2; ++j )
		  {
			lit ^= 1;
			
			for( uint32_t o = 0; o != _binaries[lit].size(); ++o )
			  {
				if( !_binaries[lit][o].IsBinary() )
				  {
					std::cout << "found non binary: " << std::endl;
					_binaries[lit][o].Print(_ca);
				  }
				assert(_binaries[lit][o].IsBinary());
			  }
			std::cout << std::endl;
		  }
	  }
  }

  void Core::DumpCNF(bool printAssignment) const
  {
	// print binary
	for( uint32_t i = 2; i < _watches.size(); ++i )
	  {
		bool unsatisfied( _assignment[i^1] );
		for( uint32_t b = 0; b != _watches[i].size(); ++ b )
		  {
			if( _watches[i][b].IsBinary() && ( _watches[i][b].GetSecondLit() > i) )
			  {
				std::cout << helper::Lit(i) << " ";
				if( printAssignment )
				  { 
					std::cout << "[" << IsAssigned(i) << "] "; 
				  }
				std::cout << helper::Lit(_watches[i][b].GetSecondLit());
				if( printAssignment )
				  { 
					uint32_t seclit = _watches[i][b].GetSecondLit();
					unsatisfied &= _assignment[(seclit)^1];

					std::cout << " [" << IsAssigned(seclit) << "] "; 
				  }

				std::cout << " 0";
				if( printAssignment && unsatisfied )
				  { std::cout << " UNSAT"; }
				std::cout << std::endl;
			  }
		  }
	  }

	// print n-nary
	for( uint32_t i = 0; i != _clauseDatabase.size(); ++i )
	  {
		uint32_t pos(0);
		bool unsatisfied(true);
		const Clause& lits = _ca[_clauseDatabase[i]];
		for( ; pos != lits.size(); ++pos )
		  {
			std::cout << helper::Lit( lits[pos] ) << " "; 
			if( printAssignment )
			  { 
				std::cout << "[" << IsAssigned(lits[pos]) << "] "; 
				unsatisfied &= _assignment[lits[pos]^1];
			  }
		  }
		std::cout << "0";
		if( printAssignment && unsatisfied )
		  { std::cout << " UNSAT"; }
		std::cout << std::endl;
	  }
	std::cout << std::endl;
  }

  void Core::PrintDatabase(bool printWatches) const
  {
	_preprocessor->PrintDatabase( false, printWatches );
  }

  void Preprocessor::PrintDatabase(bool printOccur, bool printWatches) const
  {
	std::cout << "vars: " << std::endl;
	for( uint32_t v = 1; v <= _variables; ++v )
	  {
		std::cout << v << " del: " << (_deleted[v]?"T":"F") << " dt: " << (_donttouch[v]?"T":"F") << std::endl;
	  }

	if( printOccur )
	  {
		std::cout << std::endl << "occurence lists: " << std::endl;
		for( uint32_t v = 1; v <= _variables; ++v )
		  {
			std::cout << v << ":  (" << _binaries[(v<<1)].size() << ")"  << std::endl;
			for( uint32_t b = 0; b != _binaries[v<<1].size(); ++b )
			  {
				assert( _binaries[(v<<1)][b].IsBinary() );
				std::cout << v << " " << helper::Lit(_binaries[v<<1][b].GetSecondLit()) << std::endl;
			  }
			std::cout << "occur: " << _occur[v<<1].size() << std::endl;
			PrintOccurenceList(v<<1);

			std::cout << "-" << v << ": (" << _binaries[(v<<1)^1].size() << ")" << std::endl;
			for( uint32_t b = 0; b != _binaries[(v<<1)^1].size(); ++b )
			  {
				assert( _binaries[(v<<1)^1][b].IsBinary() );
				std::cout << "-" << v << " " << helper::Lit(_binaries[(v<<1)^1][b].GetSecondLit()) << std::endl;
			  }

			std::cout << "occur: " << _occur[(v<<1)^1].size() << std::endl;
			PrintOccurenceList((v<<1)^1);
			std::cout << std::endl;
		  }
	  }

	if( printWatches )
	  {
		std::cout << std::endl << "watch lists: " << std::endl;
		for( uint32_t v = 1; v <= _variables; ++v )
		  {
			for( uint32_t literal = (v<<1); literal < ((v<<1)+2); ++literal )
			  {
				std::cout << helper::Lit(literal) << ": (" << _binaries[literal].size() << ")"  << std::endl;

				for( size_t b = 0; b != _binaries[literal].size(); ++b )
				  {
					bool satisfied(false);
					if( _binaries[literal][b].IsBinary() )
					  {
						uint32_t secondlit(_binaries[literal][b].GetSecondLit());
						if( _assignment[literal] || _assignment[secondlit] )
						  { satisfied = true; }
						std::cout << helper::Lit(literal) << " [" << IsAssigned(literal) << "] " 
								  << helper::Lit(secondlit) << " [" << IsAssigned(secondlit) << "] ";
					  }
					else if ( _binaries[literal][b].IsTernary() )
					  {
						uint32_t secondlit(_binaries[literal][b].GetSecondLit());
						uint32_t thirdlit(_binaries[literal][b].GetThirdLit());
						
						std::cout << helper::Lit(literal) << " [" << IsAssigned(literal) << "] " 
								  << helper::Lit(secondlit) << " [" << IsAssigned(secondlit) << "] "
								  << helper::Lit(thirdlit) << " [" << IsAssigned(thirdlit) << "] ";
					  }
					else
					  {
						assert( _binaries[literal][b].IsClause() );
						Clause& clause(_ca[_binaries[literal][b].GetClause()]);
						uint32_t csize(clause.size());
						for( uint32_t pos = 0; pos != csize; ++pos )
						  {
							if( _assignment[clause[pos]] )
							  { satisfied = true; }
							std::cout << helper::Lit(clause[pos]) << " [" << IsAssigned(clause[pos]) << "] ";
						  }
					  }
					if( satisfied )
					  { std::cout << "SAT"; }
					std::cout << std::endl;
				  }
			  }
		  }
	  }

	std::cout << "clausedatabase: " << std::endl;

	for( uint32_t c = 0; c != _clauseDatabase.size(); ++c )
	  {
		std::cout << c << ": ";
		_ca[_clauseDatabase[c]].Print();
	  }
	std::cout << std::endl;
	std::cout << "decisionstack: " << std::endl;
	for( uint32_t i = 1; i < _dsEndIndex; ++i )
	  {
		std::cout << i << " " << helper::Lit(_decisionStack[i]) << std::endl;
	  }

	std::cout << std::endl;
  }

  // Consistency check whether all clauses with index > lastIndex >= _variables are not used
  // Return false if this property does not hold
  bool Core::CheckMaxIndex(uint32_t lastIndex) const
	{
	  // Check whether a clause with index > lastIndex is in a watchlist of allow indices
	  uint32_t i(0);
	  for( ; i <= lastIndex; ++i )
		{
		  for( uint32_t literal = (i<<1); literal < ((i<<1)+2); ++literal )
			{
			  
			  for( uint32_t j = 0; j != _watches[literal].size(); ++j )
				{
				  // Binary clause
				  if( _watches[literal][j].IsBinary())
					{

					  if( (_watches[literal][j].GetSecondLit()>>1) > lastIndex )
						{
						  std::cout << "forbidden binary index: " << (_watches[literal][j].GetSecondLit()>>1) << std::endl;
						  return false;
						}

					}
				  // N-nary claue
				  else
					{
					  const Clause& clause = _ca[_watches[literal][j].GetClause()];
					  uint32_t size = _ca[_watches[literal][j].GetClause()].size();
					  for ( uint32_t pos = 0; pos != size; ++pos )
						{
						  if( (clause[pos]>>1) > lastIndex )
							{
							  std::cout << "forbidden n-nary index: " << (clause[pos]>>1) << std::endl;
							  return false;
							}
						}
					}
				}	  
			}
		}

	  // Check whether all watch lists with index > lastIndex is empty
	  for( ; i <= _variables; ++i )
		{
		  if( !_watches[i<<1].empty() || !_watches[(i<<1)^1].empty() )
			{ 
			  std::cout << "index not empty: " << i << std::endl;
			  
			  return false; 
			}
		}
	  return true;
	}

#ifndef NDEBUG
  void Antom::CheckGates(void) const
  {
    if( !_antomSetting->encode01 && _antomSetting->networkType == BITONIC )
	  { return; }

	//std::cout << __func__ << std::endl;
	bool error(false);
	const std::vector< uint32_t >& newModel = Model();

	for( uint32_t foo = 0; foo != _debugNetwork.size(); ++foo)
	  {
		uint32_t input1 = newModel[_debugNetwork[foo].inputs[0]];
		uint32_t input2 = newModel[_debugNetwork[foo].inputs[1]];
		uint32_t output = newModel[_debugNetwork[foo].output];

		// And case
		if ( _debugNetwork[foo].type == ANDGATE )
		  {
			// input1 = FALSE OR input2 = FALSE => output = FALSE
			if ( (input1&1) != 0 )
			  {
				if ( (output&1) == 0 )
				  {
					std::cout << "ERROR: " << (output>>1) << " should be FALSE" << std::endl
							  << "AND of: " << (input1>>1) << " and " << (input2>>1) << std::endl;
					error = true;
				  }
			  }
			if ( (input2&1) != 0 )
			  {
				if ( (output&1) == 0 )
				  {
					std::cout << "ERROR: " << (output>>1) << " should be FALSE" << std::endl
							  << "AND of: " << (input1>>1) << " and " << (input2>>1) << std::endl;
					error = true;
				  }
			  }
			// input1 = TRUE AND input2 = TRUE => output = TRUE
			if ( ( (input1&1) == 0 ) && ( (input2&1) == 0 ) )
			  {
				if ( (output&1) != 0 )
				  {
					std::cout << "ERROR: " << (output>>1) << " should be TRUE" << std::endl
							  << "AND of: " << (input1>>1) << " and " << (input2>>1) << std::endl;
					error = true;
				  }
			  }
		  }
		// Or case
		else if ( _debugNetwork[foo].type == ORGATE )
		  {
			// input1 = TRUE OR input2 = TRUE => output = TRUE
			if ( (input1&1) == 0 )
			  {
				if ( (output&1) != 0 )
				  {
					std::cout << "ERROR: " << (output>>1) << " should be TRUE" << std::endl
							  << "OR of: " << (input1>>1) << " and " << (input2>>1) << std::endl;
					error = true;
				  }
			  }
			if ( (input2&1) == 0 )
			  {
				if ( (output&1) != 0 )
				  {
					std::cout << "ERROR: " << (output>>1) << " should be TRUE" << std::endl
							  << "OR of: " << (input1>>1) << " and " << (input2>>1) << std::endl;
					error = true;
				  }
			  }
			// input1 = FALSE AND input2 = FALSE => output = FALSE
			if ( ( (input1&1) != 0 ) && ( (input2&1) != 0 ) )
			  {
				if ( (output&1) == 0 )
				  {
					std::cout << "ERROR: " << (output>>1) << " should be FALSE" << std::endl
							  << "OR of: " << (input1>>1) << " and " << (input2>>1) << std::endl;
					error = true;
				  }
			  }
		  }
		// Halfand case
		else if ( _debugNetwork[foo].type == HALFANDGATE )
		  {
			// input1 = TRUE AND input2 = TRUE => output = TRUE
			if ( ( (input1&1) == 0 ) && ( (input2&1) == 0 ) )
			  {
				if ( (output&1) != 0 )
				  {
					std::cout << "ERROR: " << (output>>1) << " should be TRUE" << std::endl
							  << "HALFAND of: " << (input1>>1) << " and " << (input2>>1) << std::endl;
					error = true;
				  }
			  }
			// output = FALSE => input1 = FALSE OR input2 = FALSE
			if( (output&1) != 0 )
			  {
				if( ((input1&1) == 0) && ((input2&1) == 0) )
				  {
					std::cout << "ERROR: " << (input1>>1) << " or " << (input2>>1) << " should be FALSE" << std::endl
							  << "HALFAND of: " << (output>>1) << std::endl;
					error = true;
				  }
			  }
		  }
		// Halfor case
		else if ( _debugNetwork[foo].type == HALFORGATE )
		  {
			assert( input2 == 0 );
			// input1 = TRUE => output = TRUE
			if ( (input1&1) == 0 )
			  {
				if ( (output&1) != 0 )
				  {
					std::cout << "ERROR: " << (output>>1) << " should be TRUE" << std::endl
							  << "HALFOR of: " << (input1>>1) << " and " << (input2>>1) << std::endl;
					error = true;
				  }
			  }

			// output = FALSE => input1 = FALSE
			if ( (output&1) != 0 )
			  {
				if ( (input1&1) == 0 )
				  {
					std::cout << "ERROR: " << (input1>>1) << " should be FALSE" << std::endl
							  << "HALFOR of: " << (output>>1) << std::endl;
					error = true;
				  }
			  }
		  }
	  }
  
	if ( error )
	  {
		for( uint32_t v = 1; v != newModel.size(); ++v )
		  {
			std::cout << helper::Lit(newModel[v]) << " ";
		  }
		std::cout << std::endl;
	  }
  
	assert( !error);
  }
#endif

  void Sorter::Print(void) const
  {
	std::cout << "current sorter: " << std::endl;
	for( uint32_t i = 0; i != _softClauses.size(); ++i )
	  {
		std::cout << "clause: ";
		for( uint32_t j = 0 ; j != _softClauses[i]->clause.size(); ++j )
		  { 
			std::cout << helper::Lit( _softClauses[i]->clause[j]) << " "; 
		  }
		std::cout << "relaxLit: " << helper::Lit(_softClauses[i]->relaxationLit) << " weight: " << _softClauses[i]->weight << std::endl;
	  }

	std::cout << "outputs: " << std::endl;
	for( uint32_t i = 0; i != size(); ++i )
	  {
		std::cout << " " << i << " " << _outputs[i] << " [" << helper::Lit(_antom->Model()[_outputs[i]]) << "] ";
	  }
	std::cout << std::endl;
  }

  void Antom::DumpBucketModel(const std::vector<uint32_t>& model)
  {
	for (int32_t whichBucket = (static_cast<int32_t>(_sorterTree.size()) - 1); whichBucket >= 0; --whichBucket)
	  {
		std::cout << std::endl << "model of sorteroutputs[" << whichBucket << "]: (";

		assert( _sorterTree[whichBucket].size() == 1 );
		int32_t pos(-1);
		for (uint32_t whichIndex = 0; whichIndex < _sorterTree[whichBucket][0]->size(); ++whichIndex)
		  {
			uint32_t currentoutput(_sorterTree[whichBucket][0]->GetOutputs()[whichIndex]); 
			std::cout << helper::Lit(model[currentoutput]);
			if ( (pos == -1) && (model[currentoutput] == (currentoutput<<1)) )
			  {
				pos = (int32_t)whichIndex;
			  }
			if (whichIndex != _sorterTree[whichBucket][0]->size() - 1)
			  { std::cout << ", "; }
		  }

		std::cout << ")" << std::endl;

		if( pos == -1 )
		  { 
			std::cout << "all positions of " << _sorterTree[whichBucket][0]->size() << " are unsat => weight = " << (_sorterTree[whichBucket][0]->size()*static_cast<size_t>(std::pow(_antomSetting->base,whichBucket))) << std::endl;
		  }
		else
		  {
			std::cout << "sat at position: " << pos << " of " << _sorterTree[whichBucket][0]->size() << " => weight = " << (uint32_t)(pos*(std::pow(_antomSetting->base,whichBucket))) << std::endl;
		  }
		std::cout << "model of sorterinputs[" << whichBucket << "]: (";

		assert( _sorterTree[whichBucket].size() == 1 );
		for (uint32_t whichIndex = 0; whichIndex < _sorterTree[whichBucket][0]->GetSoftClauses().size(); ++whichIndex)
		  {
			uint32_t currentinput(_sorterTree[whichBucket][0]->GetSoftClauses()[whichIndex]->relaxationLit); 
			std::cout << helper::Lit(model[currentinput>>1]);
			if (whichIndex != _sorterTree[whichBucket][0]->GetSoftClauses().size() - 1)
			  { std::cout << ", "; }
		  }
		std::cout << ")" << std::endl;
		std::cout << "model of sorter tares[" << whichBucket << "]: (";
		for (uint32_t whichIndex = 0; whichIndex < _sorterTree[whichBucket][0]->GetTares().size(); ++whichIndex)
		  {
			std::cout << helper::Lit(model[_sorterTree[whichBucket][0]->GetTares()[whichIndex]]);
			if (whichIndex != _sorterTree[whichBucket][0]->GetTares().size() - 1)
			  { std::cout << ", "; }
		  }
		std::cout << ")" << std::endl;
	  }
	std::cout << std::endl;
  }
}
