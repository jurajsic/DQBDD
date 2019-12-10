/********************************************************************************************
bce.cpp -- Copyright (c) 2016, Sven Reimer

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

#include "bce.h"
#include "preprocessor.h"
#include "core.h"
#include "modelrebuilder.h"

namespace antom {

  /* BCE */
  // Blocked Clause elimination based on Järvisalo, et al., 2010

  BCE::BCE(Preprocessor* prepro):
	_prepro(prepro),
	_setting(prepro->_setting),
	_rebuilder(prepro->_rebuilder),
	_ca(prepro->_ca),
	_model(prepro->_model),
	_statistics(prepro->_statistics),
	_candidates(helper::AscendingOrder<size_t>(prepro->_occurCounts))
  {}
  
  bool BCE::DoBlockedClauseElimination(void)
  {

	if( _setting->verbosity > 2 )
	  { std::cout << "c Performing BCE..." << std::endl; }

	struct rusage resources;
	getrusage(RUSAGE_SELF, &resources);
	double timeS = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	_candidates.clear();
	_candidates.resize((_prepro->_variables+1)<<1);
	_prepro->CountOccurences(false);

	for (uint32_t v = 1; v <= _prepro->_variables; ++v)
	  {
		// Skip don't touch, deleted and assigned variables
		if( !_prepro->_donttouch[v] && _prepro->IsUsed(v) )
		  {
			_candidates.insert(v<<1);
			_candidates.insert((v<<1)^1);
		  }
	  }
	uint32_t bclauses = _statistics.blockedClauses;

	while ( !_candidates.empty() )
	  {
		uint32_t checklit(_candidates.top());
		
		//		std::cout << "checklit: " << helper::Lit(checklit) << " count: " << _prepro->_occurCounts[checklit] << std::endl;
		//		_prepro->PrintCompleteLists(checklit>>1);
		if ( !CheckBlockedLit(checklit) )
		  {
			return false;
		  }
	  }

	_prepro->ClearClauseDatabase();

	getrusage(RUSAGE_SELF, &resources);
	double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	_statistics.runtime_bce += (timeC-timeS);
	
	if( _setting->verbosity > 1 )
	  {
		if (_statistics.blockedClauses != bclauses)
		  {
			std::cout << "c blocked clauses        : " << ( _statistics.blockedClauses - bclauses ) << std::endl ;
		  }
		std::cout << "c runtime bce            : " << (timeC-timeS) << "s (" << _statistics.runtime_bce << "s)" << std::endl;
	  }
	return _prepro->PropagateUnits();
  }

  void BCE::AddBlockedImplications(void)
  {
	if( _setting->verbosity > 2 )
	  { std::cout << "c Performing ABI..." << std::endl; }

	struct rusage resources;
	getrusage(RUSAGE_SELF, &resources);
	double timeS = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	uint32_t count(0);
	for (unsigned lit1 = 2; lit1 <= _prepro->_variables<<1; ++lit1)
	  {
		if (_prepro->_donttouch[lit1>>1] || !_prepro->IsUsed(lit1>>1))
		  {
			continue;
		  }
	
		for (unsigned lit2 = 2; lit2 <= _prepro->_variables<<1; ++lit2)
		  {
			if (!_prepro->IsUsed(lit2) || (lit1==lit2))
			  {
				continue;
			  }

			if( _prepro->HasBinary(lit1, lit2) )
			  {
				continue;
			  }
			
			if (CheckBlocked(lit1, lit2))
			  {
				_prepro->AddBinary(lit1,lit2,true);
				std::cout << "new blocked implication: " << helper::Lit(lit1) << " " << helper::Lit(lit2) << std::endl;
				++count;
			  }
		  }
	  }

	getrusage(RUSAGE_SELF, &resources);
	double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	if( _setting->verbosity > 1 )
	  {
		if (count != 0)
		  {
			std::cout << "c added blocked binaries : " << count << std::endl ;
			std::cout << "c runtime                : " << (timeC-timeS) << "s (" << _statistics.runtime_bce << "s)" << std::endl;
		  }
	  }
  }

  bool BCE::CheckBlockedLit(uint32_t lit)
  {
	const std::vector< Watcher >& binaries( _prepro->_binaries[lit] ); 
	const std::vector< CRef >& occur( _prepro->_occur[lit] );
	size_t bsize(binaries.size());
	size_t osize(occur.size());
	
	// Found monotone lit
	if ( (bsize+osize) == 0 )
	  { 
		// If literal is don't touch, do nothing, 
		if ( !_prepro->_donttouch[lit>>1] )
		  { 
			
			// otherwise we have detected a monotone literal
			if ( ( _prepro->_binaries[lit^1].size() + _prepro->_occur[lit^1].size() ) == 0 )
			  {
				// literal even don't care?
				_prepro->_deleted[lit>>1] = true;
				++_statistics.dontcareVariables;
				return true;
			  }
		
			_prepro->_core->AddImplication( lit^1 ); 
			++_statistics.monotoneVariables;
			_prepro->SetFlag(PF_BCE);

#ifndef NDEBUG
			bool res = _prepro->PropagateUnits(); assert(res);
#else
			_prepro->PropagateUnits();
#endif
		  }
		return true; 
	  }

	for (uint32_t i = 0; i != bsize; ++i )
	  {
		uint32_t seclit( binaries[i].GetSecondLit() );

		if ( CheckBlocked( lit, seclit ) )
		  { 
			if( !_prepro->_donttouch[seclit>>1] && !_candidates.inHeap(seclit^1) )
			  {
				_candidates.insert(seclit^1);
			  }
			RemoveBlocked(lit, binaries[i] );
			++_statistics.blockedClauses;
			--i; --bsize;
			_prepro->SetFlag(PF_BCE);
		  }
	  }

	for (uint32_t i = 0; i != osize; ++i )
	  {
		const Clause& clause = _ca[occur[i]];
		if ( CheckBlocked( lit, clause ) )
		  { 
			for (uint32_t j = 0; j != clause.size(); ++j)
			  {
				if( !_prepro->_donttouch[clause[i]>>1] && !_candidates.inHeap(clause[i]^1) )
				  {
					_candidates.insert(clause[i]^1);
				  }
			  }
			RemoveBlocked(lit, occur[i] );
			++_statistics.blockedClauses;
			--i; --osize;
			_prepro->SetFlag(PF_BCE);
		  }
	  }
	return true;
  }

  bool BCE::CheckBlocked(uint32_t lit, uint32_t secondlit) const
  {
	const std::vector< Watcher >& binaries( _prepro->_binaries[lit^1] ); 
	const std::vector< CRef >& occur( _prepro->_occur[lit^1] );
	size_t bsize(binaries.size());
	size_t osize(occur.size());

	if ( (bsize+osize) == 0 )
	  { return false; }

	uint32_t checklit = secondlit^1;
	uint32_t i = 0;

	for ( ; i != bsize; ++i )
	  {
		if ( binaries[i].GetSecondLit() != checklit )
		  { break; }
	  }

	if ( i != bsize)
	  {
		return false;
	  }


	i = 0;
	for ( ; i != osize; ++i )
	  {
		if ( _ca[occur[i]].ToDelete() )
		  { continue; }

		if ( _prepro->CountMergeClauses( _ca[occur[i]], lit, secondlit ) != 0 )
		  { break; }
	  }

	if ( i != osize )
	  { 
		return false;
	  }

	return true;

  }

  bool BCE::CheckBlocked(uint32_t lit, const Clause& clause) const
  {
	if( clause.ToDelete() )
	  { return false; }

	const std::vector< Watcher >& binaries( _prepro->_binaries[lit^1] ); 
	const std::vector< CRef >& occur( _prepro->_occur[lit^1] );
	size_t bsize(binaries.size());
	size_t osize(occur.size());

	if ( (bsize+osize) == 0 )
	  { return false; }

	uint32_t i = 0;
	for ( ; i != bsize; ++i )
	  {
		if ( _prepro->CountMergeClauses( clause, lit, binaries[i].GetSecondLit() ) != 0 )
		  { break; }
	  }

	if ( i != bsize)
	  {
		return false;
	  }


	i = 0;
	for ( ; i != osize; ++i )
	  {
		
		if ( _prepro->CountMergeClauses( clause, _ca[occur[i]], lit ) != 0 )
		  { break; }
	  }

	if ( i != osize )
	  { 
		return false;
	  }

	return true;
  }

  void BCE::RemoveBlocked(uint32_t bLit, const Watcher& binClause)
  {
	assert( binClause.IsBinary() );
	uint32_t secLit( binClause.GetSecondLit() );
	std::vector< uint32_t > blockedClause {bLit, secLit};

	_rebuilder->AddBlockedClause( bLit, blockedClause );
	_prepro->RemoveBinary(bLit,secLit);
	_prepro->RemoveBinary(secLit,bLit);

	_statistics.currentBinaryClauses -= 2;
	if (binClause.IsLearnedBinary() )
	  {
		_statistics.learnedBinary -=2;
	  }
	  
  }

  void BCE::RemoveBlocked(uint32_t bLit, CRef clause)
  {
	const Clause& bclause(_ca[clause]);
	assert(!bclause.ToDelete());
	std::vector< uint32_t > blockedClause(bclause.size(),0);
	for( uint32_t i = 0; i != bclause.size(); ++i )
	  {
		blockedClause[i] = bclause[i];
	  }
	  
	_rebuilder->AddBlockedClause( bLit, blockedClause );

	_prepro->ClearAllOccurences(clause);
  }

}
