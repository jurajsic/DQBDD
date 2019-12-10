/********************************************************************************************
hte.cpp -- Copyright (c) 2016, Sven Reimer

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

#include "hte.h"
#include "preprocessor.h"
#include "watcher.h"

namespace antom {

  /* HTE */
  // Hidden Clause Elimination based on Heule, et al., 2010

  HTE::HTE(Preprocessor* prepro) :
	_setting(prepro->_setting),
	_prepro(prepro),
	_statistics(prepro->_statistics),
	_okay(true),
	_candidates(helper::AscendingOrder<size_t>(prepro->_occurCounts)),
	_seen()
  {}

  void HTE::DoHiddenTautologyClauseElimination(void)
  {
	if( _setting->verbosity > 2 )
	  { std::cout << "c Performing HTE..." << std::endl; }

	_okay = true;

	struct rusage resources;
	getrusage(RUSAGE_SELF, &resources);
	double timeS = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	uint32_t htaut = _statistics.hiddenTautologies;
	uint32_t hsub = _statistics.hiddenSubsumptions;
	
	_seen.resize((_prepro->_variables+1)<<1,false);
	_candidates.clear();
	_candidates.resize((_prepro->_variables+1)<<1);
	
	_prepro->CountOccurences(true);

	_prepro->_methodLimit = (1<<22);
	
	for ( uint32_t v = 1; v <= _prepro->_variables; ++v )
	  {
		// Skip don't touch, deleted and assigned variables
		if ( !_prepro->IsUsed(v) || _prepro->_donttouch[v])
		  { continue; }
		
		_candidates.insert(v<<1);
		_candidates.insert((v<<1)^1);
	  }

	std::vector< uint32_t> clause;
	clause.resize(2,0);
	// Check for all binaries
	while ( !_candidates.empty() )
	  {
		uint32_t lit = _candidates.top();

		// std::cout << __func__ << " " << Lit(lit) << std::endl;
		// _prepro->printWatchList(lit);

		std::vector< Watcher >& watches ( _prepro->_binaries[lit] );
		size_t size( watches.size() );
		_prepro->_methodLimit -= size;
		if ( _prepro->_methodLimit < 0 )
		  {
			_okay = false;
			break;
		  }
		for ( uint32_t j = 0; j != size; ++j )
		  {
			//std::cout << j << " " << Lit(watches[j].getSecondLit()) << std::endl;
			clause.resize(2);
			clause[0] = lit;
			uint32_t secLit(watches[j].GetSecondLit()); 
			// Skip opposite binary clause
			if ( lit < secLit )
			  {
				//uint64_t sign(0ULL);
				uint32_t sign(0);
				clause[1] = secLit;
				_seen[lit] = true;
				_seen[secLit] = true;

				if ( CheckHiddenTautology(clause, sign) )
				  {
					if( clause.size() == 2 )
					  {
						for( uint32_t foo = 0; foo != clause.size(); ++foo)
						  {
							std::cout << helper::Lit(clause[foo]) << " ";
						  }
						std::cout << std::endl;
					  }
					// We should added something to the clause...
					assert( clause.size() > 2 );
					// Remove the hidden tautology
					_prepro->RemoveBinary( lit, secLit );
					_prepro->RemoveBinary( secLit, lit );

					_statistics.currentBinaryClauses -= 2;
					if ( watches[j].IsLearnedBinary() )
					  {
						_statistics.learnedBinary -= 2;
					  }
					--j; --size;
					_prepro->SetFlag(PF_HTE);
				  }
				// check for hidden subsumptions if clause is not a tautology
				else if (_prepro->_setting->doHse)
				  {
					// clause must be sorted for subsumption checks
					std::sort(clause.begin(), clause.end());
					if (_prepro->IsForwardSubsumedExcept(clause, sign, lit, secLit))
					  {
						// Remove the hidden subsumption
						_prepro->RemoveBinary( lit, secLit );
						_prepro->RemoveBinary( secLit, lit );

						_statistics.currentBinaryClauses -= 2;
						if ( watches[j].IsLearnedBinary() )
						  {
							_statistics.learnedBinary -= 2;
						  }
						--j; --size;

						++_statistics.hiddenSubsumptions;

						_prepro->SetFlag(PF_HTE);
					  }
				  }

				for ( uint32_t k = 0; k != clause.size(); ++k )
				  {
					_seen[clause[k]] = false;
				  }
			  }
		  }
	  }

	if (_okay)
	  {
		// Check for all non-binaries
		size_t dbsize(_prepro->_clauseDatabase.size());
		for ( uint32_t i = 0; i != dbsize; ++i )
		  {
			Clause& c( _prepro->_ca[_prepro->_clauseDatabase[i]] );
			//uint64_t sign(0ULL);
			uint32_t sign(0);
			clause.resize(c.size());

			_prepro->_methodLimit -= c.size();
			if ( _prepro->_methodLimit < 0 )
			  {
				_okay = false;
				break;
			  }

			for ( uint32_t j = 0; j != c.size(); ++j )
			  {
				clause[j] = c[j];
				_seen[clause[j]] = true;
			  }
	
			if( CheckHiddenTautology ( clause, sign ) )
			  {
				// Delete clause and its occurences
				_prepro->ClearAllOccurences( _prepro->_clauseDatabase[i] );

				_prepro->_ca.Free(_prepro->_clauseDatabase[i]);
				_prepro->_clauseDatabase[i--] = _prepro->_clauseDatabase[--dbsize];
				_prepro->_clauseDatabase.pop_back();

				_prepro->SetFlag(PF_HTE);
			  }
			// check for hidden subsumptions if clause is not a tautology
			else if (_prepro->_setting->doHse)
			  {
				// clause must be sorted for subsumption checks
				std::sort(clause.begin(), clause.end());
				if (_prepro->IsForwardSubsumedExcept(clause, sign, _prepro->_clauseDatabase[i]))
				  {
					// Delete clause and its occurences
					_prepro->ClearAllOccurences( _prepro->_clauseDatabase[i] );

					_prepro->_ca.Free(_prepro->_clauseDatabase[i]);
					_prepro->_clauseDatabase[i--] = _prepro->_clauseDatabase[--dbsize];
					_prepro->_clauseDatabase.pop_back();

					++_statistics.hiddenSubsumptions;

					_prepro->SetFlag(PF_HTE);
				  }
			  }

			for ( uint32_t j = 0; j != clause.size(); ++j )
			  {
				_seen[clause[j]] = false;
			  }
		  }
	  }

	getrusage(RUSAGE_SELF, &resources);
	double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	_statistics.runtime_hte += (timeC-timeS);

	if( _setting->verbosity > 1 )
	  {
		if (_statistics.hiddenTautologies != htaut)
		  {
			std::cout << "c hidden tautologies     : " << ( _statistics.hiddenTautologies - htaut ) << std::endl ;
		  }
		if (_statistics.hiddenSubsumptions != hsub)
		  {
			std::cout << "c hidden subsumptions    : " << ( _statistics.hiddenSubsumptions - hsub ) << std::endl ;
		  }
		std::cout << "c runtime hte            : " << (timeC-timeS) << "s (" << _statistics.runtime_hte << "s)" << std::endl;
	  }
  }

 
  //  bool HTE::CheckHiddenTautology(std::vector< uint32_t >& clause, uint64_t& sign)
  bool HTE::CheckHiddenTautology(std::vector< uint32_t >& clause, uint32_t& sign)
  {

	size_t csize(clause.size());
	bool binary(false);
	if (csize == 2)
	  {
		binary = true;
	  }

	assert( clause[0] < clause[1] );	

	for (uint32_t i = 0; i != csize; ++i)
	  {
		uint32_t l0 = clause[i];

		//sign |= 1ULL << (l0 % 64); 
		sign |= 1 << (l0 % 32); 

		_prepro->_methodLimit -= _prepro->_binaries[l0].size();

		for (uint32_t j = 0; j != _prepro->_binaries[l0].size(); ++j)
		  {
			uint32_t blit( _prepro->_binaries[l0][j].GetSecondLit() );
			
			// Skip currently considered binary clause and learned binaries
			if (_prepro->_binaries[l0][j].IsLearnedBinary() || 
				 ( binary && ( ( (blit == clause[0]) && (l0 == clause[1]) )
							 || ( (blit == clause[1]) && (l0 == clause[0]) ) ) 
				  ))
			  {
				continue;
			  }

			// We found a hidden tautology!
			if (_seen[blit])
			  {
				++_statistics.hiddenTautologies;
				return true;
			  }
			else if (!_seen[blit^1])
			  {
				clause.push_back( blit^1 );
				++csize;
				_seen[blit^1] = true;
			  }
		  }
	  }
	return false;
  }
}
