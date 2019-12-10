/********************************************************************************************
bva.cpp -- Copyright (c) 2016, Sven Reimer

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

#include "bva.h"
#include "preprocessor.h"
#include "core.h"

namespace antom {

  /* BVA */
  // Based on "Automated Reencoding of Boolean formulas" Manthey, et al., 2012

  BVA::BVA(Preprocessor* prepro):
	_prepro(prepro),
	_setting(prepro->_setting),
	_queue(helper::DescendingOrder<size_t>(prepro->_occurCounts)),
	_seen(),
	_seenmlit(),
	_mlits(),
	_mlitscurrent(),
	_mcls(),
	_potential(),
	_clausesToRemove(),
	_touched(),
	_ca(prepro->_ca),
	_statistics(prepro->_statistics),
	_okay(prepro->_okay)
  {}


  uint32_t BVA::GetLeastOccuring(const WatcherFull& c)
  {
	// Mark mlits
	for( uint32_t i = 0; i != _mlits.size(); ++i )
	  {
		_seen[_mlits[i].lit1] = true;
		_seen[_mlits[i].lit2] = true;
	  }
	_seen[0] = false;

	uint32_t smallestLit(0);
	uint32_t leastOccuring(~0);

	_prepro->_methodLimit -= c.size();
	for( uint32_t i = 0; i != c.size(); ++i )
	  {
		if ( _seen[c[i]] )
		  { continue; }
		
		uint32_t litcost = _prepro->OccurenceCount(c[i]);
		if( litcost < leastOccuring )
		  {
				smallestLit = c[i];
				leastOccuring = litcost;
		  }
	  }

	// Unmark mlits
	for( uint32_t i = 0; i != _mlits.size(); ++i )
	  {
		_seen[_mlits[i].lit1] = false;
	  }
	_seen[0] = false;

	return smallestLit;
  }

  BVA::LitPair BVA::GetMostOccuringPotential(uint32_t& largest)
  {
	largest = 0;

	// Sort by "lit" in PotentialClause
	std::sort(_potential.begin(), _potential.end());

	_prepro->_methodLimit -= _potential.size();

	LitPair most(0);
	LitPair last(0);
	uint32_t num(0);

	for( uint32_t i = 0; i != _potential.size(); ++i )
	  {
		if ( last != _potential[i].litpair )
		  {
			if( num >= largest )
			  {
				largest = num;
				most = last;
			  }
			last = _potential[i].litpair;
			num = 1;
		  }
		else
		  {
			++num;
		  }
	  }

	if( num >= largest )
	  {
		largest = num;
		most = last;
	  }

	return most;
  }

  void BVA::RemoveDuplicates(void)
  {
	if ( _mcls.size() == 1 )
	  { return; }

	_prepro->_methodLimit -= 2*_mcls.size();
	std::sort( _mcls.begin(), _mcls.end(), WatchedSorterFull() );

	uint32_t i(0);
	uint32_t j(0);
	for ( ;i+1 < _mcls.size(); ++i )
	  {
		const WatcherFull& prev = _mcls[j];
		const WatcherFull& next = _mcls[i+1];
		if( prev.size() != next.size() )
		  {
			_mcls[j+1] = _mcls[i+1];
			++j;
			continue;
		  }

		bool toDelete(true);
		
		_prepro->_methodLimit -= prev.size();
		for ( uint32_t i = 0; i != prev.size(); ++i )
		  {
			if( prev[i] != next[i] )
			  {
				toDelete = false;
				break;
			  }
		  }

		if ( !toDelete )
		  {
			_mcls[j+1] = _mcls[i+1];
			++j;
		  }
	  }
	_mcls.resize(_mcls.size() - (i-j) );
  }

  BVA::LitPair BVA::DiffClauses(const WatcherFull& c1, const WatcherFull& c2) 
  {
	uint32_t differences(0);
	LitPair difflit(0);

	uint32_t c2size(c2.size());
	for( uint32_t i = 0; i != c2.size(); ++i )
	  {
		_seen[c2[i]] = true;
	  }

	uint32_t c1size(c1.size());
	_prepro->_methodLimit -= (c2size+c1size);

	for( uint32_t i = 0; i != c1.size(); ++i )
	  {
		if( !_seen[c1[i]] )
		  {
			if( differences == 0 )
			  {
				difflit.lit1 = c1[i];
			  }
			else
			  {
				difflit.lit2 = c1[i];
			  }
			++differences;
			if( differences > 2 )
			  { break; }
		  }
	  }

	for( uint32_t i = 0; i != c2size; ++i )
	  {
		_seen[c2[i]] = false;
	  }

  if( (differences >= 1) && (differences <= 2 ) )
	{ 
	  return difflit;
	}
  else
	{
	  return LitPair(0);
	}
  }

  bool BVA::SmallerFormula(uint32_t occurences) const
  {
	// In first run, we need at least two occurences
	if( _mlits.size() == 1 )
	  {
		return (occurences >=2);
	  }

    int32_t orgsize = FormulaSize( (int32_t)_mlits.size(), (int32_t)_mcls.size() );
    int32_t newsize = FormulaSize( (int32_t)_mlits.size()+1, (int32_t)occurences );

	if( newsize <= 0 )
	  { return false; }

	if( newsize < orgsize )
	  { return false; }

	return true;
  }

  int32_t BVA::FormulaSize(int32_t litsize, int32_t clssize) const
  {
	return (litsize*clssize-litsize-clssize);
  }

  uint32_t BVA::GetNewVariable(void)
  {
	uint32_t index = _prepro->_variables+1;
	++_statistics.usedVariables;
	_prepro->_core->SetMaxIndex(index);
	// Insert new variable for core
	_prepro->_core->_varOrder[_prepro->_core->_varGroup[index]]->resize(index);
	_prepro->_core->_varOrder[_prepro->_core->_varGroup[index]]->insert(index);
	_prepro->UpdateDataStructures();
	_seen.resize((index+1)<<1,false);
	_seenmlit.resize((index+1)<<1,false);

	return index;
  }

  bool BVA::PerformBVA(void) 
  {
	_touched.clear();
	uint32_t x = GetNewVariable();
	_queue.resize((x+1)<<1);
	++_statistics.bvaVariables;
	std::vector< uint32_t > clause;

	size_t addedlits(0);
	size_t removedlits(0);
	// add one resolvent for mlit
	for ( uint32_t i = 0; i != _mlits.size(); ++i )
	  {
		clause.clear();
		assert( _mlits[i].lit1 != 0 );
		clause.push_back(_mlits[i].lit1);
		
		if( _mlits[i].lit2 != 0 )
		  {
			clause.push_back( _mlits[i].lit2 );
			// Ensure that clause is sorted
			if( clause[1] < clause[0] )
			  {
				std::swap(clause[0],clause[1]);
			  }
		  }
		clause.push_back(x<<1);
		addedlits += clause.size();
		if( !_prepro->AddClausePrepro(clause) )
		  { return false; }
		_touched.Touch(clause);
	  }


	// add other resolvent for mcls
	for ( uint32_t i = 0; i != _mcls.size(); ++i )
	  {
		clause.clear();
		const WatcherFull& c = _mcls[i];

		uint32_t wsize(c.size());
		clause.resize(wsize);

		uint32_t j = 0;
		for( uint32_t i = 0; i != wsize; ++i )
		  {
			if( c[i] != c.wlit )
			  {
				clause[j] = c[i];
				++j;
			  }
		  }
		// Ensure that the new variable is added back, so that clause is automatically sorted
		clause[j] = (x<<1)^1;
		addedlits += clause.size();
		if ( !_prepro->AddClausePrepro(clause) )
		  { return false; }
		_touched.Touch(clause);
	  }

	std::sort(_prepro->_binaries[(x<<1)].begin(),_prepro->_binaries[(x<<1)].end(), WatchedSorter() );
	std::sort(_prepro->_binaries[(x<<1)^1].begin(),_prepro->_binaries[(x<<1)^1].end(), WatchedSorter() );

	// Now remove old clauses
    CollectClausesForRemove();

	std::vector< uint32_t > toRemove;
	for ( uint32_t j = 0; j != _mlits.size(); ++j )
	  {
		uint32_t removelit( _mlits[j].lit1 );
		uint32_t removelit2( _mlits[j].lit2 );
		for ( uint32_t i = 0; i != _clausesToRemove.size(); ++i )
		  {
			toRemove.clear();
			toRemove.push_back(removelit);
			if ( removelit2 != 0 )
			  {
				toRemove.push_back(removelit2);
			  }
			for( uint32_t k = 0; k != _clausesToRemove[i].size(); ++k )
			  {
				toRemove.push_back(_clausesToRemove[i][k]);
			  }

			removedlits += toRemove.size();
			if( toRemove.size() == 2 )
			  {
				_prepro->RemoveBinary(toRemove[0], toRemove[1]);
				_prepro->RemoveBinary(toRemove[1], toRemove[0]);
				_statistics.currentBinaryClauses -= 2;
			  }
			else 
			  {
				CRef c = GetBVAClause( toRemove );
				assert( !_ca[c].IsLearned() );
				_prepro->ClearAllOccurences(c);
			  }
		  }
	  }

	assert( addedlits < removedlits );
	_statistics.bvaLiterals += static_cast<uint32_t>(removedlits-addedlits);

	UpdateQueue();

	return true;
  }

  void BVA::CollectPotentialClauses(uint32_t nextlit)
  {
	_potential.clear();
	LitPair nextpair(nextlit);
	// Check each clause for potential clauses
	for ( uint32_t i = 0; i != _mcls.size(); ++i )
	  {
		const WatcherFull& c( _mcls[i] );
		const uint32_t cSize( c.size() );

		uint32_t lmin = GetLeastOccuring( c );

		if( lmin == 0 )
		  { continue; }

		_mlitscurrent = _mlits;
		_prepro->_methodLimit -= _mlitscurrent.size();
		// Mark current lits as already seen
		for( uint32_t l = 0; l != _mlitscurrent.size(); ++l )
		  {
			_seenmlit[_mlitscurrent[l].hash(_seenmlit.size())] = true;
		  }
		
		_prepro->_methodLimit -= _prepro->OccurenceCount(lmin,false);
		for( uint32_t j = 0; j != _prepro->_binaries[lmin].size(); ++j )
		  {
			if( _prepro->_binaries[lmin][j].IsLearnedBinary() )
			  { continue; }
			const WatcherFull& d = WatcherFull( lmin, _prepro->_binaries[lmin][j], _ca );

			const uint32_t dSize( d.size() );

			if( ( c != d ) && 
				( (cSize == dSize) || ( (cSize+1 == dSize) && _setting->bvaTwoLitDiff ) ) && 
				( DiffClauses(c,d) == nextpair ) )
			  {

				if( (cSize+1) == dSize )
				  { 
					std::cout << "two lit diff " << std::endl;
					c.Print();
					d.Print();
				  }

				const LitPair difflits = DiffClauses(d,c);
				if( !_seenmlit[difflits.hash(_seenmlit.size())] )
				  {
					_potential.push_back(PotentialClause(difflits, c ) );
					_seenmlit[difflits.hash(_seenmlit.size())] = true;
					_mlitscurrent.push_back(difflits);
				  }
			  }
		  }

		for( uint32_t j = 0; j != _prepro->_occur[lmin].size(); ++j )
		  {
			Clause& clause = _ca[_prepro->_occur[lmin][j]];
			if( clause.IsLearned() )
			  { continue; }
			const WatcherFull& d = WatcherFull( lmin, Watcher(_prepro->_occur[lmin][j], clause[0], NNARY ), _ca );
			const uint32_t dSize( d.size() );

			if( ( c != d ) && 
				( (cSize == dSize) || ( (cSize+1 == dSize) && _setting->bvaTwoLitDiff ) ) && 
				( DiffClauses(c,d) == nextpair ) )
			  {

				if( (cSize+1) == dSize )
				  { 
					std::cout << "two lit diff" << std::endl;
					c.Print();
					d.Print();
				  }
				LitPair difflits = DiffClauses(d,c);

				if( !_seenmlit[difflits.hash(_seenmlit.size())] )
				  {
					_potential.push_back(PotentialClause(difflits, c ) );
					_seenmlit[difflits.hash(_seenmlit.size())] = true;
					_mlitscurrent.push_back(difflits);
				  }
			  }
		  }

		for( uint32_t l = 0; l != _mlitscurrent.size(); ++l )
		  {
			_seenmlit[_mlitscurrent[l].hash(_seenmlit.size())] = false;
		  }
	  }
  }

  void BVA::UpdateQueue(void)
  {
	const std::vector< uint32_t >& touchlist = _touched.GetList();
	for ( uint32_t i = 0; i != touchlist.size(); ++i )
	  {
		uint32_t var(touchlist[i]);
		if ( _queue.inHeap(var<<1) )
		  {
			_prepro->OccurenceCount(var<<1, false);			
			_queue.update(var<<1);
		  }
		if ( _queue.inHeap((var<<1)^1) )
		  {
			_prepro->OccurenceCount((var<<1)^1, false);			
			_queue.update((var<<1)^1);
		  }
	  }
	_touched.clear();
  }

  void BVA::CollectClausesForRemove(void)
  {
	_clausesToRemove.clear();
	std::vector< uint32_t > clause;
	
	for( uint32_t i = 0; i != _mcls.size(); ++i )
	  {
		clause.clear();
		
		for( uint32_t j = 0; j != _mcls[i].size(); ++j )
		  {
			if( _mcls[i].wlit != _mcls[i][j] )
			  {			
				clause.push_back( _mcls[i][j] );
			  }
		  }
		_clausesToRemove.push_back(clause);
	  }
  }

  CRef BVA::GetBVAClause(const std::vector< uint32_t >& clause)
  {
	CRef ref(0);
	uint32_t minoccur(0);
	size_t numberofminoccur(std::numeric_limits<size_t>::max());
	for( uint32_t i = 0; i != clause.size(); ++i )
	  {
		_seen[clause[i]] = true;
		if( numberofminoccur > _prepro->_occur[clause[i]].size() )
		  {
			minoccur = clause[i];
			numberofminoccur = _prepro->_occur[clause[i]].size();
		  }
	  }

	assert( minoccur != 0 );

	uint32_t i = 0;
	assert(_prepro->_occur[minoccur].size() != 0);
	for( ; i != _prepro->_occur[minoccur].size(); ++i )
	  {
		ref = _prepro->_occur[minoccur][i];
		const Clause& check = _ca[ref];
		
		uint32_t j = 0;
		for( ; j != check.size(); ++j )
		  {
			if( !_seen[check[j]] )
			  {
				break;
			  }
		  }

		if ( j == check.size() )
		  {
			break;
		  }
	  }

	assert( i != _prepro->_occur[minoccur].size() );

	for( uint32_t i = 0; i != clause.size(); ++i )
	  {
		_seen[clause[i]] = false;
	  }
	return ref;
  }

  bool BVA::BoundedVariableAddition(void)
  {
	if( _setting->verbosity > 2 )
	  { std::cout << "c Performing bounded variable addition..." << std::endl; }

	struct rusage resources;
	getrusage(RUSAGE_SELF, &resources);
	double timeS = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	uint32_t bvavars = _statistics.bvaVariables;
	uint32_t bvalits = _statistics.bvaLiterals;
	_okay = true;

	_prepro->_methodLimit = (1<<25);

	_queue.clear();
	_queue.resize((_prepro->_variables+1)<<1);
	_seen.resize((_prepro->_variables+1)<<1,false);
	_seenmlit.resize((_prepro->_variables+1)<<1,false);

	_prepro->CountOccurences(false);

	// Collect and sort candidates
	for ( uint32_t v = 1; v <= _prepro->_variables; ++v )
	  {
		uint32_t poslit(v<<1);
		uint32_t neglit((v<<1)^1);
		if( !_prepro->IsUsed(v) )
		  { continue; }
		_queue.insert(poslit);
		_queue.insert(neglit);
	  }

	while ( !_queue.empty() )
	  {
		if( _prepro->_methodLimit < 0 )
		  { break; }

		uint32_t nextlit(_queue.top());

		// No occurence? Then we can stop here
		if ( _prepro->_occurCounts[nextlit] == 0 )
		  { break; }

		//		std::cout << "nextlit: " << Lit(nextlit) << " " << _prepro->_occurCounts[nextlit] << std::endl;
		_mcls.clear();
		_mlits.clear();
		_mlits.push_back(LitPair(nextlit));
		

		// Fill current clause lists
		for ( uint32_t i = 0; i != _prepro->_binaries[nextlit].size(); ++i )
		  {
			// Consider static clauses only
			if( !_prepro->_binaries[nextlit][i].IsLearnedBinary() )
			  {
				_mcls.push_back(WatcherFull(nextlit, _prepro->_binaries[nextlit][i], _ca ) );
			  }
		  }
		for ( uint32_t i = 0; i != _prepro->_occur[nextlit].size(); ++i )
		  {
			// Consider static clauses only
			if( !_ca[_prepro->_occur[nextlit][i]].IsLearned() )
			  {
				_mcls.push_back(WatcherFull(nextlit, Watcher(_prepro->_occur[nextlit][i],nextlit,NNARY), _ca) );
			  }
		  }

		RemoveDuplicates();

		_prepro->_methodLimit -= _mcls.size();
		
		while ( true )
		  {
			if ( _prepro->_methodLimit < 0 )
			  { break; }
			CollectPotentialClauses(nextlit);
			
			uint32_t numoccur(0);
			LitPair lmax = GetMostOccuringPotential(numoccur);

			//std::cout << "most occuring: " << Lit(lmax.lit1) << " " << Lit(lmax.lit2) << std::endl;

			// TODO, does lmax^1 == nextlit happens? If so, we can perform selfsubsumption

			// Formula shrinks by adding "lmax"?
			if ( SmallerFormula(numoccur) )
			  {
				_mlits.push_back(lmax);
				_mcls.clear();

				_prepro->_methodLimit -= _potential.size();
				// update "_mcls" with potential clauses
				for ( uint32_t j = 0; j != _potential.size(); ++j )
				  {
					if( _potential[j].litpair == lmax )
					  {
						_mcls.push_back(_potential[j].clause);
					  }
				  }
			  }
			else
			  { break; }
		  }

		// We found our candidate for BVA
		if ( FormulaSize((int32_t)_mlits.size(), (int32_t)_mcls.size()) > 0 )
		  {
			_prepro->SetFlag(PF_BVA);

			if ( !PerformBVA() )
			  { 
				_okay = false;
				break;
			  }
		  }
	  }

	getrusage(RUSAGE_SELF, &resources);
	double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
	
	_statistics.runtime_bva += (timeC-timeS);

	if ( _setting->verbosity > 1 )
	  {
		std::cout << "c runtime bva            : " << (timeC-timeS) << "s (" << _statistics.runtime_bva << "s)" << std::endl;
		if ( _statistics.bvaVariables != bvavars )
		  {
			std::cout << "c bvas                   : " << (_statistics.bvaVariables-bvavars) << std::endl;
			std::cout << "c literal reduction      : " << (_statistics.bvaLiterals-bvalits) << std::endl;

		  }
	  }
	return _okay;
  }
}
