/********************************************************************************************
upla.cpp -- Copyright (c) 2014-2016, Sven Reimer

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

#include "upla.h"

#include "preprocessor.h"
#include "core.h"
#include "modelrebuilder.h"

namespace antom {

  UPLA::UPLA(Preprocessor* prepro):
	_core(prepro->_core),
	_prepro(prepro),
	_setting(prepro->_setting),
	_rebuilder(prepro->_rebuilder),
	_nextcandidates(),
	_newEquivalences(),
	_posimplications(),
	_negimplications(),
	_activity(),
	_heap(helper::DescendingOrder<uint32_t>(_activity)),
	_dsEndIndex(prepro->_dsEndIndex),
	_dsImplIndex(prepro->_dsImplIndex),
	_assignment(prepro->_assignment),
	_ca(prepro->_ca),
	_statistics(prepro->_statistics),
	_okay(prepro->_okay)
  {
	_nextcandidates.resize(prepro->_variables+1,false);
  }

  // Backtracks to decision literal "declit". 
  void UPLA::Backtrack(uint32_t declit)
  {
	// Initialization.
	uint32_t lit(_prepro->_decisionStack[--_dsEndIndex]);
	
	// This loop assumes that we have a dummy assignment at the first 
	// position of the decision stack, which has been assigned on decision level 0. 
	while( lit != declit )
	  {
		// Undo the current assignment.
		_assignment[lit] = false;

		_prepro->_forcing[lit>>1].ClearReason();
		  
		// Get the next assignment to be undone.
		lit = _prepro->_decisionStack[--_dsEndIndex];
	  } 

	// Undo the assignment of "declit"
	_assignment[lit] = false;

	assert(_prepro->_forcing[lit>>1].NoReason());
  
	// Update "_dsImplIndex".
	_dsImplIndex = _dsEndIndex;
	return;
  } 


  // Adds a decision/implication to the temp upla decision stack.
  void UPLA::AddUnit(uint32_t lit)
  {    
	// The variable corresponding to "lit" has to be undefined.
	assert(!_assignment[lit] && !_assignment[lit ^ 1]); 
	assert( lit > 1 );

	// Update "_assignment".
	_assignment[lit] = true;

	_prepro->_forcing[lit>>1] = Reason();
    
	// Push "lit" onto the decision stack.    
	_prepro->_decisionStack[_dsEndIndex++] = lit; 
  }

  /* End: some helper fucntions */

  // Performs "Unit Propagation Lookahead". 
  // Returns FALSE if the formula is unsatisfiable, otherwise TRUE will be returned.
  bool UPLA::DoUpla(bool incremental)
  {
	if( _setting->verbosity > 2 )
	  { std::cout << "c Performing UPLA..." << std::endl; }

	struct rusage resources;
	getrusage(RUSAGE_SELF, &resources);
	double timeS = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
	
	_prepro->UpdateWatches(false);

	_prepro->_methodLimit = (1<<24);

	// Initialization.
	_nextcandidates.resize( _prepro->_variables+1,false );
	_activity.resize(_prepro->_variables+1, 0);
	_heap.resize(_prepro->_variables);

	_newEquivalences.clear();

	uint32_t constants = _statistics.uplaConstantVariables;
	uint32_t uplaequivs = _statistics.uplaEquivalentVariables;

	bool bcpres1(true);
	bool bcpres2(true);

	_okay = true;

	_heap.clear();
	
	uint32_t low(1);
	if( incremental )
	  { low = _prepro->_firstPreVarIndex; }

	for(uint32_t v = low; v <= _prepro->_variables; ++v )
	  {
		if( _prepro->_deleted[v] )
		  { continue; }

		_activity[v] = _prepro->OccurenceCount(v<<1) + _prepro->OccurenceCount((v<<1)^1);

		if( _activity[v] == 0 )
		  { 
			continue; 
		  }

		_heap.insert(v);
	  } 

	uint32_t tmpimplindex(_dsImplIndex);

	// Perform UPLA as long as mandatory implications can be found. 
	while (  !_heap.empty() )
	  {
		uint32_t tmppos = _dsEndIndex;

		while ( !_heap.empty() )
		  {
			// Initialization.
			uint32_t lit(_heap.top()<<1);

			if( _prepro->_methodLimit < 0 )
			  { break; }

			_nextcandidates[lit>>1] = false;
			_posimplications.clear();
			_negimplications.clear();

			_prepro->_methodLimit -= _activity[lit>>1];

			// Current variable still unassigned?
			if (_assignment[lit] || _assignment[lit ^ 1])
			  { continue; }

			// Add "lit^1" as a decision to the decision stack.
			AddUnit( lit^1 );
		  
			// Initialization.
			uint32_t pos(_dsEndIndex);

			// What about the effects of this decision?	
			bcpres1 = _prepro->_core->Deduce();

			// store all implications, triggered by decision of "lit^1"
			for(uint32_t i = pos; i < _dsEndIndex; ++i)
			  {
				_negimplications.push_back( _prepro->_decisionStack[i] ); 
			  }

			// Backtrack.
			Backtrack(lit^1);

			// If there are no conflict and no implications for "lit^1", we don't have to consider "lit"
			if( bcpres1 && _negimplications.empty() )
			  { continue;}

			// Add "lit" as a decision to the decision stack.
			AddUnit( lit );

			// What about the effects of this decision?		  
			bcpres2 = _core->Deduce();
				
			// store all implications, triggered by decision of "lit"
			for(uint32_t i = pos; i < _dsEndIndex; ++i)
			  { 
				_posimplications.push_back( _prepro->_decisionStack[i] ); 
			  }

			// Backtrack.
			Backtrack(lit);

			// Both assumption leads to a conflict => formula UNSAT
			if( !bcpres1 && !bcpres2 )
			  { 
				if( _setting->verbosity > 1 )
				  { std::cout << "c UPLA: both polarities of " << (lit>>1) << " are conflicting!" << std::endl; }
				_okay = false;
				break;
			  }
			
			// Assumption "lit" leads to a conflict => imply "lit^1"
			if( bcpres1 && !bcpres2 )
			  { 	
				// Add "real" implication and perform "real" deduction
				_core->AddImplication( lit^1 );

				++_statistics.uplaConstantVariables;

				if ( !_core->Deduce() )
				  { 
					_okay = false;
					break;
				  }

				_prepro->SetFlag(PF_UPLA);
			  }
			// Assumption "lit^1" leads to a conflict => imply "lit"
			else if( !bcpres1 && bcpres2 )
			  { 
				// Add "real" implication and perform "real" deduction
				_core->AddImplication( lit );

				++_statistics.uplaConstantVariables;

				if ( !_core->Deduce() )
				  { 
					_okay = false;
					break;
				  }

				_prepro->SetFlag(PF_UPLA);
			  }
			// No conflict, but maybe we have learnt some constants or equivalences
			else
			  {
				size_t possize (_posimplications.size());
				for( size_t i = 0; i != possize; ++i )
				  {
					uint32_t posimpl(_posimplications[i]);

					if( _assignment[posimpl] || _assignment[posimpl^1] || _prepro->_donttouch[posimpl>>1] ) 
					  { continue; }

					size_t negsize(_negimplications.size());
					_prepro->_methodLimit -= negsize;
					for( size_t j = 0; j != negsize; ++j )
					  {
						uint32_t negimpl(_negimplications[j]);

						// New Constant
						if( posimpl == negimpl )
						  {
							_prepro->SetFlag(PF_UPLA);

							++_statistics.uplaConstantVariables;

							// Constant already assigned to opposite polarity? => problem UNSAT
							if( _assignment[(posimpl)^1] )
							  { 
								_okay = false;
								break;
							  }

							_core->AddImplication( posimpl );
						  }
						// New Equivalence
						else if( posimpl == (negimpl^1) )
						  {
							// Add the binaries, which are responsible for the equivalence

							uint32_t firstlit(lit);
							uint32_t secondlit(posimpl);

							if( _prepro->_donttouch[secondlit>>1] )
							  { 
								if( !_prepro->_donttouch[firstlit>>1] )
								  { 
									secondlit = lit;
									firstlit = posimpl;
								  }
								// Both variables don't touch? Do nothing...
								else
								  { continue; }
							  }

							bool firstbin = _prepro->AddBinary(firstlit, secondlit^1, false);
							bool secondbin = _prepro->AddBinary(firstlit^1, secondlit, false);

							// if both binaries already exists, we have to do nothing...
							if( firstbin || secondbin )
							  {
								_newEquivalences.push_back(std::make_pair(lit, posimpl));
							  }
						  }
					  }
				  }

				if ( !_core->Deduce() )
				  { 
					_okay = false;
					break;
				  }
			  }
		  }

		if( _prepro->_methodLimit < 0 || !_okay )
		  { break; }

		// Every candidate is proceeded, calculate candidates for next round

		// Proceed occurence list of assigned literals for next upla round
		for( uint32_t j = tmppos; j < _dsEndIndex; ++j )
		  {
			uint32_t litOnStack(_prepro->_decisionStack[j]);

			// clear reasons for newly introduced implications on level 0 since corresponding reason clauses will be removed

			_prepro->_forcing[litOnStack>>1].ClearReason();
			const std::vector< CRef >& occur( _prepro->_occur[ litOnStack^1 ] );

			size_t osize( occur.size() );

			_prepro->_methodLimit -= osize;

			for( size_t c = 0; c != osize; ++c )
			  {
				Clause& clause = _ca[occur[c]];
				uint32_t csize( clause.size() );
					
				for( uint32_t opos = 0; opos != csize; ++opos )
				  {
					uint32_t clit = clause[opos];

					if( !_nextcandidates[ clit>>1 ] && !_assignment[clit] && !_assignment[clit^1] )
					  {
						_heap.insert(clit>>1);
						_nextcandidates[clit>>1] = true;
					  }
				  }
			  }
		  }
	  }

	// clear candidates flag for possibly unproceeded variables
	while ( !_heap.empty() )
	  {
		// Initialization.
		uint32_t var(_heap.top());
		_nextcandidates[var] = false;
	  }
 
	// Set Implication Pointer for next Simplification
	_dsImplIndex = tmpimplindex;

	// Order of variables changes during UPLA => reorder
	for( size_t c = 0 ; c != _prepro->_clauseDatabase.size(); ++c )
	  {	_ca[_prepro->_clauseDatabase[c]].Sort(); }

	// Remove n-nary and ternary from watch list
	for( uint32_t v = 1; v <= _prepro->_variables; ++v )
	  {
		if( _prepro->_binaries[v<<1].empty() && _prepro->_binaries[(v<<1)^1].empty() )
		  {	continue; }

		for( uint32_t literal = (v<<1); literal < ((v<<1)+2); ++literal )
		  {	
			std::vector<Watcher>& watches( _prepro->_binaries[literal] );

			size_t size( watches.size() );
			for( size_t i = 0; i != size; ++i )
			  {
				if( !watches[i].IsBinary() )
				  {
					watches[i--] = watches[--size];
					watches.pop_back();
				  }
			  }
			std::sort( watches.begin(), watches.end(), WatchedSorter() );
		  }
	  }

	// Propagate units
	_prepro->PropagateUnits();

	if ( _okay )
	  {
		// Now add the new equivalences
		for (size_t i = 0; i != _newEquivalences.size(); ++i)
		  {
			uint32_t firstlit = _newEquivalences[i].first; 
			uint32_t secondlit = _newEquivalences[i].second;

			if (_prepro->_deleted[firstlit>>1])
			  { 			
				// Both variables already deleted? Or deleted variable already replaced by "secondlit"?
				uint32_t replacingLit(_rebuilder->IsReplaced(firstlit));

				// Do nothing if either:
				// Both lits already eliminated or
				// Equivalance already performed
				if (_prepro->_deleted[secondlit>>1] || ( replacingLit == secondlit ))
				  { 
					continue; 
				  }

				// Change replace lit
				firstlit = replacingLit;
				assert( firstlit != 0 );
				assert( !_prepro->_deleted[firstlit>>1]); 
			  }
			else if (_prepro->_deleted[secondlit>>1])
			  {
				uint32_t replacingLit(_rebuilder->IsReplaced(secondlit));

				// "secondlit" already replaced by "firstlit"?
				if (replacingLit == firstlit)
				  { 
					continue; 
				  }
				// Change replace lit
				secondlit = replacingLit; 
				assert( secondlit != 0 );
			  }

			uint32_t toReplace(secondlit);
			uint32_t replace(firstlit);

			if (_prepro->_donttouch[toReplace>>1] )
			  { 
				if (!_prepro->_donttouch[replace>>1] )
				  { 
					toReplace = replace;
					replace = secondlit;
				  }
				// Both variables don't touch? Do nothing...
				else
				  { continue; }
			  }
			
			// We may have already deleted the literals
			if (!_prepro->_deleted[replace>>1])
			  {
				_prepro->RemoveBinary( replace^1, toReplace );
				_prepro->RemoveBinary( replace, toReplace^1 );
			  }

			if (!_prepro->_deleted[toReplace>>1])
			  { 
				_prepro->RemoveBinary( toReplace, replace^1 );
				_prepro->RemoveBinary( toReplace^1, replace );
			  }

			assert( !_prepro->_donttouch[toReplace>>1] );

			++_statistics.uplaEquivalentVariables;

			if (!_prepro->ReplaceVariable( toReplace, replace ))
			  { 
				_okay = false;
				break;
			  }
		  }
	  }

	getrusage(RUSAGE_SELF, &resources);
	double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
	
	_statistics.runtime_upla += (timeC-timeS);
	
	if( _setting->verbosity > 1 )
	  {
		std::cout << "c runtime upla           : " << timeC-timeS << "s (" << _statistics.runtime_upla << "s)" << std::endl;
		if ( (_statistics.uplaConstantVariables - constants) || (_statistics.uplaEquivalentVariables - uplaequivs) )
		  {
			std::cout << "c upla constants         : " << (_statistics.uplaConstantVariables - constants) << std::endl;
			std::cout << "c upla equivalences      : " << (_statistics.uplaEquivalentVariables - uplaequivs) << std::endl;
		  }
	  }

	return _prepro->PropagateUnits();      
  }
}
