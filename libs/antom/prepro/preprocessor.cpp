/********************************************************************************************
preprocessor.cpp -- Copyright (c) 2014-2016, Sven Reimer

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

#include "preprocessor.h"
#include "core.h"
#include "bce.h"
#include "bva.h"
#include "hte.h"
#include "upla.h"
#include "modelrebuilder.h"

#include <list>

/* For all routines:
   Add an extra loop, if something changes after unitPropagation, consider the subroutine again
   
   So far:
   At the end of every subroutine unitPropagation is performend, the subroutine does not take the result of the propagation into account
 */

namespace antom 
{
  // Constructor
  Preprocessor::Preprocessor(Core* core) :
	_core( core ),	 
	_control(core->_control),
	_setting(core->_setting),
	_bce(nullptr),
	_bva(nullptr),
	_hte(nullptr),
	_upla(nullptr),
	_rebuilder(nullptr),
	_methodLimit(0),
	_emptyClause(core->_emptyClause),
	_variables(core->_variables),
	_dsImplIndex(core->_dsImplIndex),
	_dsEndIndex(core->_dsEndIndex),
	_decisionLevel(core->_decisionLevel),
	_assignment(core->_assignment),
	_level(core->_level),
	_decisionStack(core->_decisionStack),
	_model(core->_model),
	_forcing(core->_forcing),
	_deleted(core->_deleted),
	_clauseDatabase(core->_clauseDatabase),
	_binaries(core->_watches),
	_ca(core->_ca),
	// initialize data structures
	_occur(),
	_occurCounts(),
	_occurenceCount(),
	_donttouch(),
	_varCandidates(),
	_binCandidates(helper::DescendingOrder<size_t>(_occurCounts)),
	_firstPreVarIndex(1),
	_firstPreClauseIndex(1),
	_lastImplIndex(1),
	_statistics(core->_statistics),
	_pFlag(PF_NOTHING),
	_okay(true)
  {
	assert(_control != nullptr);
	assert(_setting != nullptr);

	_rebuilder = new ModelRebuilder(this);
	_bce = new BCE(this);
	_bva = new BVA(this);
	_hte = new HTE(this);
	_upla = new UPLA(this);
  }

  Preprocessor::~Preprocessor(void)
  {
	delete _rebuilder;
	delete _bce;
	delete _bva;
	delete _hte;
	delete _upla;
  }

  inline bool Preprocessor::OccurenceSorter::operator()(uint32_t l1, uint32_t l2) const
  {
	return ( occurSizes[l1] < occurSizes[l2]) ;
  }

  // Simplifies the current CNF formula by performing some preprocessing steps.
  // Returns ANTOM_UNSAT if the formula is unsatisfiable, otherwise ANTOM_UNKNOWN
  uint32_t Preprocessor::Preprocess(PreproType type) 
  {
	// If the solver is not on decision level 0, we might have a problem.
	assert(_decisionLevel == 0);
	assert(!_control->GetTimeOutReached());

	if (_emptyClause )
	  { return ANTOM_UNSAT; }

	switch( type )
	  {
	  case INPROCESS :
		++_statistics.inprocessings; 
		if (_setting->verbosity > 1 )
		  { std::cout << "c Inprocessing No. " << _statistics.inprocessings << std::endl; }
		break;
	  case PREPROCESS :
		if (_setting->verbosity > 1 )
		  { std::cout << "c Start preprocessing..." << std::endl; }
		break;
	  case INCREMENTAL :
		if (_setting->verbosity > 1 )
		  { std::cout << "c Start incremental preprocessing..." << std::endl; }
		break;
	  default : 
		assert(false);
		break;
	  }

	double timeC(0.0);
	struct rusage resources;
	getrusage(RUSAGE_SELF, &resources);
	double timeS = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
	_core->_control->SetStartTime(timeS);

	uint32_t loops(0);
	uint32_t result(ANTOM_UNKNOWN);

	// Some statistics
	uint32_t tcons = _statistics.constantVariables;
	uint32_t tequiv = _statistics.equivalentVariables;
	uint32_t tuplacons = _statistics.uplaConstantVariables;
	uint32_t tuplaequiv = _statistics.uplaEquivalentVariables;
	uint32_t tconsbysat = _statistics.constantVariablesBySAT;
	uint32_t tresvars = _statistics.resolvedVariables;
	uint32_t treslits = _statistics.resolvedLiterals;
	uint32_t tblockedclauses = _statistics.blockedClauses;
	uint32_t thiddentautologies = _statistics.hiddenTautologies;
	uint32_t thiddensubsumptions = _statistics.hiddenSubsumptions;
	uint32_t tmonolits = _statistics.monotoneVariables;
	uint32_t tdclits = _statistics.dontcareVariables;
	uint32_t tsubsumptions = _statistics.subsumptions;
	uint32_t tselfsubsumptions = _statistics.selfSubsumptions;
	uint32_t tbvavars = _statistics.bvaVariables;
	uint32_t tbvalits = _statistics.bvaLiterals;
	uint32_t tvivifySubsumptions = _statistics.vivifySubsumptions;
	uint32_t tvivifyUnits = _statistics.vivifyUnits;
	uint32_t tvivifyDiff = _statistics.vivifyDiff;
	uint32_t tunitprop = _statistics.unitPropagations;

	bool doitagain( true );
	_pFlag = PF_ALL;
	bool done(false);

	bool fast(_dsImplIndex != _dsEndIndex);

	// Extract binaries, build occurence lists
	PreparePreprocessing(type);

	// Always start over in incremental mode
	if (_setting->incrementalMode )
	  { _lastImplIndex = 1; }

	_dsImplIndex = _lastImplIndex;

	if (!PropagateUnits() ) 
	  { 
		_control->SetDone(); 
		result = ANTOM_UNSAT; 
		done = true; 
	  }
	if (!done && !FastPreprocess(fast, type==INCREMENTAL) ) 
	  { 
		_control->SetDone(); 
		result = ANTOM_UNSAT; 
		done = true; 
	  }
	doitagain = fast;

	if (!done )
	  {

		// main preprocessing loop
		do {
		  ++loops;

		  if (_setting->verbosity > 1 )
			{
			  std::cout << "c pre/inpro loop #" << loops << std::endl;
			}

		  doitagain = false;

		  // clear vivify part;
		  UnsetFlag(PF_VIVIFY);

		  // Perform Vivification?
		  if (_setting->doVivification )
			{
			  // no changes since last vivification? -> end loop
			  if (_pFlag == PF_NOTHING )
				{ break; }
			  
			  if (!Vivify() )
				{ _control->SetDone(); result = ANTOM_UNSAT; break; }

			  bool fast = GetFlag(PF_VIVIFY);
			  if (!FastPreprocess(fast, type==INCREMENTAL) ) 
				{ _control->SetDone(); result = ANTOM_UNSAT; break; }
			  doitagain |= fast;

			  if (_control->ReachedLimits() )
				{ 
				  getrusage(RUSAGE_SELF, &resources); 
				  timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
				  
				  _statistics.runtime_preprocessing += (timeC-timeS);
				  return ANTOM_UNKNOWN;
				}
			}


		  // clear upla part; 0x7b = 11111011
		  UnsetFlag(PF_UPLA);

		  // Perform UPLA?
		  if (_setting->doUpla )
			{
			  // no changes since last upla? -> end loop
			  if (_pFlag == PF_NOTHING )
				{ break; }

			  if (!_upla->DoUpla(type==INCREMENTAL) ) 
				{ _control->SetDone(); result = ANTOM_UNSAT; break; }

			  bool fast = GetFlag(PF_UPLA);
			  if (!FastPreprocess(fast, type==INCREMENTAL) ) 
				{ _control->SetDone(); result = ANTOM_UNSAT; break; }
			  doitagain |= fast;

			  if (_control->ReachedLimits() )
				{ 
				  getrusage(RUSAGE_SELF, &resources); 
				  timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
				  
				  _statistics.runtime_preprocessing += (timeC-timeS);
				  return ANTOM_UNKNOWN;
				}
			}

		  // clear hte part; 0x7f = 01111111
		  UnsetFlag(PF_HTE);

		  // Perform HTE?
		  if (_setting->doHte )
			{
			  // no changes since last hte? -> end loop
			  if (_pFlag == PF_NOTHING )
				{ break; }

			  _hte->DoHiddenTautologyClauseElimination();

			  if (_control->ReachedLimits() )
				{ 
				  getrusage(RUSAGE_SELF, &resources); 
				  timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
				  
				  _statistics.runtime_preprocessing += (timeC-timeS);
				  return ANTOM_UNKNOWN;
				}
			}

		  // clear bce part; 0xbf = 10111111
		  UnsetFlag(PF_BCE);

		  if (_setting->doBce )
			{
			  // no changes since last bce? -> end loop
			  if (_pFlag == PF_NOTHING )
				{ break; }

			  if (!_bce->DoBlockedClauseElimination() )
				{ _control->SetDone(); result = ANTOM_UNSAT; break; }

			  bool fast = GetFlag(PF_BCE);
			  if (!FastPreprocess(fast, type==INCREMENTAL) ) 
				{ _control->SetDone(); result = ANTOM_UNSAT; break; }
			  doitagain |= fast;
			}
		  
		  // clear variable elimination part; 0xf7 = 11110111
		  UnsetFlag(PF_VAR_ELIM);

		  // Perform variable elimination?
		  if (_setting->doVarElimination )
			{
			  // no changes since last variable elimination? -> end loop
			  if (_pFlag == PF_NOTHING )
				{ break; }

			  if (!VarElimination(type==INCREMENTAL) ) 
				{ _control->SetDone(); result = ANTOM_UNSAT; break; }

			  bool fast = GetFlag(PF_VAR_ELIM);
			  if (!FastPreprocess(fast, type==INCREMENTAL) ) 
				{ _control->SetDone(); result = ANTOM_UNSAT; break; }
			  doitagain |= fast;

			  if (_control->ReachedLimits() )
				{ 
				  getrusage(RUSAGE_SELF, &resources); 
				  timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
				  
				  _statistics.runtime_preprocessing += (timeC-timeS);
				  return ANTOM_UNKNOWN;
				}
			}

		  // clear subsumption part; 0xfd = 11111101
		  UnsetFlag(PF_SUBSUMPTION);

		  // Perform full subsumption check?
		  if (_setting->doSubsumption )
			{
			  // no changes since last subsumption? -> end loop
			  if (_pFlag == PF_NOTHING )
				{ break; }

			  if (!FullSubsumptionCheck(type==INCREMENTAL) ) 
				{ _control->SetDone(); result = ANTOM_UNSAT; break;}

			  bool fast = GetFlag(PF_SUBSUMPTION);
			  if (!FastPreprocess(fast, type==INCREMENTAL) ) 
				{ _control->SetDone(); result = ANTOM_UNSAT; break;}
			  doitagain |= fast;

			  if (_control->ReachedLimits() )
				{ 
				  getrusage(RUSAGE_SELF, &resources); 
				  timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
				  
				  _statistics.runtime_preprocessing += (timeC-timeS);
				  return ANTOM_UNKNOWN;
				}
			}

		  // clear bva part; 0xdf = 11011111
		  UnsetFlag(PF_BVA);

		  if (_setting->doBva )
			{
			  // no changes since last variable elimination? -> end loop
			  if (_pFlag == PF_NOTHING )
				{ break; }

			  if (!_bva->BoundedVariableAddition() )
				{ _control->SetDone(); result = ANTOM_UNSAT; break; }

			  bool fast = GetFlag(PF_BVA);
			  if (!FastPreprocess(fast, type==INCREMENTAL) ) 
				{ _control->SetDone(); result = ANTOM_UNSAT; break; }
			  doitagain |= fast;

			  if (_control->ReachedLimits() )
				{ 
				  getrusage(RUSAGE_SELF, &resources); 
				  timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
				  
				  _statistics.runtime_preprocessing += (timeC-timeS);
				  return ANTOM_UNKNOWN;
				}

			}
		}
		while( doitagain && (loops < _setting->maxLoops) && !_core->_control->Done() );
	  }

	if (type==INCREMENTAL )
	  {
		_firstPreVarIndex = _variables+1;
		_firstPreClauseIndex = static_cast<uint32_t>(_clauseDatabase.size())+1;
	  }

	getrusage(RUSAGE_SELF, &resources); 
	timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	_statistics.runtime_preprocessing += (timeC-timeS);
	
	if (_setting->verbosity > 1 )
	  {
		std::cout << "c -----------------------------" << std::endl;
		switch( type )
		  {
		  case PREPROCESS:
			std::cout << "c After Preprocessing " << std::endl
					  << "c prepro loops           : " << loops << std::endl;
			break;
		  case INPROCESS:
			std::cout << "c After Inprocessing " << std::endl
					  << "c inpro loops            : " << loops << std::endl;
			break;
		  case INCREMENTAL:
			std::cout << "c After Incremental Preprocessing " << std::endl
					  << "c prepro loops           : " << loops << std::endl;
			break;
		  default:
			assert(false);
			break;
		  }

		std::cout << "c binary constants       : " << (_statistics.constantVariables - tcons) << " (" << _statistics.constantVariables << ")" <<std::endl
				  << "c binary equivalences    : " << (_statistics.equivalentVariables - tequiv) << " (" << _statistics.equivalentVariables << ")" << std::endl
				  << "c upla constants         : " << (_statistics.uplaConstantVariables - tuplacons) << " (" << _statistics.uplaConstantVariables << ")" <<  std::endl
				  << "c upla equivalences      : " << (_statistics.uplaEquivalentVariables - tuplaequiv) << " (" << _statistics.uplaEquivalentVariables << ")" <<  std::endl
				  << "c constants by SAT check : " << (_statistics.constantVariablesBySAT - tconsbysat) << " (" << _statistics.constantVariablesBySAT << ")" <<  std::endl
				  << "c variable eliminations  : " << (_statistics.resolvedVariables - tresvars) << " (" << _statistics.resolvedVariables << ")" <<  std::endl
				  << "c elim literal reduction : " << (_statistics.resolvedLiterals - treslits) << " (" << _statistics.resolvedLiterals << ")" <<  std::endl
				  << "c blocked clauses        : " << (_statistics.blockedClauses - tblockedclauses) << " (" << _statistics.blockedClauses << ")" <<  std::endl
				  << "c hidden tautologies     : " << (_statistics.hiddenTautologies - thiddentautologies) << " (" << _statistics.hiddenTautologies << ")" <<  std::endl
				  << "c hidden subsumptions    : " << (_statistics.hiddenSubsumptions - thiddensubsumptions) << " (" << _statistics.hiddenSubsumptions << ")" <<  std::endl
				  << "c monotone variables     : " << (_statistics.monotoneVariables - tmonolits) << " (" << _statistics.monotoneVariables << ")" <<  std::endl
				  << "c dc variables           : " << (_statistics.dontcareVariables - tdclits) << " (" << _statistics.dontcareVariables << ")" <<  std::endl
				  << "c subsumed clauses       : " << (_statistics.subsumptions - tsubsumptions) << " (" << _statistics.subsumptions << ")" <<  std::endl
				  << "c selfsubsumed literals  : " << (_statistics.selfSubsumptions - tselfsubsumptions) << " (" << _statistics.selfSubsumptions << ")" <<  std::endl
				  << "c bva variables          : " << (_statistics.bvaVariables - tbvavars) << " (" << _statistics.bvaVariables << ")" <<  std::endl
				  << "c bva literal reduction  : " << (_statistics.bvaLiterals - tbvalits) << " (" << _statistics.bvaLiterals << ")" <<  std::endl
				  << "c vivify subsumptions    : " << (_statistics.vivifySubsumptions - tvivifySubsumptions) << " (" << _statistics.vivifySubsumptions << ")" << std::endl
				  << "c vivify units           : " << (_statistics.vivifyUnits - tvivifyUnits) << " (" << _statistics.vivifyUnits << ")" << std::endl
				  << "c vivify lit reduction   : " << (_statistics.vivifyDiff - tvivifyDiff) << " (" << _statistics.vivifyDiff << ")" << std::endl
				  << "c unit propagations      : " << (_statistics.unitPropagations - tunitprop) << " (" << _statistics.unitPropagations << ")" << std::endl
				  << "c pre/inpro time         : " << (timeC-timeS) << "s (" << _statistics.runtime_preprocessing << "s)" << std::endl;
		
	  }

	for (uint32_t v = 1; v <= _variables; ++v )
	  {
		std::vector<CRef>().swap(_occur[v<<1]);
		std::vector<CRef>().swap(_occur[(v<<1)^1]);
	  }

	// Rebuild watches
	UpdateWatches();

	// garbage collect after watchupdates, so that watches are allocated together
	_core->CheckGarbage();

	// Update implication index for next inprocessing
	_lastImplIndex = _dsEndIndex;

	// Everything went fine.
	return result;
  }

  // Some fast preprocessing routines, which can be called often
  bool Preprocessor::FastPreprocess(bool& didsomething, bool incremental)
  {
	assert( _dsImplIndex == _dsEndIndex );

	DetectMonotone(incremental);
	if (!FindBinaryConsAndEquiv(didsomething) ) { return false; }

	return true;
  }

  void Preprocessor::ExtendModel(void)
  {
	_rebuilder->ExtendModel();
  }

  void Preprocessor::ClearRestoreData(uint32_t begin, uint32_t end)
  {
	_rebuilder->ClearRestoreData(begin, end);
  }
	
  // Resets the Preprocessor. The following status flags/variables remain untouched:
  // * The references of the core
  // * Information about replaced variables and clauses 
  // * Don't touch variables
  // * Overall statistics 
  // * DoPreprocessing and DoInprocessing flags
  void Preprocessor::Reset(void)
  {
	InstanceReset();
	_setting->ResetPrepro();
  }

  void Preprocessor::InstanceReset(void)
  {
	// Reset all variable related data structures.
	for (uint32_t i = 2; i < _occur.size(); ++i)
	  {
		std::vector< CRef >().swap( _occur[i] );
	  }

	_lastImplIndex = 1;
	_statistics.ResetPrepro();
	_pFlag = PF_NOTHING;
	_okay = true;
  }

  /* Begin: some helper functions */

  // Updates all data structures depending on the number of variables to be able to handle "_variables" variables.

  void Preprocessor::UpdateDataStructures(void)
  {
	_rebuilder->_replacedBy.resize(_variables+1,0);
	_donttouch.resize(_variables+1, false);

	_varCandidates.resize(_variables+1);
	
	// Initialization.
	uint32_t max((_variables << 1) + 2);	
	_occurCounts.resize(max,0);	
	  
	_occur.resize(max);
  }

  // Adds (lit1 + lit2) to binaries. 
  // If (lit1 + lit2) already exists, do nothing#
  // Returns false, if binary already exists
  bool Preprocessor::AddBinary(uint32_t lit1, uint32_t lit2, bool learned)
  {
	// Check if new binary is already in watch list
	uint32_t firstlit(lit1);
	uint32_t seclit(lit2);

	if (_binaries[seclit].size() < _binaries[firstlit].size() )
	  {
		firstlit = lit2;
		seclit = lit1;
	  }

	size_t f( 0 );
	size_t fsize( (uint32_t)_binaries[firstlit].size() );
	 
	for (; f != fsize; ++f )
	  {
		// We potentially add binaries in upla, there are possible non-binaries in watchlist
		if (_binaries[firstlit][f].IsBinary() && _binaries[firstlit][f].GetSecondLit() == seclit )
		  { 
			// If we add not non-learned binary, and we found a learned -> change learned to unlearned
			// in all other cases the status of the found binary is untouched
			if (!learned && _binaries[firstlit][f].IsLearnedBinary() )
			  { 
				_binaries[firstlit][f].SetLearnedBinary(false);
				// search for corresponding binary
				uint32_t j(0);
				size_t jsize( (uint32_t)_binaries[seclit].size() );
				for (; j != jsize; ++j )
				  {
					if (_binaries[seclit][j].IsBinary() && _binaries[seclit][j].GetSecondLit() == firstlit )
					  {
						_binaries[seclit][j].SetLearnedBinary(false);
						break;
					  }
				  }
				assert( j < jsize );
			  }
			break; 
		  }

	  }
	// Add new binary
	if (f == fsize )
	  {
		_binaries[lit1].push_back(Watcher(lit2,learned));
		_binaries[lit2].push_back(Watcher(lit1,learned));
		++_statistics.currentBinaryClauses += 2;

		return true;
	  }

	return false;
  }

  // Remove the "pos"th binary entry in the binary list of "lit"
  // "pos" and "size" will be updated
  void Preprocessor::RemoveBinary(uint32_t lit, uint32_t& pos, uint32_t& size)
  {
	_binaries[lit][pos--] = _binaries[lit][--size];
	_binaries[lit].pop_back();
	--_statistics.currentBinaryClauses;
  }

  // Remove the binary (lit1 + lit2) from binary list of "lit1"
  void Preprocessor::RemoveBinary(uint32_t lit1, uint32_t lit2)
  {
	size_t i( 0 );
	size_t size( _binaries[lit1].size() );

	for (; i != size; ++i )
	  {
		if (_binaries[lit1][i].GetSecondLit() == lit2 )
		  {
			break;
		  }
	  }

	assert( i != size );

	// copy remaining watches
	size_t j = i;
	++i;
	for (; i != size; ++i, ++j )
	  {
		_binaries[lit1][j] = _binaries[lit1][i];
	  }
	_binaries[lit1].resize(j);
  }

  void Preprocessor::RemoveWatch( uint32_t lit, CRef c)
  {
	size_t j = 0;
	size_t size = _binaries[lit].size();
	for (; j != size; ++j)
	  {
		if( !_binaries[lit][j].IsClause() )
		  {
			continue;
		  }
		
		if (_binaries[lit][j].GetClause() == c)
		  {
			_binaries[lit][j--] = _binaries[lit][--size];
			_binaries[lit].pop_back();
			break;
		  }
	  }
	assert( j != size );
  }

  bool Preprocessor::HasBinary(uint32_t lit1, uint32_t lit2) const
  {
	for (uint32_t i = 0; i != _binaries[lit1].size(); ++i)
	  {
		if (lit2 == _binaries[lit1][i].GetSecondLit())
		  {
			return true;
		  }
	  }
	return false;
  }

  // Removes "clause" from occurence list of "lit"
  void Preprocessor::RemoveFromOccurenceList(uint32_t lit, CRef cr)
  {
	size_t osize( _occur[ lit ].size() );
	size_t f( 0 );
	for (; f != osize; ++f )
	  {
		if (cr == _occur[lit][f] )
		  {
			_occur[lit][f--] = _occur[lit][--osize];
			_occur[lit].pop_back();
			break;
		  }
	  }

	assert( f != osize );
  }

  // Removes every occurence of "clause"
  // Mark "clause" as "to deleted"
  void Preprocessor::ClearAllOccurences(CRef cr)
  {
	Clause& clause( _ca[cr] );
	//std::cout << __func__ << std::endl;
	//clause.print();
	// In case length was equal to "2", we have deleted the third literal in advance
	if (clause.size() <= 3 )
	  { --_statistics.currentTernaryClauses; }
	else
	  { --_statistics.currentNaryClauses; }

	uint32_t size(clause.size());
	for (uint32_t pos = 0; pos != size; ++pos)
	  {
		RemoveFromOccurenceList( clause[pos], cr );
	  }
	
	// Mark clause as to be deleted
	clause.MarkForDeletion();
  }

  void Preprocessor::ClearClauseDatabase(void)
  {
	size_t csize(_clauseDatabase.size());
	for (size_t c = 0; c != csize; ++c )
	  {
		// delete marked clauses permanatly
		if (_ca[_clauseDatabase[c]].ToDelete() )
		  {
			_ca.Free(_clauseDatabase[c]);
			_clauseDatabase[c--] = _clauseDatabase[--csize];
			_clauseDatabase.pop_back();
		  }
	  }
  }

  void Preprocessor::CountOccurences(bool countLearned)
  {
	_occurCounts.clear();
	_occurCounts.resize( (_variables+1)<<1, 0 );

	for (uint32_t v = 1; v <= _variables; ++v )
	  {

		uint32_t poslit(v<<1);
		uint32_t neglit((v<<1)^1);
		_occurCounts[poslit] = OccurenceCount(poslit, countLearned);
		_occurCounts[neglit] = OccurenceCount(neglit, countLearned);
	  }
  }

  void Preprocessor::CountBinaryOccurences(void)
  {
	_occurCounts.clear();
	_occurCounts.resize( _variables+1, 0 );
	for (uint32_t v = 1; v <= _variables; ++v )
	  {
		uint32_t poslit(v<<1);
		uint32_t neglit((v<<1)^1);
		_occurCounts[v] = _binaries[poslit].size() + _binaries[neglit].size();
	  }
  }

  // Returns number occurences of "lit" (including binaries)
  uint32_t Preprocessor::OccurenceCount(uint32_t lit, bool countLearned) const
  {
	if (countLearned )
	  {
		return static_cast<uint32_t>(_binaries[lit].size() + _occur[lit].size());
	  }

	uint32_t count(0);
	for (uint32_t i = 0; i != _binaries[lit].size(); ++i )
	  {
		if (!_binaries[lit][i].IsLearnedBinary() )
		  {
			++count;
		  }
	  }
	for (uint32_t i = 0; i != _occur[lit].size(); ++i )
	  {
		if (!_ca[_occur[lit][i]].IsLearned() )
		  {
			++count;
		  }
	  }
	return count;
  }

  uint32_t Preprocessor::GetClauseSize(const Watcher& watcher) const
  {
	switch (watcher.GetType())
	  {
	  case BINARY : return 2; break;
	  case TERNARY : return 3; break;
	  case NNARY : return _ca[watcher.GetClause()].size(); break;
	  default : assert(false); return 0; break;
	  }
  }

  // Removes "literal" from clause
  // Eventually introduces new binary and mark old n-nary for deletion
  // Returns "true" if n-nary clause is deleted, "false" otherwisse
  bool Preprocessor::StrengthenClause(CRef cr, uint32_t literal, bool keepoccurence)
  {
	Clause& clause( _ca[cr] );
	bool result(false);
	uint32_t pos( 0 );
	uint32_t qos( 0 );
	uint32_t sign( 0 );
	uint32_t size( clause.size() );

	// Remove "literal" from clause
	for (; pos != size; ++pos )
	  {
		if (clause[ pos ] != literal )
		  {
			clause[ qos ] = clause[ pos ];
			//sign |= 1ULL << (clause[qos] % 64);
			sign |= 1 << (clause[qos] % 32);
			++qos;
		  }
	  }

	assert( pos != qos );

	// Is "clause" before deletion of "literal" a ternary clause?
	// Then we got a binary
	if (clause.size() == 3 )
	  {
		assert(qos == 2 );
		// Add the new binary
		AddBinary(clause[0], clause[1], clause.IsLearned());

		clause.SetLength( 2 );
		ClearAllOccurences(cr);
		result = true;
	  }
	else
	  {
		_ca.FreeLiterals(pos-qos);
		clause.SetSign(sign);
		assert( qos == (size-1) );
		clause.SetLength( qos );
		if (qos == 3 )
		  { ++_statistics.currentTernaryClauses; --_statistics.currentNaryClauses; }
	  }

	if (!keepoccurence )
	  { RemoveFromOccurenceList( literal, cr ); }

	return result;
  }

  // 1. Extract all binaries of watchlist and push them in extra data structure
  // (this can also be done during "addClause", if preprocessing is enabled)
  // 2. Creates occurence lists
  // Assumes that there are no duplicated binaries in "_watches"
  void Preprocessor::PreparePreprocessing(PreproType type)
  {
	// Allocate memory for preprocessing
	UpdateDataStructures();

	// Reserve memory for varheaps
	for (uint32_t v = 0; v != _core->_noOfVarGroups; ++v )
	  { 
		_core->_varOrder[v]->resize(_variables); 
	  }

	_statistics.usedVariables = 0;
	_statistics.currentBinaryClauses = 0;
	_statistics.currentTernaryClauses = 0;
	_statistics.currentNaryClauses = 0;

	// Clear binary and n-nary occurence lists
	// Remove n-nary and ternary from watch list
	for (uint32_t v = 1; v <= _variables; ++v )
	  {
		if( !_occur[v<<1].empty() )
		  {
			std::cout << "occur of " << v << " not empty" << std::endl;
		  }
		assert( _occur[v<<1].empty() );
		assert( _occur[(v<<1)^1].empty() );

		if (_binaries[v<<1].empty() && _binaries[(v<<1)^1].empty() )
		  {	continue; }

		for (uint32_t literal = (v<<1); literal < ((v<<1)+2); ++literal )
		  {	
			std::vector<Watcher>& watches( _binaries[literal] );

			size_t size( watches.size() );

			for (size_t i = 0; i != size; ++i )
			  {
				if (!watches[i].IsBinary() )
				  {
					// Remove bon-binary
					watches[i--] = watches[--size];
					watches.pop_back();
				  }
				else
				  {
					++_statistics.currentBinaryClauses;
				  }
			  }
	
			std::sort( watches.begin(), watches.end(), WatchedSorter() );
		  }
	  }

	// Build occurence lists
	size_t size(_clauseDatabase.size());
	for (size_t c = 0; c != size; ++c )
	  {
		// Get a pointer to the next clause.
		CRef cr( _clauseDatabase[c] );
		Clause& clause( _ca[cr] );

		if (clause.size() == 3 )
		  { ++_statistics.currentTernaryClauses; }
		else
		  { ++_statistics.currentNaryClauses; }

		assert( !clause.ToDelete() );

		// sort clause... needed for subsumption checks
		clause.Sort();
		// calc signatur for subsumption checks
		clause.CalcSign();

		uint32_t size( clause.size());
		for (uint32_t pos = 0; pos != size; ++pos ) 
		  {
			_occur[clause[pos]].push_back(cr);
		  }
	  }

	// Count used variables
	for (uint32_t v = 1; v <= _variables; ++v )
	  {
		uint32_t poslit(v<<1);
		uint32_t neglit((v<<1)^1);
		if (_binaries[poslit].empty() && _binaries[neglit].empty() && _occur[poslit].empty() && _occur[neglit].empty() )
		  {	continue; }
		++_statistics.usedVariables;
		_deleted[v] = false;
	  }	

	if (_setting->verbosity > 1 )
	  {
		std::cout << "c ----------------------------" << std::endl;

		switch( type )
		  {
		  case PREPROCESS:
			std::cout << "c Before Preprocessing " << std::endl;
			break;
		  case INPROCESS:
			std::cout << "c Before Inprocessing " << std::endl;
			break;
		  case INCREMENTAL:
			std::cout << "c Before Incremental Preprocessing " << std::endl;
			break;
		  default:
			assert(false);
			break;
		  }

		_statistics.PrintClauseStats();
	  }
  }

  // Update all watchlists
  void Preprocessor::UpdateWatches(bool showstats)
  {
	// Now update watches for n-nary clauses
	size_t size = _clauseDatabase.size();
	for (size_t c = 0; c != size; ++c )
	  {
		CRef cr = _clauseDatabase[c];
		assert( !_ca[cr].ToDelete() );

		_core->AttachClause(cr, false);
	  }

	if (showstats && _setting->verbosity > 1 )
	  {
		_statistics.PrintClauseStats();
	  }
  }

  bool Preprocessor::CopyBinaryList(uint32_t toReplace, uint32_t replace)
  {
	// Proceed positive occurences of binaries first
	size_t bsize( _binaries[toReplace].size() );

	for (size_t i = 0; i != bsize; ++i )
	  {
		uint32_t seclit( _binaries[toReplace][i].GetSecondLit() );
		
		// Find constant...
		if (replace == seclit )
		  {
			++_statistics.constantVariables;

			// Variable already assigned in opposite polarity?
			if (_assignment[seclit^1] )
			  { return false; }
			
			if (!_assignment[seclit] )
			  { _core->AddImplication( seclit ); }

			// Remove the dual binary clause of seclit
			RemoveBinary( seclit, toReplace );
		  }
		// Skip the responsible clause for the equivalence
		else if ((replace^1) != seclit )
		  { 
			// Add new binary
			AddBinary( replace, seclit, _binaries[toReplace][i].IsLearnedBinary() );

			// Remove the dual binary clause of seclit
			RemoveBinary( seclit, toReplace );
		  }
	  }

	// Clear binary list of "toReplace"
	_binaries[toReplace].clear();
	return true;
  }

  void Preprocessor::CopyOccurenceList(uint32_t toReplace, uint32_t replace)
  {
	// Consider the positive occurences of the n-nary clauses
	size_t osize( _occur[toReplace].size() );

	for (size_t i = 0; i != osize; ++i )
	  {
		CRef cr = _occur[toReplace][i];
		Clause& clause = _ca[cr];

		// clause already marked as deleted?
		if (clause.ToDelete() )
		  { continue; }

		int32_t occ = clause.HasLitPolarity(replace);

		// After replacement clause would consists of ( l1 + ... + replace + ~replace ) 
		// Clause is satisfied after replacement -> delete occurences of clause and mark for deletion
		if (occ == -1 )
		  {

			ClearAllOccurences(cr);
			// Clause is also deleted in current vector -> update values
			--i; --osize;
		  }
		// Replaced literal is already in clause -> remove "toReplace" from clause
		else if (occ == +1 )
		  {
			// Keep occurence of clause, since it will be deleted anyway afterwards
			StrengthenClause( cr, toReplace, true );
		  }
		// just replace
		else
		  {
			uint32_t size( clause.size() );
			
			for (uint32_t pos = 0; pos != size; ++pos )
			  {
				if (clause[pos] == toReplace )
				  {
					clause[pos] = replace;
					break;
				  }
			  }

			clause.Sort();
			// Recalc the signature
			clause.CalcSign();

			// Update occurence list of "replace"
			_occur[replace].push_back(cr);			
		  }
	  }

	// Delete the old occurence list
	_occur[toReplace].clear();
  }

  // Replace every occurence of "toReplace" with "replace"
  // Refresh occurence lists
  // Preserves model for replaced variable
  bool Preprocessor::ReplaceVariable(uint32_t toReplace, uint32_t replace)
  {
	//std::cout << __FUNCTION__ << " toReplace: " << helper::Lit(toReplace) << " replace: " << replace << std::endl;

	assert( !_donttouch[toReplace>>1] );
	// If "toReplace" and "replace" are the same varible, we'll run into trouble
	assert( (toReplace>>1) != (replace>>1) );

	_statistics.currentBinaryClauses -= static_cast<uint32_t>( _binaries[toReplace].size() + _binaries[toReplace^1].size() );

	// Copy Binaries and occurences
	if (!CopyBinaryList(toReplace,replace) )
	  { return false; }
	if (!CopyBinaryList(toReplace^1,replace^1) )
	  { return false; }
	CopyOccurenceList(toReplace,replace);
	CopyOccurenceList(toReplace^1,replace^1);

	_rebuilder->AddVarEquivalence(toReplace, replace);

	// Mark "toReplace" as deleted
	_deleted[toReplace>>1] = true;

	return PropagateUnits();
  }

  // Delete all occurences and clauses of "var" in database
  void Preprocessor::DeleteVariable(uint32_t var)
  {
	assert( !_donttouch[var] );

	// Delete binaries
	uint32_t poslit(var<<1);
	uint32_t neglit((var<<1)^1);

	_statistics.resolvedLiterals += static_cast<uint32_t>(_binaries[poslit].size()<<1);
	_statistics.resolvedLiterals += static_cast<uint32_t>(_binaries[neglit].size()<<1);
	for (size_t i = 0; i != _binaries[poslit].size(); ++i )
	  {
		RemoveBinary( _binaries[poslit][i].GetSecondLit(), poslit );
	  }
	for (size_t i = 0; i != _binaries[neglit].size(); ++i )
	  {
		RemoveBinary( _binaries[neglit][i].GetSecondLit(), neglit );
	  }

	_statistics.currentBinaryClauses -= static_cast<uint32_t>( _binaries[poslit].size() + _binaries[neglit].size() );
	_binaries[poslit].clear();
	_binaries[neglit].clear();
	
	while( !_occur[poslit].empty() )
	  { 
		_statistics.resolvedLiterals += _ca[*_occur[poslit].begin()].size();
		ClearAllOccurences( *_occur[poslit].begin() ); 
	  }

	while( !_occur[neglit].empty() )
	  { 
		_statistics.resolvedLiterals += _ca[*_occur[neglit].begin()].size();
		ClearAllOccurences( *_occur[neglit].begin() ); 
	  }

	// Mark "var" as deleted
	_deleted[var] = true;
	assert( _forcing[var].NoReason());
	--_statistics.usedVariables;

	// Push fake assignment on stack
	if (!_assignment[var<<1] && !_assignment[(var<<1)^1])
	  { 
		_core->AddImplication(var<<1); 
		++_dsImplIndex;
	  }
  }

  // Merges two n-nary clauses with common literal "reslit"
  // Store new clause in "newClause"
  void Preprocessor::MergeClauses(Clause& clause1, Clause& clause2, uint32_t reslit, std::vector< uint32_t >& newClause) const
  {
	newClause.clear();
	newClause.reserve( clause1.size() + clause2.size() - 2 );
  
	uint32_t negreslit( reslit^1 );

	uint32_t clause1size( clause1.size() );
	uint32_t clause2size( clause2.size() );

	for (uint32_t i = 0; i != clause1size; ++i )
	  {
		assert( !_assignment[clause1[i]] && !_assignment[clause1[i]^1] );
		if (clause1[i] != reslit )
		  { newClause.push_back(clause1[i]); }
	  }
	for (uint32_t i = 0; i != clause2size; ++i )
	  {
		assert( !_assignment[clause2[i]] && !_assignment[clause2[i]^1] );
		if (clause2[i] != negreslit )
		  { newClause.push_back(clause2[i]); }
	  }
	
	std::sort( newClause.begin(), newClause.end() );

	uint32_t target( 0 );
	
	for (size_t i = 1; i != newClause.size(); ++i )
	  {
		if (newClause[i] == (newClause[target]^1) )
		  {
			/* tautologic clause */
			newClause.clear();
			return;
		  }
		/* skip duplicate literals */
		else if (newClause[i] != newClause[target] )
		  {
			++target;
			if (target != i )
			  { 
				newClause[target] = newClause[i]; 
			  }
		  }
	  }

	// Cut off last literals
	if ((target + 1) < newClause.size() )
	  { 
		newClause.resize( target + 1 ); 
	  }    
  }

  // Merges a n-nary clauses "clause" with a binary (reslit + otherlit) and common literal "reslit"
  // Store new clause in "newClause"
  void Preprocessor::MergeClauses(Clause& clause, uint32_t reslit, uint32_t otherlit, std::vector< uint32_t >& newClause) const
  {
	newClause.clear();
	newClause.reserve( clause.size() );
  
	assert( !_assignment[otherlit] && !_assignment[otherlit^1] );

	bool pushed( false );
	uint32_t size( clause.size() );

	for (uint32_t i = 0; i != size; ++i )
	  {
		assert( !_assignment[clause[i]] && !_assignment[clause[i]^1] );

		// Tautological clause?
		if (clause[i] == (otherlit^1) )
		  { newClause.clear(); return; }
		// Put "otherlit" into new clause?
		else if (!pushed && ( clause[i] > otherlit ) )
		  { 
			newClause.push_back(otherlit); 
			pushed = true;
		  }

		// Skip resolvent
		if (clause[i] != (reslit^1) ) 
		  { 
			if (clause[i] != otherlit )
			  { newClause.push_back(clause[i]); }
		  }
	  }

	// Reached last element without adding "otherlit"? 
	// Then push "otherlit" add end of new clause
	if (!pushed )
	  { newClause.push_back(otherlit); }
  }

  // Counts the literals of merge of clause "clause1" and "clause2"
  // Return 0, if result is a tautology
  uint32_t Preprocessor::CountMergeClauses(const Clause& clause1, const Clause& clause2, uint32_t reslit) const
  {
	uint32_t negreslit( reslit^1 );

	uint32_t i( 0 );
	uint32_t j( 0 );
	uint32_t resultsize( 0 );
	uint32_t skipped( 0 );
	  
	uint32_t c1size(clause1.size());
	uint32_t c2size(clause2.size());

	while( i < c1size && j < c2size )
	  {
		if(clause1[i] == 0 || clause2[j] == 0)
		  {
			std::cout << "ooopsi " << i << " "  << j << std::endl;
			clause1.Print();
			clause2.Print();
		  }
		assert (clause1[i] != 0);
		assert (clause2[j] != 0);
		// skip resolvent literal
		if ((clause1[i] == reslit) )
		  { ++i; ++skipped; continue;}
		if ((clause2[j] == negreslit) )
		  { ++j; ++skipped; continue;}

		// tautological clause
		if (clause1[i] == (clause2[j]^1) )
		  { return 0; } 
		else if (clause1[i] < clause2[j] )
		  { ++i; ++resultsize; }
		else if (clause1[i] > clause2[j] )
		  { ++j; ++resultsize; }
		// count duplicated literal only once
		else
		  { ++i; ++j; ++resultsize; }
	  }
	
	assert( skipped >0 );
	// If both resolvents are visited, "(c1size-i) + (c2size-j)" is the correct size.
	// If the resolvent in one clause is not visited, we have to adjust the size.
	resultsize += ((c1size-i) + (c2size-j) - (2-skipped) );
	
	return resultsize;
  }

  // Counts the literals of merge of clause "clause" and "(reslit + otherlit)"
  // Return 0, if result is a tautology
  uint32_t Preprocessor::CountMergeClauses(const Clause& clause, uint32_t reslit, uint32_t otherlit) const
  {
	uint32_t i(0);
	// Count from "1" for "otherlit"
	int32_t resultsize(1);

	assert( !clause.ToDelete() );

	uint32_t csize(clause.size());

	while( i < csize )
	  {
		if(clause[i] == 0 )
		  {
			std::cout << "ooopsi " << i << std::endl;
			clause.Print();
		  }
		assert (clause[i] != 0);
		// skip resolvent literal
		if ((clause[i] == reslit) )
		  { ++i; ++resultsize; }
		// tautological clause
		else if (clause[i] == (otherlit^1) )
		  { return 0; } 
		else if (clause[i] < otherlit )
		  { ++i; ++resultsize; }
		else if (clause[i] > otherlit )
		  { break; }
		// Skip duplicated literal
		else
		  { ++i; }
	  }

	resultsize += (csize-i)-1;

	return resultsize;
  }

  // Estimates costs for variable elimination of "var"
  int32_t Preprocessor::EstimateCosts(uint32_t var) 
  {
	int32_t opv = (int)OccurenceCount(var<<1);
	int32_t onv = (int)OccurenceCount((var<<1)^1);

	if (opv == 0 && onv == 0 ) 
	  {
		// Unused variable

		// Return some incredibly high cost 
		return 1000000000;
	  }
    
	int32_t spv = ( (int)_binaries[var<<1].size() )<<1;

	for (uint32_t i = 0; i != _occur[var<<1].size(); ++i )
	  { spv += _ca[_occur[var<<1][i]].size(); }

	int32_t snv = ( (int)_binaries[(var<<1)^1].size() )<<1;

	for (uint32_t i = 0; i != _occur[(var<<1)^1].size(); ++i )
	  { snv += _ca[_occur[(var<<1)^1][i]].size(); }

	assert( opv + onv > 0 );  
	assert( spv >= (opv<<1) );
	assert( snv >= (onv<<1) );
      
	int32_t cost = onv * (spv - opv ) + opv * ( snv - onv ) - spv - snv;
    
	return cost;
  }

  // Estimate costs for variable elimination of "var" with more accurate "countMergeClauses"
  int32_t Preprocessor::EstimateCosts2(uint32_t var) const 
  {
	assert( var <= _variables );
	uint32_t poslit(var<<1);
	uint32_t neglit((var<<1)^1);
	size_t posbinsize(_binaries[poslit].size());
	size_t negbinsize(_binaries[neglit].size());
	size_t posoccursize(_occur[poslit].size());
	size_t negoccursize(_occur[neglit].size());

	// Return ridiculously high cost if variable occurs more than 20 times 
	// I.e. do not take this variable into account
	if ((posbinsize + negbinsize + posoccursize + negoccursize) > 20 )
	  { return 10000000; }

	// decrement binary lits 
	int32_t totalliterals = -static_cast<int32_t>(( _binaries[poslit].size() + _binaries[neglit].size() )<<1);
	
	for (uint32_t i = 0; i != posoccursize; ++i )
	  { 
		Clause& clausepos( _ca[_occur[poslit][i]] );
		assert( !clausepos.ToDelete()); 

		// TODO un/count learnt clauses
		// decrement pos n-nary lits
		totalliterals -= clausepos.size(); 

		if (!clausepos.IsLearned() )
		  { 
			// increment neg binary, pos n-nary lits
			for (uint32_t j = 0; j != negbinsize; ++j )
			  {
				const Watcher& negbin( _binaries[neglit][j]);
				if (!negbin.IsLearnedBinary() )
				  {
					uint32_t tmpcount = CountMergeClauses( clausepos, poslit, negbin.GetSecondLit() );
					totalliterals += (int)tmpcount;
				  }
			  }

			// increment pos n-nary, neg n-nary lits
			for (uint32_t k = 0; k != negoccursize; ++k )
			  {
				Clause& clauseneg( _ca[_occur[neglit][k]] );
				if (!clauseneg.IsLearned() )
				  {
					uint32_t tmpcount = CountMergeClauses( clausepos, clauseneg, poslit );
					totalliterals += (int)tmpcount;
				  }
			  }
		  }
	  }

	for (uint32_t i = 0; i != negoccursize; ++i )
	  { 
		Clause& clause( _ca[_occur[neglit][i]] );
		assert( !clause.ToDelete() ); 

		// decrement neg n-nary lits
		totalliterals -= clause.size();

		if (!clause.IsLearned() )
		  { 
			// increment pos binary, neg n-nary lits
			for (uint32_t j = 0; j != posbinsize; ++j )
			  {
				const Watcher& posbin(_binaries[poslit][j]);
				if (!posbin.IsLearnedBinary() )
				  {
					uint32_t tmpcount = CountMergeClauses( clause, neglit, posbin.GetSecondLit() );
					totalliterals += (int)tmpcount;
				  }
			  }
		  }
	  }

	for (uint32_t i = 0; i != posbinsize; ++i )
	  {
		const Watcher& posbin(_binaries[poslit][i]);
		for (uint32_t j = 0; j != negbinsize; ++j )
		  {
			const Watcher& negbin( _binaries[neglit][j]);
			// increment pos+neg binary lits
			if (!posbin.IsLearnedBinary() && !negbin.IsLearnedBinary() && (posbin.GetSecondLit() != (negbin.GetSecondLit()^1)) )
			  { 
				totalliterals += 2;  
			  }
 
		  }
	  }
	return totalliterals;
  }

  // Estimate costs for variable elimination of "var" with more accurate "countMergeClauses"
  // This method only counts literals from static clauses
  void Preprocessor::EstimateCosts3(VarCandidate& varcand) const 
  {
	uint32_t var( varcand.GetVariable() );
	assert( var <= _variables );
	uint32_t poslit(var<<1);
	uint32_t neglit((var<<1)^1);
	size_t posbinsize(_binaries[poslit].size());
	size_t negbinsize(_binaries[neglit].size());
	size_t posoccursize(_occur[poslit].size());
	size_t negoccursize(_occur[neglit].size());

	int32_t totalliterals(0);

	// Return ridiculously high cost if variable occurs more than 20 times 
	// I.e. do not take this variable into account
	
	if ((posbinsize + negbinsize + posoccursize + negoccursize) > 20 )
	  { 
		varcand.SetTooLarge(true);
		return; 
	  }

	for (uint32_t i = 0; i != posoccursize; ++i )
	  { 
		Clause& clausepos( _ca[_occur[poslit][i]] );
		assert( !clausepos.ToDelete()); 

		if (!clausepos.IsLearned() )
		  { 
			// decrement pos n-nary lits
			totalliterals -= clausepos.size(); 
			 

			// increment neg binary, pos n-nary lits
			for (uint32_t j = 0; j != negbinsize; ++j )
			  {
				const Watcher& negbin( _binaries[neglit][j]);
				if (!negbin.IsLearnedBinary() )
				  {
					uint32_t tmpcount = CountMergeClauses( clausepos, poslit, negbin.GetSecondLit() );
					totalliterals += (int)tmpcount;
				  }
			  }

			// increment pos n-nary, neg n-nary lits
			for (uint32_t k = 0; k != negoccursize; ++k )
			  {
				Clause& clauseneg( _ca[_occur[neglit][k]] );
				if (!clauseneg.IsLearned() )
				  {
					uint32_t tmpcount = CountMergeClauses( clausepos, clauseneg, poslit );
					totalliterals += (int)tmpcount;
				  }
			  }
		  }
	  }

	for (uint32_t i = 0; i != negoccursize; ++i )
	  { 
		const Clause& clause( _ca[_occur[neglit][i]] );
		assert( !clause.ToDelete() ); 

		if (!clause.IsLearned() )
		  { 
			// decrement neg n-nary lits
			totalliterals -= clause.size();
			// increment pos binary, neg n-nary lits
			for (uint32_t j = 0; j != posbinsize; ++j )
			  {
				const Watcher& posbin(_binaries[poslit][j]);
				if (!posbin.IsLearnedBinary() )
				  {
					uint32_t tmpcount = CountMergeClauses( clause, neglit, posbin.GetSecondLit() );
					totalliterals += (int)tmpcount;
				  }
			  }
		  }
	  }

	for (uint32_t i = 0; i != posbinsize; ++i )
	  {
		const Watcher& posbin(_binaries[poslit][i]);
		
		if (!posbin.IsLearnedBinary() )
		  {
			totalliterals -= 2;

			for (uint32_t j = 0; j != negbinsize; ++j )
			  {
				const Watcher& negbin( _binaries[neglit][j]);
				// increment pos+neg binary lits
				if (!negbin.IsLearnedBinary() && (posbin.GetSecondLit() != (negbin.GetSecondLit()^1)) )
				  { 
					totalliterals += 2;  
				  }
				
			  }
		  }
	  }

	for (uint32_t i = 0; i != negbinsize; ++i )
	  {
		const Watcher& negbin(_binaries[neglit][i]);
		if (!negbin.IsLearnedBinary() )
		  {
			totalliterals -= 2;
		  }
	  }

	//	std::cout << __func__ << var << " " << totalliterals << std::endl;

	varcand.SetRecalc(false);
	varcand.SetCosts(totalliterals);
  }

  // Adds a clause to database
  // Do _not_ update watch lists
  // A binary clause is put into the binary datastructure of preprocessor
  // Updates occurence lists for n-nary clauses
  // This method is only used in the preprocessor 
  bool Preprocessor::AddClausePrepro(std::vector<uint32_t>& clause, bool updateWatches, uint32_t lbd)
  {
	assert( !_emptyClause );

	// Are we really on decision level 0?
	assert(_decisionLevel == 0); 

	// If "clause" is empty, we might have a problem.
	assert(!clause.empty()); 

	/*
	std::cout << __func__ << " ";
	for (uint32_t c = 0; c < clause.size(); ++c)
	  {
		std::cout << Lit(clause[c]) << " [" << isAssigned(clause[c]) << "] ";
	  }
	std::cout << std::endl;
	*/
		
	// Assumes that clause is already sorted

	// "Shifted" variable indices are assumed to be greater 1.
	assert(clause.front() > 1); 
	assert( (clause.back() >> 1) <= _variables );

	// Initialization.
	size_t stop(clause.size());
	uint32_t lit(0);
	uint32_t size(0); 

	// Check whether "clause" is already satisfied or represents a tautological clause.  
	// By the way, search for multiple copies of the same literal and literals evaluating to FALSE. 
	for (size_t c = 0; c < stop; ++c)
	  {
		// Get the next literal.
		uint32_t l(clause[c]);

		// "clause" satisfied by "l"? Do we have a tautological clause?
		if (_assignment[l] || (l ^ 1) == lit)
		  { clause.resize(0); return true; }

		// Do we have to take the current literal into account?
		if (!_assignment[l ^ 1] && l != lit)
		  { clause[size++] = l; lit = l; }
	  }

	// Do we have an empty clause? CNF formula unsatisfiable?
	if (size == 0)
	  { _emptyClause = true; return false; }
      
	// Resize "clause" (necessary for the multi-threaded mode to work correctly).
	clause.resize(size); 

	// Do we have a unit clause?
	if (size == 1)
	  {
		// Push the unit literal as a "fake" implication onto the decision stack.
		_core->AddImplication(clause[0]);

		// Everything went fine.
		return true; 
	  }

	// Do we have a binary clause?
	if (size == 2)
	  {	 
		uint32_t wl0(clause[0]);
		uint32_t wl1(clause[1]); 
		_core->IncreaseActivity(wl0>>1);
		_core->IncreaseActivity(wl1>>1);

		AddBinary( wl0, wl1, false );

		// Everything went fine.
		return true; 
	  }
	else if (size == 3 )
	  {
		++_statistics.currentTernaryClauses;
	  }
	else
	  {
		++_statistics.currentNaryClauses;
	  }

	// Initialization.
	CRef cr(_ca.Alloc(clause, lbd, _core->_activityInc, size));

	_ca[cr].CalcSign();
	
	// update occurence lists
	for (uint32_t l = 0; l < size; ++l)
	  { 
		_core->IncreaseActivity(clause[l]>>1);
		_occur[clause[l]].push_back(cr); 
	  }

	// Finally, add the new clause to the clause database.
	// Update Watch lists if necessary
	if (updateWatches)
	  {
		_core->AttachClause(cr);
	  }
	else
	  {
		_clauseDatabase.push_back(cr);
	  }

	// Everything went fine.
	return true;
  }

  /* End: some helper functions */

  // Propagate new units
  // Delete all satisfied clauses and their occurence	
  // Strengthen Clauses
  // Returns "false" if a contradiction occurs (formula is UNSAT), otherwise true
  bool Preprocessor::PropagateUnits(void)
  {
	assert( _decisionLevel == 0 );

	if (!_okay )
	  { 
		return false;
	  }

	uint32_t units(0);

	// first deduce all collected implications
	while( _dsImplIndex != _dsEndIndex )
	  {
		// Remove every clause with "poslit" and strengthen every clause with "neglit"
		uint32_t poslit( _decisionStack[_dsImplIndex] );
		uint32_t neglit( _decisionStack[_dsImplIndex]^1 );

		_forcing[poslit>>1].ClearReason();

		_statistics.currentBinaryClauses -= static_cast<uint32_t>( _binaries[poslit].size() + _binaries[neglit].size() );

		++_dsImplIndex;
		
		// Has to be assigned, since it's part of the decision stack.
		assert( _assignment[poslit] || _assignment[neglit] ); 

		if (_binaries[poslit].empty() && _binaries[neglit].empty() && _occur[poslit].empty() && _occur[neglit].empty() )
		  {
			// Mark variable as deleted
			_deleted[poslit>>1] = true;
			_level[poslit>>1] = 0;
			continue; 
		  }

		++units;
		++_statistics.unitPropagations;
		--_statistics.usedVariables;

		//std::cout << __func__ << " " << helper::Lit(poslit) << std::endl;
		//PrintCompleteLists(poslit>>1);
		//_core->dumpCNF(true);

		// Proceed positive binary occurences
		for (size_t i = 0; i != _binaries[poslit].size(); ++i)
		  {
			uint32_t seclit( _binaries[poslit][i].GetSecondLit() );
			RemoveBinary( seclit, poslit );
		  }

		// Clear binary list
		_binaries[poslit].clear();

		// Proceed positive n-nary clauses
		while( !_occur[poslit].empty() ) 
		  {
			CRef cr = *_occur[poslit].rbegin();
			Clause& clause( _ca[cr] );
			_occur[poslit].pop_back();
				
			uint32_t csize(clause.size());

			if (csize == 3 )
			  {
				--_statistics.currentTernaryClauses;
			  }
			else
			  {
				--_statistics.currentNaryClauses;
			  }
			
			// Update all other occurence lists
			for (uint32_t pos = 0; pos != csize; ++pos )
			  {
				if (clause[pos] != poslit )
				  { RemoveFromOccurenceList( clause[pos], cr ); }
			  }
			// Mark clause as to deleted
			clause.MarkForDeletion();
		  }		
		
		// Now strengthen the clauses with negative occurences


		// First consider binary clauses
		size_t binsize(_binaries[neglit].size()); 
		for (size_t i = 0; i != binsize; ++i)
		  {
			uint32_t seclit( _binaries[neglit][i].GetSecondLit() );

			if (_assignment[seclit^1] )
			  { 
				if (_setting->verbosity > 1 )
				  {
					std::cout << "c contradiction in unit propagation with " << helper::Lit(seclit) << " -> UNSAT" << std::endl;
				  }
				_okay = false;
				return false; 
			  }
		
			if (!_assignment[seclit] )
			  { 
				_core->AddImplication(seclit); 
			  }

			RemoveBinary( seclit, neglit );
		  }
		
		// Clear binary list
		_binaries[neglit].clear();

		// Now consider n-nary clauses with negative occurences
		while( !_occur[neglit].empty() ) 
		  {
			CRef cr = *_occur[neglit].rbegin();
			Clause& clause( _ca[cr] );
			_occur[neglit].pop_back();

			// Clause already deleted
			if (clause.ToDelete() )
			  { continue; }

			StrengthenClause( cr, neglit, true );
		  }

		// Mark variable as deleted
		_deleted[poslit>>1] = true;
		_level[poslit>>1] = 0;
	  }
	
	// At this point all implications are performed 
	// And all occurence lists are up to date
	// Now, delete marked clauses from clausedatabase

	ClearClauseDatabase();

	if ((_setting->verbosity > 2) && ( units != 0 ) )
	  { 
		std::cout << "c unit propagations      : " << units << std::endl; 
	  } 
	return true;
  }

  // Detect trivial monotone/pure literals
  void Preprocessor::DetectMonotone(bool incremental)
  {
	uint32_t monolits =  _statistics.monotoneVariables;
	uint32_t dclits =  _statistics.dontcareVariables;

	bool doitagain( false );

	uint32_t low(1);
	if (incremental )
	  { low = _firstPreVarIndex; }

	do { 
	  doitagain = false;
	  for (uint32_t v = low; v <= _variables; ++v )
		{
		  // Skip deleted variables and variables which are already assigned
		  if (!IsUsed(v) || _donttouch[v] )
			{ continue; }

		  bool dc( false );
		  uint32_t lit( 0 );

		  // Positive occurences of "v"?
		  if (_binaries[v<<1].empty() && _occur[v<<1].empty() )
			{
			  lit = (v<<1)^1;
			  doitagain = true;
			  dc = true;
			} 

		  // Negative occurences of "v"?
		  if (_binaries[(v<<1)^1].empty() && _occur[(v<<1)^1].empty() )
			{
			  lit = v<<1;
			  // Neither positive nor negative occurence of "v"? 
			  // => "v" is not used anymore
			  if (dc && !_donttouch[v] )
				{
				  _deleted[v] = true;
				  lit = 0;
				  ++_statistics.dontcareVariables;
				}
			  doitagain = true;
			}
 
		  // "v" is monotone with polarity of "lit"
		  if (lit != 0 )
			{ 
			  _core->AddImplication( lit ); 
			  ++_statistics.monotoneVariables;
			}
		}

#ifndef NDEBUG 
	  bool res = PropagateUnits(); assert(res);
#else
	  PropagateUnits();
#endif
	} while( doitagain );

	if (_setting->verbosity > 1 && (_statistics.monotoneVariables != monolits) )
	  {
		std::cout << "c monotone lits          : " << ( _statistics.monotoneVariables - monolits ) << std::endl 
				  << "c dc lits                : " << ( _statistics.dontcareVariables - dclits ) << std::endl;
	  }
	return;
  }

  // Search within binaries for two cases:
  //  1. [(a + b) * (a + ~b)] => a is constant, imply a
  // 2a. [(a + ~b) * (~a + b )] => a and b are equivalent, replace all occurences of a with b (or b with a)
  // 2b. [(a + b) * (~a + ~b )] => a and ~b are equivalent, replace all occurences of a with ~b (or b with ~a)
  bool Preprocessor::FindBinaryConsAndEquiv(bool& didsomething)
  {
	_methodLimit = 1<<19;

	uint32_t cons( _statistics.constantVariables );
	uint32_t equiv( _statistics.equivalentVariables );

	CountBinaryOccurences();

	_binCandidates.clear();
	_binCandidates.resize((_variables+1)<<1);

	for (uint32_t v = 1; v <= _variables; ++v )
	  {
		if (!_deleted[v] && _occurCounts[v] > 0 ) 
		  { 
			_binCandidates.insert(v);
		  }
	  }

	while( !_binCandidates.empty() )
	  {
		uint32_t v = _binCandidates.top();

		uint32_t plit( v<<1 );
		uint32_t nlit( plit^1 );

	  STARTOVER_VAR:

		//std::cout << v << " " << _methodLimit << std::endl;

		if (_methodLimit < 0 )
		  { break; }

		size_t psize(_binaries[plit].size());
		size_t nsize(_binaries[nlit].size());

		if (_deleted[v] || (psize==0) || (nsize=0) ) 
		  { continue; }

		_methodLimit -= (psize*nsize);

		for (size_t p = 0; p < psize; ++p )
		  {
			uint32_t seclit( _binaries[plit][p].GetSecondLit() );

			for (size_t n = 0; n < nsize; ++n)
			  {
				// Do we have a constant?
				if (_binaries[nlit][n].GetSecondLit() == seclit )
				  {
					didsomething = true;
					++_statistics.constantVariables;

					// Variable should be unassigned, since we delete immediatly the clauses of implications
					assert( !_assignment[seclit] && !_assignment[seclit^1] );

					// Add constant to decisionstack
					_core->AddImplication(seclit);

					// Propagate the new implication for the constant variable
					// Eliminates automatically the defining binaries for the constant
					if (!PropagateUnits() )
					  { return false; }

					// Start over with this variable
					goto STARTOVER_VAR;
				  }
				// Do we have an equivalence?
				else if (_binaries[nlit][n].GetSecondLit() == (seclit^1) )
				  {
					uint32_t replacelit( nlit );

					// If "seclit" is don't touch, swap with "replacelit" so that the don't touch variable will not be removed

					if (_donttouch[seclit>>1] )
					  { 
						if (!_donttouch[replacelit>>1] )
						  { 
							replacelit = seclit;
							seclit = nlit;
						  }
						// Both variables don't touch? Do nothing...
						else
						  { continue; }
					  }

					didsomething = true;
					++_statistics.equivalentVariables;

					// Remove binary clauses defining the equivalence
					RemoveBinary(replacelit, seclit^1);
					RemoveBinary(replacelit^1, seclit);

					// Replace the occurence of "seclit" with "replacelit"
					if (!ReplaceVariable(seclit, replacelit ) )
					  { 
						if (_setting->verbosity > 1 )
						  { std::cout << "c contradiction after replacement of " << helper::Lit(seclit) << " with " << helper::Lit( nlit ) << " -> UNSAT" << std::endl; }
						return false; 
					  }

					// Update "psize" and "nsize"
					psize = _binaries[plit].size();
					nsize = _binaries[nlit].size();
					  
					// Start over with this variable
					goto STARTOVER_VAR;
				  }
			  }
		  }
	  }

	if (_setting->verbosity > 1 && ( ( _statistics.constantVariables != cons ) || (_statistics.equivalentVariables != equiv ) ) )
	  {
		std::cout << "c binary cons            : " << (_statistics.constantVariables - cons) << std::endl
				  << "c binary equiv           : " << (_statistics.equivalentVariables - equiv) << std::endl;
	  }
	
	return PropagateUnits();
  }

  bool Preprocessor::FullSubsumptionCheck(bool incremental)
  {
	if (_setting->verbosity > 2 )
	  { std::cout << "c Performing subsumption checks..." << std::endl; }

	struct rusage resources;
	getrusage(RUSAGE_SELF, &resources);
	double timeS = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	uint32_t tmpsub(_statistics.subsumptions);
	uint32_t tmpselfsub(_statistics.selfSubsumptions);

	_methodLimit = (1<<23);
	
	uint32_t low(1);
	if (incremental )
	  { low = _firstPreVarIndex; }

	for (uint32_t i = low; i <= _variables; ++i )
	  {
		if (_deleted[i] )
		  { continue; }

		uint32_t poslit = i<<1;
		uint32_t neglit = (i<<1)^1;

		// First check if the binaries subsume n-nary clauses
		size_t bsize( _binaries[poslit].size() );
		for (size_t j = 0; j != bsize; ++j )
		  {
			if (_methodLimit < 0 )
			  { break; }
			CheckBinSub( poslit, _binaries[poslit][j] );
		  }

		bsize = _binaries[neglit].size();
		for (size_t j = 0; j != bsize; ++j )
		  {
			if (_methodLimit < 0 )
			  { break; }
			CheckBinSub( neglit, _binaries[neglit][j] );
		  }
	  }

	uint32_t startindex = 0;
	if (incremental )
	  { startindex = _firstPreVarIndex; }

	// Now check for every n-nary, if the clause is self-subsuming or subsumed by another n-nary clause
	for (uint32_t c = startindex; c != _clauseDatabase.size(); ++c )
	  {
		if (_methodLimit < 0 )
		  { break; }
		CRef cr = _clauseDatabase[c];
		Clause& clause = _ca[cr];

		// Already marked as deleted?
		if (clause.ToDelete() )
		  { continue; }
		
		// Subsumption Check
		if (IsSubsumed( clause, clause.Sign(), clause.IsLearned() ) )
		  { 
			++_statistics.subsumptions;
			ClearAllOccurences( cr ); 
		  }
	  }

	getrusage(RUSAGE_SELF, &resources);
	double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	_statistics.runtime_subsumption += (timeC-timeS);

	if (_setting->verbosity > 1 )
	  {
		std::cout << "c runtime subsumption    : " << (timeC-timeS) << "s (" << _statistics.runtime_subsumption << "s)" << std::endl;
	  }

	if (( (_statistics.subsumptions - tmpsub) > 0 ) || ( (_statistics.selfSubsumptions - tmpselfsub) > 0 ) )
	  { 
		SetFlag(PF_SUBSUMPTION);
		
		if (_setting->verbosity > 1 )
		  {
			std::cout << "c subsumptions           : " << (_statistics.subsumptions - tmpsub) << std::endl
					  << "c selfsubsumptions       : " << (_statistics.selfSubsumptions - tmpselfsub) << std::endl;
		  }
	  }

	return PropagateUnits();
  }

  // Mark all n-nary clauses which are subsumed by the binary (i1 + i2)
  // Perform self-subsumption, with one binary clause
  void Preprocessor::CheckBinSub(uint32_t i1, Watcher& binaryclause)
  {
	uint32_t i2 = binaryclause.GetSecondLit();
	// Get Signature of binary clause
	//uint64_t sign( 0 );
	//sign |= 1ULL << (i1 % 64);
	//sign |= 1ULL << (i2 % 64);
	uint32_t sign( 0 );
	sign |= 1 << (i1 % 32);
	sign |= 1 << (i2 % 32);
	
	uint32_t checklit( i1 );
	uint32_t seclit( i2 );

	if (_occur[seclit].size() < _occur[checklit].size() )
	{
	  checklit = i2;
	  seclit = i1;
	}

	size_t size( _occur[checklit].size() );

	for (size_t i = 0; i != size; ++i)
	  {
		CRef cr(_occur[checklit][i]);
		Clause& clause(_ca[cr]);

		if (( ( sign & ~clause.Sign() ) != 0 ) || clause.ToDelete() )
		  { continue; }

		uint32_t csize(clause.size());

		_methodLimit -= csize;

		for (uint32_t pos = 0; pos != csize; ++pos )
		  {
			// "clause" is self-subsuming with ( i1 + i2 )
			if (( clause[pos] == (seclit^1) ) && !_donttouch[seclit>>1] )
			  { 
				// If we have deleted the n-nary clause, we have to update the status variables
				if (StrengthenClause( cr, seclit^1, false ) )
				  { --size; --i; }

				++_statistics.selfSubsumptions;
			  }
			// "clause" is subsumed by ( i1 + i2 )
			else if (clause[pos] == seclit )
			  { 
				++_statistics.subsumptions;

				for (uint32_t p = 0; p != clause.size(); ++p)
				  {
					if (clause[p] != checklit )
					  { 
						_methodLimit -= _occur[clause[p]].size();
						RemoveFromOccurenceList( clause[p], cr ); 
					  }
				  }

				if (clause.size() == 3 )
				  { --_statistics.currentTernaryClauses; }
				else
				  { --_statistics.currentNaryClauses; }
					
				clause.MarkForDeletion();

				// Is case a learned clause subsumes a static one, we have to change the role
				// of the learned clause to "static"
				if (binaryclause.IsLearnedBinary() && !clause.IsLearned() )
				  {
					// "binaryclause" is related to "i1"
					binaryclause.SetLearnedBinary(false);
					// search for corresponding binary of "i2"
					size_t j = 0;
					size_t i2size = _binaries[i2].size();
					_methodLimit -= i2size;
					for (; j != i2size; ++j )
					  {
						if (_binaries[i2][j].GetSecondLit() == i1 )
						  {
							_binaries[i2][j].SetLearnedBinary(false);
							break;
						  }
					  }
					assert( j != i2size );
				  }

				assert( size == _occur[checklit].size() );
				_occur[checklit][i--] = _occur[checklit][--size];
				_occur[checklit].pop_back();
				break;
			  }
		  }
	  }
  }

  template< class T>
  uint32_t Preprocessor::GetSubsumerCandidate(const T& c)
  {
	uint32_t candidate(c[0]);
	size_t costs( _occur[c[0]].size() );

	_methodLimit -= c.size();

	for (uint32_t i = 1; i != c.size(); ++i )
	  {
		if (_occur[c[i]].size() < costs )
		  {
			costs = _occur[c[i]].size();
			candidate = c[i];
		  }
	  }
	return candidate;	
  }

  // Return true if c1 subsumes c2
  template< class T >
  //bool Preprocessor::Subsumes(const Clause& c1, const T& c2, uint64_t c2sign)
  bool Preprocessor::Subsumes(const Clause& c1, const T& c2, uint32_t c2sign)
  {
	size_t c1size( c1.size() );
	//uint64_t c1sign( c1.Sign() );
	size_t c1sign( c1.Sign() );
	size_t c2size( c2.size() );

	assert( !c1.ToDelete() );

	/* skip if:
	   clause has same or larger size
	   clause's signature doesn't fit
	*/
	if (( c1size >= c2size ) ||
		( ( c1sign & ~c2sign ) != 0 )
		)
	  { return false; }
	
	size_t i1(0);
	size_t i2(0);

	_methodLimit -= (c1size+c2size);
	
	// Check if "c1" subsumes "c2"
	while( ( i1 != c1size ) && 
		   ( i2 != c2size ) )
	  {
		if (c1[i1] == c2[i2] )
		  {
			++i1;
			++i2;
		  }
		else if (c1[i1] > c2[i2] )
		  {
			++i2;
		  }
		else
		  {
			return false;
		  }
	  }
	
	if (i1 == c1size )
	  {
		return true;
	  }
	return false;
  }

  // Checks if "clause" is subsumed by another n-nary clause
  template<class T>
  //bool Preprocessor::IsSubsumed(const T& clause, uint64_t sign, bool islearned)
  bool Preprocessor::IsSubsumed(const T& clause, uint32_t sign, bool islearned)
  {
	uint32_t lit( GetSubsumerCandidate(clause) );
	
	const std::vector< CRef >& occur( _occur[lit] );
	size_t osize( occur.size() );

	_methodLimit -= osize;

	for (size_t j = 0; j != osize; ++j )
	  {
		Clause& checkClause(_ca[occur[j]]);

		if (Subsumes( checkClause, clause, sign ) )
		  {
			
			// In case a learned clause subsumes the new clause, we have to declare the learned clause as static
			if (checkClause.IsLearned() && !islearned )
			  {
				checkClause.SetLBD(1);
			  }
			return true;
		  }
	  }
	return false;
    }

  // Checks if "clause" is subsumed by another n-nary clause
  template<class T>
  //bool Preprocessor::IsSubsumedExcept(const T& clause, uint64_t sign, bool islearned, CRef except)
  bool Preprocessor::IsSubsumedExcept(const T& clause, uint32_t sign, bool islearned, CRef except)
  {
	uint32_t lit( GetSubsumerCandidate(clause) );
	
	const std::vector< CRef >& occur( _occur[lit] );
	size_t osize( occur.size() );

	_methodLimit -= osize;

	for (size_t j = 0; j != osize; ++j )
	  {
		if (occur[j] == except)
		  { continue; }

		Clause& checkClause(_ca[occur[j]]);

		if (Subsumes( checkClause, clause, sign ) )
		  {
			
			// In case a learned clause subsumes the new clause, we have to declare the learned clause as static
			if (checkClause.IsLearned() && !islearned )
			  {
				checkClause.SetLBD(1);
			  }
			return true;
		  }
	  }
	return false;
    }

  //bool Preprocessor::CheckAllBinaries(const std::vector< uint32_t >& clause, uint64_t sign)
  bool Preprocessor::CheckAllBinaries(const std::vector< uint32_t >& clause, uint32_t sign)
  {
	size_t csize( clause.size() );
	
	// Check binaries first
	for (size_t i = 0; i != csize; ++i )
	  {
		size_t bsize( _binaries[clause[i]].size() );
		_methodLimit -= bsize;
		for (size_t j = 0; j != bsize; ++j )
		  {
			uint32_t seclit( _binaries[clause[i]][j].GetSecondLit());
			//uint64_t binsign( ( 1ULL << (seclit % 64) ) | ( 1ULL << ( clause[i] % 64) ) );
			uint32_t binsign( ( 1 << (seclit % 32) ) | ( 1 << ( clause[i] % 32 ) ) );

			// subsumption candidate!
			if (( binsign & ~sign ) == 0 )
			  {
				uint32_t k( 0 );
				while( ( k < csize ) && ( clause[k] < seclit ) )
				  { ++k; }
				// Subsumption found!
				if (( k < csize ) && ( clause[k] == seclit ) )
				  { 
					// In case a learned clause subsumes the new clause, we have to declare the learned clause as static
					if (_binaries[clause[i]][j].IsLearnedBinary() )
					  {
						_binaries[clause[i]][j].SetLearnedBinary(false);
						size_t k = 0;
						size_t seclitsize = _binaries[seclit].size();
						for (; k != seclitsize; ++k )
						  {
							if (_binaries[seclit][k].GetSecondLit() == clause[i] )
							  {
								_binaries[seclit][k].SetLearnedBinary(false);
								break;
							  }
						  }
						assert( k != seclitsize );
					  }

					return true; 
				  }
			  }
		  }
	  }
	return false;
  }

	//bool Preprocessor::CheckAllBinariesExcept(const std::vector< uint32_t >& clause, uint64_t sign, uint32_t lit1, uint32_t lit2)
	bool Preprocessor::CheckAllBinariesExcept(const std::vector< uint32_t >& clause, uint32_t sign, uint32_t lit1, uint32_t lit2)
  {
	size_t csize( clause.size() );
	
	// Check binaries first
	for (size_t i = 0; i != csize; ++i )
	  {
		size_t bsize( _binaries[clause[i]].size() );
		_methodLimit -= bsize;
		for (size_t j = 0; j != bsize; ++j )
		  {
			uint32_t seclit( _binaries[clause[i]][j].GetSecondLit());

			if (((lit1 == clause[i]) && (lit2 == seclit)) || ((lit1 == seclit) && (lit2 == clause[i])))
			  { continue; }

			//uint64_t binsign( ( 1ULL << (seclit % 64) ) | ( 1ULL << ( clause[i] % 64) ) );
			uint32_t binsign( ( 1 << (seclit % 32) ) | ( 1 << ( clause[i] % 32) ) );

			// subsumption candidate!
			if (( binsign & ~sign ) == 0 )
			  {
				uint32_t k( 0 );
				while( ( k < csize ) && ( clause[k] < seclit ) )
				  { ++k; }
				// Subsumption found!
				if (( k < csize ) && ( clause[k] == seclit ) )
				  { 
					// In case a learned clause subsumes the new clause, we have to declare the learned clause as static
					if (_binaries[clause[i]][j].IsLearnedBinary() )
					  {
						_binaries[clause[i]][j].SetLearnedBinary(false);
						size_t k = 0;
						size_t seclitsize = _binaries[seclit].size();
						for (; k != seclitsize; ++k )
						  {
							if (_binaries[seclit][k].GetSecondLit() == clause[i] )
							  {
								_binaries[seclit][k].SetLearnedBinary(false);
								break;
							  }
						  }
						assert( k != seclitsize );
					  }

					return true; 
				  }
			  }
		  }
	  }
	return false;
  }

  // Is "clause" subsumed by some clause in database?
  // Assumes that "clause" is sorted
  //bool Preprocessor::IsForwardSubsumed(const std::vector< uint32_t >& clause, uint64_t sign)
	bool Preprocessor::IsForwardSubsumed(const std::vector< uint32_t >& clause, uint32_t sign)
  {
	// Check with binary clauses
	if (CheckAllBinaries(clause, sign) )
	  { return true; }

	// Check n-naries
	return IsSubsumed(clause, sign, false);
  }

  // Is "clause" subsumed by some clause in database?
  // Assumes that "clause" is sorted
  //bool Preprocessor::IsForwardSubsumedExcept(const std::vector< uint32_t >& clause, uint64_t sign, CRef except)
  bool Preprocessor::IsForwardSubsumedExcept(const std::vector< uint32_t >& clause, uint32_t sign, CRef except)
  {
	// Check with binary clauses
	if (CheckAllBinaries(clause, sign) )
	  { return true; }

	// Check n-naries
	return IsSubsumedExcept(clause, sign, false, except);
  }

  // Is "clause" subsumed by some clause in database?
  // Assumes that "clause" is sorted
  //bool Preprocessor::IsForwardSubsumedExcept(const std::vector< uint32_t >& clause, uint64_t sign, uint32_t lit1, uint32_t lit2)
	bool Preprocessor::IsForwardSubsumedExcept(const std::vector< uint32_t >& clause, uint32_t sign, uint32_t lit1, uint32_t lit2)
  {
	// Check with binary clauses
	if (CheckAllBinariesExcept(clause, sign, lit1, lit2) )
	  { return true; }

	// Check n-naries
	return IsSubsumed(clause, sign, false);
  }

  // Performs variable elimination
  bool Preprocessor::VarElimination(bool incremental)
  {
	if (_setting->verbosity > 2 )
	  { std::cout << "c Performing variable elimination..." << std::endl; }

	struct rusage resources;
	getrusage(RUSAGE_SELF, &resources);
	double timeS = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	uint32_t resvars(_statistics.resolvedVariables);
	uint32_t reslits(_statistics.resolvedLiterals);

	bool recalcall(false);
	bool restart(false);

	_okay = true;

	_methodLimit = 1<<20;

	//printDatabase();

	// Initialize varcandidates before first variable elimination
	if (_varCandidates[1].GetVariable() != 1 )
	  {
		for (uint32_t v = 1; v <= _variables; ++v )
		  {	_varCandidates[v] = VarCandidate(v); }
	  }

	std::list< VarCandidate* > varcost;

	uint32_t low(1);
	if (incremental )
	  { low = _firstPreVarIndex; }

	for (uint32_t v = low; v <= _variables; ++v )
	  {
		// Skip deleted and don't touch variables
		if (!_deleted[v] && !_donttouch[v] )
		  {
			--_methodLimit;
			EstimateCosts3( _varCandidates[v] );
			  
			varcost.push_back( &_varCandidates[v] );
		  }
	  }

	varcost.sort(VarCandidateSorter());
	
	while( true )
	  {

#ifndef NDEBUG
		uint32_t reslitsinternal(_statistics.resolvedLiterals);
#endif

		if (varcost.empty() )
		  { break; }

		if (_methodLimit < 0)
		  { break; }

		if (recalcall )
		  {
			recalcall = false;
			varcost.sort(VarCandidateSorter());
			std::list<VarCandidate*>::iterator s = varcost.end();
			--s;
			for (; s != varcost.end(); --s)
			  {
				if (!(*s)->ToRecalc() )
				  { break; }

				--_methodLimit;
				EstimateCosts3(**s);
				varcost.sort(VarCandidateSorter());
			  }
		  }

		VarCandidate nxtcandidate;

		for (std::list<VarCandidate*>::iterator s = varcost.begin(); s != varcost.end(); ++s )
		  {
			// found a candidate for var Elimination!
			uint32_t var( (*s)->GetVariable() );
			int32_t costs ( (*s)->GetCosts() );
			
			if (!_deleted[var] && !(*s)->IsTooLarge() && ( costs <= _setting->varIncrease ) )
			  { 
				if ((*s)->ToRecalc() )
				  {
					--_methodLimit;
					EstimateCosts3(**s);
					if ((*s)->IsTooLarge() || ( (*s)->GetCosts() > _setting->varIncrease ) )
					  {
						continue;
					  }
				  }
				nxtcandidate = **s;
				varcost.erase(s);
				break;
			  }
		  }

		uint32_t resvar(nxtcandidate.GetVariable());
		// No candidate found, end loop
		if (resvar == 0) 
		  { 
			if (restart )
			  {
				restart = false;
				recalcall = true;
				continue;
			  }
			break;
		  }

		restart = true;

		++_statistics.resolvedVariables;

		SetFlag(PF_VAR_ELIM);

		std::vector< std::vector< uint32_t > > resolvents;

		//std::cout << "eliminate: " << resvar << " " << nxtcandidate.GetCosts() << std::endl;
		//printCompleteLists(resvar);
		
		// Store some original clauses for model preservation
		_rebuilder->AddVarElimination(resvar);

		uint32_t poslit = resvar<<1;
		uint32_t neglit = (resvar<<1)^1;

		std::vector<uint32_t> newClause( 2 );

		_methodLimit -= _binaries[poslit].size()*_binaries[neglit].size();
		_methodLimit -= _binaries[poslit].size()*_occur[neglit].size();
		_methodLimit -= _binaries[neglit].size()*_binaries[poslit].size();
		_methodLimit -= _binaries[neglit].size()*_occur[poslit].size();

		/* Now start resolution steps */
		for (size_t i1 = 0; i1 != _binaries[poslit].size(); ++i1 )
		  {

			// Skip learned
			if (_binaries[poslit][i1].IsLearnedBinary())
			  { 
				//std::cout << "binresolvent is learned" << std::endl;
				continue; 
			  }

			uint32_t posseclit( _binaries[poslit][i1].GetSecondLit() );

			assert( !_assignment[posseclit] && !_assignment[posseclit^1] );

			// Resolve binary clauses of "poslit" and "neglit"
			for (size_t i2 = 0; i2 != _binaries[neglit].size(); ++i2 )
			  {
				// Skip learned
				if (_binaries[neglit][i2].IsLearnedBinary())
				  { 
					//std::cout << "binresolvent is learned" << std::endl;
					continue; 
				  }
				  
				uint32_t negseclit( _binaries[neglit][i2].GetSecondLit() );
				
				assert( !_assignment[negseclit] && !_assignment[negseclit^1] );

				// Resolved tautology? => Skip!
				if (posseclit == (negseclit^1) )
				  { continue; }
				// Resolved (a+b) and (~a+b)? => b is constant
				else if (posseclit == negseclit )
				  {
					newClause.resize( 1 );
					newClause[0] = posseclit;
					
					resolvents.push_back(newClause);
					newClause.resize( 2 );
				  }
				// Ordinary resolution step
				else
				  {
					newClause[0] = posseclit;
					newClause[1] = negseclit;

					resolvents.push_back(newClause);
					newClause.resize( 2 );
				  }
			  }

			// Resolve binary clauses of "poslit" and n-nary clauses of "neglit"
			for (size_t i2 = 0; i2 != _occur[neglit].size(); ++i2 )
			  {
				Clause& c2( _ca[_occur[neglit][i2]] );
				assert( !c2.ToDelete() );

				if (c2.IsLearned() )
				  { continue; }

				MergeClauses( c2, poslit, posseclit, newClause );
				
				if (!newClause.empty() )
				  { 
					resolvents.push_back(newClause); 
				  }
			  }
			newClause.resize( 2 );
		  }

		for (size_t i1 = 0; i1 != _occur[poslit].size(); ++i1 )
		  {

			Clause& c1( _ca[_occur[poslit][i1]] );
			assert( c1.Lbd() != 0 );

			// Skip learned
			if (c1.IsLearned() )
			  { continue; }

			// Resolve binary clauses of "neglit" and n-nary clauses of "poslit"
			size_t binsize(_binaries[neglit].size());
			for (size_t i2 = 0; i2 != binsize; ++i2 )
			  {
				// Skip learned
				if (_binaries[neglit][i2].IsLearnedBinary())
				  { 
					continue; 
				  }

				uint32_t negseclit( _binaries[neglit][i2].GetSecondLit() );
				assert( !_assignment[negseclit] && !_assignment[negseclit^1] );

				MergeClauses( c1, neglit, negseclit, newClause );

				if (!newClause.empty() )
				  { 
					resolvents.push_back(newClause); 
				  }
				newClause.resize( 2 );
			  }

			// Resolve n-nary clauses of "neglit" and "poslit"
			binsize = _occur[neglit].size();
			for (size_t i2 = 0; i2 != binsize; ++i2 )
			  {
				Clause& c2(_ca[_occur[neglit][i2]]);
				assert( !c2.ToDelete() );

				if (c2.IsLearned() )
				  { continue; }

				MergeClauses( c1, c2, poslit, newClause );

				if (!newClause.empty() )
				  { 
					resolvents.push_back(newClause); 
				  }
			  }
		  }

		// First delete resolved variables from database
		DeleteVariable(resvar);

		// Now add temporary saved clauses
		for (size_t c = 0; c != resolvents.size(); ++c )
		  { 
			// Special treatment for units and binary 
			if (resolvents[c].size() > 2 )
			  {
				// Calc signature
			    uint32_t sign( 0 );
				for (uint32_t s = 0; s != resolvents[c].size(); ++s )
				  { 
					//sign |= 1ULL << (resolvents[c][s] % 64); 
					sign |= 1 << (resolvents[c][s] % 32); 
					// Recalc costs for involved literals
					_varCandidates[resolvents[c][s]>>1].SetRecalc(true);

				  }

				std::sort(resolvents[c].begin(), resolvents[c].end());
					 
				// Clause is subsumed by existing clause?
				if (!IsForwardSubsumed( resolvents[c], sign ) )
				  {
					_statistics.resolvedLiterals -= static_cast<uint32_t>(resolvents[c].size());
					// Add clause without initilizing watchlists
					if (!AddClausePrepro(resolvents[c]) )
					  { 
						if (_setting->verbosity > 1 )
						  { 
							std::cout << "c contradiction with resolved clause -> UNSAT" << std::endl; 
						  }
						_okay = false;
						break;
					  }
				  }
			  }
			// resolvent is binary or unit... no subsumption check needed (is performed implicitly in "addClausePrepro")
			else 
			  { 
				if (_varCandidates[resolvents[c][0]>>1].GetCosts() <= _setting->varIncrease )
				  {
					_varCandidates[resolvents[c][0]>>1].SetRecalc(true);
				  }
				if (resolvents[c].size() > 1 )
				  { 
					_varCandidates[resolvents[c][1]>>1].SetRecalc(true);
				  }

				_statistics.resolvedLiterals -= static_cast<uint32_t>(resolvents[c].size());
				if (!AddClausePrepro(resolvents[c]) )
				  {
					if (_setting->verbosity > 1 )
					  { 
						std::cout << "c contradiction with resolved unit binary clause -> UNSAT" << std::endl; 
					  }
					_okay = false;
					break;
				  }
			  }
		  }

		assert( (int)reslitsinternal <= (int)_statistics.resolvedLiterals );

		if (!_okay )
		  { break; }

		// Did we add new unit(s)?
		// First propagate these units before next resolution step
		if (_dsImplIndex != _dsEndIndex )
		  { 
			recalcall = true;
			if (!PropagateUnits() )
			  { 
				_okay = false;
				break;
			  }
		  }
	  }

	//printDatabase();

	assert( reslits <= _statistics.resolvedLiterals );
	getrusage(RUSAGE_SELF, &resources);
	double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	_statistics.runtime_varElim += (timeC-timeS);

	if (_setting->verbosity > 1)
	  {
		std::cout << "c runtime varElimination : " << (timeC-timeS) << "s (" << _statistics.runtime_varElim << "s)" << std::endl;
		if (_statistics.resolvedVariables != resvars) 
		  {
			std::cout << "c variable eliminations  : " << (_statistics.resolvedVariables-resvars) << std::endl;
			std::cout << "c literal reduction elim : " << (_statistics.resolvedLiterals-reslits) << std::endl;
		  }
	  }
	return PropagateUnits();
  }

  bool Preprocessor::Vivify(void)
  {
	// TODO: sort clause
	// TODO: candidate lists...
	if (_setting->verbosity > 2 )
	  { std::cout << "c Performing Vivification..." << std::endl; }

	struct rusage resources;
	getrusage(RUSAGE_SELF, &resources);
	double timeS = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	_methodLimit = 1<<19;
	
	// Update watch clauses
	UpdateWatches(false);

	CountOccurences(true);

	VarHeap<helper::DescendingOrder, size_t>* candidates = new VarHeap<helper::DescendingOrder, size_t>( helper::DescendingOrder<size_t>(_occurCounts));

	uint32_t subsumptions(_statistics.vivifySubsumptions);
	uint32_t units(_statistics.vivifyUnits);
	uint32_t diff(_statistics.vivifyDiff);
	uint32_t tmpimplindex(_dsImplIndex);

	size_t size(_clauseDatabase.size());
	for( size_t i = 0; i != size; ++i )
	  {
		if( _methodLimit < 0 )
		  {
			break;
		  }
		CRef cr = _clauseDatabase[i];
		const Clause& clause(_ca[cr]);
		
		// Skip if we use special treatment of ternary clauses
		if( _setting->useTernary && clause.size() == 3 )
		  {
			continue;
		  }

		candidates->clear();
		candidates->resize((_variables+1)<<1);
						
		std::vector< uint32_t > assumptions;
		std::vector< uint32_t > newclause;
		newclause.reserve(clause.size());
		bool skipclause(false);
		uint32_t checksize(0);
		// Try to add subsumed clauses, leading to a conflict
		for (uint32_t j = 0; j != clause.size(); ++j)
		  {
			if ( _assignment[clause[j]] )
			  {
				skipclause = true;
				break;
			  }
			else if ( !_assignment[clause[j]^1] )
			  {
				candidates->insert(clause[j]);
				++checksize;
			  }
		  }

		if (skipclause)
		  {
			continue;
		  }

		// temporally remove watches
		RemoveWatch(clause[0], cr);
		RemoveWatch(clause[1], cr);

		while(!candidates->empty())
		  {
			uint32_t nxtAssumption(candidates->top());
			assumptions.push_back(nxtAssumption^1);
			_methodLimit -= _binaries[nxtAssumption].size();
		  }

		
		// formula is conflicting together with assumptions
		if (!_core->DeduceAssumptions(assumptions))
		  {
			newclause = _core->GetFailedAssumptions();
		  }
		
		// Add watches again if clause was not subsumed
		if (newclause.empty() || (newclause.size()>=checksize))
		  {
			_binaries[clause[0]].push_back(Watcher(cr,clause[1],NNARY));
			_binaries[clause[1]].push_back(Watcher(cr,clause[0],NNARY));
		  }
		// Otherwise clear occurencelists and add new clause to database
		else
		  {
			++_statistics.vivifySubsumptions;

			// Flip all assumptions
			for (uint32_t n = 0; n != newclause.size(); ++n)
			  {
				newclause[n] ^= 1;
			  }

			// std::cout << "add new clause: " << std::endl;
			// for (uint32_t foo = 0; foo != newclause.size(); ++foo)
			//   {
			//  	std::cout << helper::Lit(newclause[foo]) << " ";
			//   }
			// std::cout << std::endl;
			
			_statistics.vivifyDiff += (checksize-static_cast<uint32_t>(newclause.size()));
			size_t lbd = clause.Lbd();
			if (lbd > newclause.size())
			  {
				lbd = newclause.size();
			  }
			if (newclause.size() == 1)
			  {
				++_statistics.vivifyUnits;
			  }
			else
			  {
				// update occurence lists
				// in case of unit, occurences are cleared by "PropagateUnits"
				_methodLimit -= clause.size();
				ClearAllOccurences(cr);
			  }
			
			SetFlag(PF_VIVIFY);
			// Now add the new subsumped clause to database
			if (!AddClausePrepro(newclause, true, static_cast<uint32_t>(lbd)) )
			  {
				delete candidates;
				return false;
			  }
		  }
	  }

	// Remove n-nary and ternary from watch list
	for( uint32_t v = 1; v <= _variables; ++v )
	  {
		if( _binaries[v<<1].empty() && _binaries[(v<<1)^1].empty() )
		  {	continue; }

		for( uint32_t literal = (v<<1); literal < ((v<<1)+2); ++literal )
		  {	
			std::vector<Watcher>& watches( _binaries[literal] );

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

	_control->ResetDone();

	getrusage(RUSAGE_SELF, &resources);
	double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	_statistics.runtime_vivify += (timeC-timeS);

	if (_setting->verbosity > 1)
	  {
		std::cout << "c runtime vivify         : " << (timeC-timeS) << "s (" << _statistics.runtime_vivify << "s)" << std::endl;
		if (_statistics.vivifySubsumptions != subsumptions) 
		  {
			std::cout << "c vivify subsumptions    : " << (_statistics.vivifySubsumptions-subsumptions) << std::endl;
			std::cout << "c vivify units           : " << (_statistics.vivifyUnits-units) << std::endl;
			std::cout << "c vivify diff literals   : " << (_statistics.vivifyDiff-diff) << std::endl;
		  }
	  }
				  
	_dsImplIndex = tmpimplindex;
	delete candidates;
	return PropagateUnits();
  }

  void Preprocessor::SetFlag(PreproFlag flag)
  {
	_pFlag = _pFlag|flag;
  }

  void Preprocessor::UnsetFlag(PreproFlag flag)
  {
	_pFlag = _pFlag&(~flag);
  }

  bool Preprocessor::GetFlag(PreproFlag flag) const 
  {
	return _pFlag&flag;
  }
}
