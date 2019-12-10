/********************************************************************************************
solverstate.cpp -- Copyright (c) 2013-2017, Tobias Schubert, Sven Reimer

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

#include "antombase.h"
#include "core.h"

#include "helper.h"
#include "statistics.h"
#include "watcher.h"
#include "reason.h"

#include <cstring>

namespace antom 
{

  // TODO: also store/delete/save pre-/inprostatus
  // Deletes the status of the SAT solving core saved before by "saveStatus()". 
  void Core::DeleteStatus(void)
  {
	for( uint32_t v = 0; v != _noOfVarGroups; ++v )
	  { _varOrder[v]->deleteStatus(); }

	_solverState.clauseDatabase.clear();
	for (uint32_t v = 1; v <= _solverState.variables; ++v)
	  {
		std::vector<Watcher>().swap(_solverState.watches[v << 1]); 
		std::vector<Watcher>().swap(_solverState.watches[(v << 1) ^ 1]); 
	  }
	std::vector<uint32_t>().swap(_solverState.level);  
	std::vector<double>().swap(_solverState.activity); 
	std::vector<unsigned char>().swap(_solverState.polarity);
	std::vector<unsigned char>().swap(_solverState.assignment);    
	std::vector<uint32_t>().swap(_solverState.model);
	std::vector<uint32_t>().swap(_solverState.decisionStack);
	std::vector<uint32_t>().swap(_solverState.dl2ds);   
	std::vector<uint32_t>().swap(_solverState.varGroup);    
	std::vector<Reason>().swap(_solverState.forcing);
  }

  // Saves the current status of the SAT solving core. 
  // The following variables/vectors are not stored:
  // -- "_control",
  // -- "_id",
  // -- "_progress", 
  // -- "_binaryConflictingClause". 
  void Core::SaveStatus(void)
  {
	// First of all, reset what we have saved so far (wrt. to "watches" and "clauseDatabase"). 

	_solverState.clauseDatabase.clear();
	for (uint32_t v = 1; v <= _solverState.variables; ++v)
	  {
		std::vector<Watcher>().swap(_solverState.watches[v << 1]); 
		std::vector<Watcher>().swap(_solverState.watches[(v << 1) ^ 1]); 
	  }

	// Now, simplify the clause database.
	Simplify(false); 

	// Store the variable ordering.
	for( uint32_t v = 0; v != _noOfVarGroups; ++v )
	  { 
		_varOrder[v]->saveStatus(); 
		_solverState.modeDS[v] = _modeDS[v];
		_solverState.decVarSign[v] = _decVarSign[v];
	  }

	//std::memcpy(_solverState.modeDS.data(), _modeDS.data(), _modeDS.size() * sizeof(_modeDS[0]));

	// Store our status variables/flags.
	_solverState.stats                = _statistics;
	_solverState.basicOps             = _basicOps;
	_solverState.decisionLevel        = _decisionLevel;
	_solverState.dsEndIndex           = _dsEndIndex;
	_solverState.dsImplIndex          = _dsImplIndex;
	_solverState.incVarActivity       = _incVarActivity;
	_solverState.emptyClause          = _emptyClause;
	_solverState.startPtr             = _startPtr;
	_solverState.endPtr               = _startPtr;

	// Initialization.
	uint32_t max((_variables << 1) + 2);

	// Resize the various vectors.
	
	_solverState.modeDSForVariables.resize(_variables+1,CACHEANDTOGGLE);
	//std::memcpy(_solverState.modeDSForVariables.data(), _modeDSForVariables.data(), _modeDSForVariables.size() * sizeof(_modeDSForVariables[0]));
	_solverState.level.resize(_variables + 1, 0);
	//std::memcpy(_solverState.level.data(), _level.data(), _level.size() * sizeof(_level[0]));
	_solverState.activity.resize(_variables + 1, 0.0);
	//std::memcpy(_solverState.activity.data(), _activity.data(), _activity.size() * sizeof(_activity[0]));
	_solverState.polarity.resize(_variables + 1, true);
	//std::memcpy(_solverState.polarity.data(), _polarity.data(), _polarity.size() * sizeof(_polarity[0]));
	_solverState.model.resize(_variables + 1, 0);
	//std::memcpy(_solverState.model.data(), _model.data(), _model.size() * sizeof(_model[0]));
	_solverState.decisionStack.resize(_variables + 1, 0);
	//std::memcpy(_solverState.decisionStack.data(), _decisionStack.data(), _decisionStack.size() * sizeof(_decisionStack[0]));
	_solverState.dl2ds.resize(_variables + 1, 1);
	//std::memcpy(_solverState.dl2ds.data(), _dl2ds.data(), _dl2ds.size() * sizeof(_dl2ds[0]));
	_solverState.varGroup.resize(_variables + 1, 1);
	//std::memcpy(_solverState.varGroup.data(), _varGroup.data(), _varGroup.size() * sizeof(_varGroup[0]));
	_solverState.forcing.resize(_variables + 1);
	//std::memcpy(_solverState.forcing.data(), _forcing.data(), _forcing.size() * sizeof(_forcing[0]));
	_solverState.assignment.resize(max, false);
	//std::memcpy(_solverState.assignment.data(), _assignment.data(), _assignment.size() * sizeof(_assignment[0]));
	
	_solverState.watches.resize(max);

	// Store all vectors 1/2. 

	
	for (uint32_t v = 1; v <= _variables; ++v)
	  {
		_solverState.modeDSForVariables[v] = _modeDSForVariables[v];
		_solverState.level[v]              = _level[v]; 
		_solverState.activity[v]           = _activity[v]; 
		_solverState.polarity[v]           = _polarity[v]; 
		_solverState.model[v]              = _model[v]; 
		_solverState.decisionStack[v]      = _decisionStack[v];
		_solverState.dl2ds[v]              = _dl2ds[v]; 
		_solverState.varGroup[v]           = _varGroup[v]; 
		_solverState.forcing[v]            = _forcing[v]; 
	  }
	

	// Store all vectors 2/2.
	for (uint32_t p = 2; p < max; ++p)
	  { 
		uint32_t w(0); 
		_solverState.assignment[p] = _assignment[p]; 
		std::vector<Watcher>& watches(_watches[p]); 
		while (w < watches.size() && watches[w].IsBinary())
		  { _solverState.watches[p].push_back(watches[w]); ++w; }
	  }

	_solverState.ca.AllocRegion(_ca.size());
	// copy Core region to 
	_ca.CopyTo(_solverState.ca);
	size_t size(_clauseDatabase.size());
	for (size_t c = 0; c < size; ++c)
	  {
		_solverState.clauseDatabase.push_back(_clauseDatabase[c]);
	  }
  }

  // Restores the status of the SAT solving core saved before by "saveStatus()".
  void Core::RestoreStatus(void)
  {
	_clauseDatabase.clear();      
	for (uint32_t v = 1; v <= _variables; ++v)
	  {
		std::vector<Watcher>().swap(_watches[v << 1]); 
		std::vector<Watcher>().swap(_watches[(v << 1) ^ 1]); 
	  }

	// Restore the variable ordering.
	for( uint32_t v = 0; v != _noOfVarGroups; ++v )
	  { 
		_varOrder[v]->restoreStatus();
		_modeDS[v] = _solverState.modeDS[v];
		_decVarSign[v] = _solverState.decVarSign[v];
	  }

	//std::memcpy(_modeDS.data(), _solverState.modeDS.data(), _modeDS.size() * sizeof(_modeDS[0]));

	_statistics           = _solverState.stats;
	// Restore our status variables/flags.
	_basicOps             = _solverState.basicOps;
	_variables            = _solverState.variables;
	_decisionLevel        = _solverState.decisionLevel;
	_dsEndIndex           = _solverState.dsEndIndex;
	_dsImplIndex          = _solverState.dsImplIndex;
	_incVarActivity       = _solverState.incVarActivity;
	_emptyClause          = _solverState.emptyClause;
	_startPtr             = _solverState.startPtr;
	_endPtr               = _solverState.endPtr;

	// Initialization.
	uint32_t max((_variables << 1) + 2);

	// Resize the various vectors.
	_modeDSForVariables.resize(_variables+1,CACHEANDTOGGLE);
	//std::memcpy(_modeDSForVariables.data(), _solverState.modeDSForVariables.data(), _modeDSForVariables.size() * sizeof(_modeDSForVariables[0]));
	_level.resize(_variables + 1, 0);
	//std::memcpy(_level.data(), _solverState.level.data(), _level.size() * sizeof(_level[0]));
	_activity.resize(_variables + 1, 0.0);
	//std::memcpy(_activity.data(), _solverState.activity.data(), _activity.size() * sizeof(_activity[0]));
	_polarity.resize(_variables + 1, true);
	//std::memcpy(_polarity.data(), _solverState.polarity.data(), _polarity.size() * sizeof(_polarity[0]));
	_model.resize(_variables + 1, 0);
	//std::memcpy(_model.data(), _solverState.model.data(), _model.size() * sizeof(_model[0]));
	_decisionStack.resize(_variables + 1, 0);
	//std::memcpy(_decisionStack.data(), _solverState.decisionStack.data(), _decisionStack.size() * sizeof(_decisionStack[0]));
	_dl2ds.resize(_variables + 1, 1);
	//std::memcpy(_dl2ds.data(), _solverState.dl2ds.data(), _dl2ds.size() * sizeof(_dl2ds[0]));
	_varGroup.resize(_variables + 1, 0);
	//std::memcpy(_varGroup.data(), _solverState.varGroup.data(), _varGroup.size() * sizeof(_varGroup[0]));
	_forcing.resize(_variables + 1);
	//std::memcpy(_forcing.data(), _solverState.forcing.data(), _forcing.size() * sizeof(_forcing[0]));
	_assignment.resize(max, false);
	//	std::memcpy(_assignment.data(), _solverState.assignment.data(), _assignment.size() * sizeof(_assignment[0]));
	_watches.resize(max); 

	// Restore all vectors 1/2.
	for (uint32_t v = 1; v <= _variables; ++v)
	  {
		_modeDSForVariables[v]   = _solverState.modeDSForVariables[v];
		_level[v]                = _solverState.level[v]; 
		_activity[v]             = _solverState.activity[v]; 
		_polarity[v]             = _solverState.polarity[v]; 
		_model[v]                = _solverState.model[v]; 
		_decisionStack[v]        = _solverState.decisionStack[v];
		_dl2ds[v]                = _solverState.dl2ds[v]; 
		_varGroup[v]             = _solverState.varGroup[v]; 
		_forcing[v]              = _solverState.forcing[v]; 
	  }

	// Store all vectors 2/2.
	for (uint32_t p = 2; p < max; ++p)
	  { 
		_assignment[p] = _solverState.assignment[p]; 
		std::vector<Watcher>& watches(_solverState.watches[p]); 
		size_t wSize(watches.size());
		for (size_t w = 0; w < wSize; ++w)
		  { _watches[p].push_back(watches[w]); }
	  }

	// Restore the clause database.
	uint32_t size = (uint32_t)_solverState.clauseDatabase.size();
	for (uint32_t c = 0; c < size; ++c)
	  {
		_ca.AllocRegion(_solverState.ca.size());
		_solverState.ca.CopyTo(_ca);
		AttachClause(_solverState.clauseDatabase[c]);
	  }
  }

  // Resets the SAT solving core. The following status flags/variables remain untouched:
  // -- The SAT solving threads unique ID number: "_id".
  // -- The pointer to the "Control" object: "_control".
  // -- The number of variables for which memory has been reserved: "_variables".
  void Core::Reset(void)
  {
	InstanceReset();

	// Reset the variable order.
	// Note: do not touch _noOfVarGroups
	for (uint32_t v = 0; v != _noOfVarGroups; ++v )
	  { 
		_modeDS[v]     = CACHEANDTOGGLE;
	  }
	
	// Reset some more status flags.
	_setting->ResetCore();
  }

  void Core::InstanceReset(void)
  {
	// Reset the clause database.
	_ca.Reset();
	_clauseDatabase.clear();
	_candidates.clear();

	// Reset the variable order.
	// Note: do not touch _noOfVarGroups
	for (uint32_t v = 0; v != _noOfVarGroups; ++v )
	  { 
		_varOrder[v]->clear(); 
		_decVarSign[v] = false;
	  }

	if (_variables > 0)
	  {
		std::fill(_level.begin(), _level.begin()+_variables+1, 0);
		std::fill(_activity.begin(), _activity.begin()+_variables+1, 0.0);
		std::fill(_polarity.begin(), _polarity.begin()+_variables+1, true);
		std::fill(_varGroup.begin(), _varGroup.begin()+_variables+1, 0);
		std::fill(_modeDSForVariables.begin(), _modeDSForVariables.begin()+_variables+1, CACHEANDTOGGLE);
		std::fill(_decisionStack.begin(), _decisionStack.begin()+_variables+1, 0);
		std::fill(_dl2ds.begin(), _dl2ds.begin()+_variables+1, 1);
		std::fill(_model.begin(), _model.begin()+_variables+1, 0);
		std::fill(_deleted.begin(), _deleted.begin()+_variables+1, false);
		std::fill(_forcing.begin(), _forcing.begin()+_variables+1, Reason());
		std::fill(_assignment.begin(), _assignment.begin()+((_variables+1)<<1), false);
	  }

	// Reset all variable related data structures.
	for (uint32_t v = 1; v <= _variables; ++v)
	  {
		//_level[v]              = 0;
		//_activity[v]           = 0.0;
		//_polarity[v]           = true;
		//_varGroup[v]           = 0;
		//_modeDSForVariables[v] = CACHEANDTOGGLE;
		//_decisionStack[v]      = 0;
		//_dl2ds[v]              = 1;
		//_model[v]              = 0;
		//_deleted[v]            = false;
		//_forcing[v].ClearReason();
		
		uint32_t pLit(v << 1);
		uint32_t nLit((v << 1) ^ 1);
		
		//_assignment[pLit]      = false;
		//_assignment[nLit]      = false;

#ifdef CLEARWATCHES
		_watches[pLit].clear();
		_watches[nLit].clear();
#else
		std::vector<Watcher>().swap(_watches[pLit]); 
		std::vector<Watcher>().swap(_watches[nLit]);
#endif
	  }
	_statistics.ResetCore();

	// Reset some more status flags.
	_emptyClause     = false;
	_startPtr        = 1;
	_endPtr          = 1;
	_activityInc     = 1;
	_conflict.ClearReason();
	_newUnits        = 0;
	_incVarActivity  = 1.0;
	_decisionLevel   = 0;
	_dsEndIndex      = 1;
	_dsImplIndex     = 1;
	_conflictLiteral = 0;
	_basicOps        = 0; 
	// reset also "_variables", the capacity is stored in _capacity
	// The value for capacity should be not touched, since the data structure sizes are unchanged
	_variables       = 0;
	_learnedClausesLimit = 20000;
  }
}
