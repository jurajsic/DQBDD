 
/********************************************************************************************
antombase.cpp -- Copyright (c) 2013-2016, Tobias Schubert, Sven Reimer

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

// Include antom related headers.
#include "antombase.h"

#include "helper.h"
#include "core.h"
#include "preprocessor.h"
#include "control.h"

//#include "solverproxy.h"
//#include "circuit.h"

#include <sys/resource.h>

namespace antom
{
  // Some (hopefully) promising configurations. 
  const uint32_t configs(10); 
  const DecisionStrategy c_decisionStrategy[configs] = {     CACHEANDTOGGLE,     CACHE,     CACHEANDTOGGLE,     CACHEANDTOGGLE,     ALWAYSFALSE,     CACHE,     ALWAYSFALSE,    ALWAYSTRUE,     CACHEANDTOGGLE,    CACHE  };
  const RestartStrategy c_restartStrategy[configs]  = { LUBY, LUBY,  GLUCOSE, GLUCOSE, LUBY, LUBY, GLUCOSE, GLUCOSE, GLUCOSE, LUBY };
  const uint32_t c_unitFactor[configs]       = {     8,     5,     8,     5,     6,     8,     6,    10,    12,    9 }; 
  const double c_decayFactor[configs]            = {  1.05,  1.10,  1.05,  1.10,  1.20,  1.05,  1.50,  1.30,  1.15, 1.03 }; 
  const bool c_lhbr[configs]                     = {  true,  true,  true,  true, false, false,  true, false, false, true };

  // Constructor. 
  AntomBase::AntomBase() :
	_antomSetting(new Settings),
    _control(new Control(_antomSetting)),
    _core(1),
    _preprocessor(1),
	//_solver(NULL),
	//_circuit(new Circuit),
    _sID(0),
    _result(1),
	_globalPropertyTrigger(),
	_stacksize(0),
	_satSolverCalls(0),
	_lastModel()
#ifndef NDEBUG
	,
	_debugNetwork()
#endif
  { 
    // Consistency check.
    assert(_antomSetting->threads != 0); 

    // Set the number of threads to be used in parallel.
#ifdef PARALLEL
    omp_set_nu_antomSetting->threads(_antomSetting->threads);
#endif

	_core[0] = new Core(_control, _antomSetting, 0);
	_preprocessor[0] = new Preprocessor( _core[0] );
	_core[0]->SetPreprocessor(_preprocessor[0]);
	
	_globalPropertyTrigger.resize(1,0);
  }
 
  // Destructor.
  AntomBase::~AntomBase(void) 
  {
    for (uint32_t i = 0; i < _antomSetting->threads; ++i)
      {
		delete _core[i];
		delete _preprocessor[i];
      }

	delete _control;
	delete _antomSetting;
	//delete _circuit;
	//delete _solver;
  }

  void AntomBase::SetThreads(uint32_t threads)
  {
	if (threads <= _antomSetting->threads)
	  {
		return;
	  }

	_result.resize(threads);
	_core.resize(threads);
	_preprocessor.resize(threads);
	// Initialize the SAT solving cores.
    for (uint32_t i = _antomSetting->threads; i < threads; ++i)
      {
		_core[i] = new Core(_control, _antomSetting, i);
		_preprocessor[i] = new Preprocessor( _core[i] );
		_core[i]->SetPreprocessor(_preprocessor[i]);
	  }
  }

  // Returns the ID of the thread that was able to solve the CNF.
  uint32_t AntomBase::SolvingThread(void) const 
  { return _sID; } 

  // Returns the number of used variables
  uint32_t AntomBase::Variables(void) const 
  { return _core[_sID]->Variables(); }

  // Returns the current number of clauses within the clause database. 
  uint32_t AntomBase::Clauses(void) const 
  { return _core[_sID]->Clauses(); }

  // Returns the current number of literals within the clause database. 
  uint32_t AntomBase::Literals(void) const 
  { return _core[_sID]->Literals(); }

  // Returns the number of decisions made so far.
  uint32_t AntomBase::Decisions(void) const 
  { return _core[_sID]->Decisions(); }

  // Returns the number of BCP operations performed so far.
  uint32_t AntomBase::Bcps(void) const 
  { return _core[_sID]->Bcps(); }

  // Returns the number of implications found so far. 
  uint64_t AntomBase::Implications(void) const 
  { return _core[_sID]->Implications(); }

  // Returns the number of conflicts encountered so far.
  uint32_t AntomBase::Conflicts(void) const 
  { return _core[_sID]->Conflicts(); }

  // Returns the number of restarts performed so far.
  uint32_t AntomBase::Restarts(void) const 
  { return _core[_sID]->Restarts(); }

  uint32_t AntomBase::BlockedRestarts(void) const 
  { return _core[_sID]->BlockedRestarts(); }

  // Returns the number of database simplifications performed so far.
  uint32_t AntomBase::Simplifications(void) const 
  { return _core[_sID]->Simplifications(); }

  // Returns the number of binary clauses deduced due to "Lazy Hyper Binary Resolution".
  uint32_t AntomBase::Lhbr(void) const 
  { return _core[_sID]->Lhbr(); }

  uint32_t AntomBase::UsedVariables(void ) const 
  { return _core[_sID]->UsedVariables(); }

  uint32_t AntomBase::CurrentBinaryClauses(void) const 
  { return _core[_sID]->CurrentBinaryClauses(); }

  uint32_t AntomBase::CurrentTernaryClauses(void) const 
  { return _core[_sID]->CurrentTernaryClauses(); }

  uint32_t AntomBase::CurrentNaryClauses(void) const 
  { return _core[_sID]->CurrentNaryClauses(); }

  // Returns the number of unit clauses deduced due to conflict analysis.
  uint32_t AntomBase::LearntUnitClauses(void) const 
  { return _core[_sID]->LearntUnitClauses(); }

  // Returns the number of binary clauses deduced due to conflict analysis.
  uint32_t AntomBase::LearntBinaryClauses(void) const 
  { return _core[_sID]->LearntBinaryClauses(); }

  // Returns the number of ternary clauses deduced due to conflict analysis.
  uint32_t AntomBase::LearntTernaryClauses(void) const 
  { return _core[_sID]->LearntTernaryClauses(); }

  uint32_t AntomBase::MinimizedLiterals(void) const
  { return _core[_sID]->MinimizedLiterals(); }

  // Returns the number of synchronizations performed so far. 
  uint32_t AntomBase::Synchronizations(void) const 
  { return _core[_sID]->Synchronizations(); }

  // Returns the number of inprocessings steps during solving main routine
  uint32_t AntomBase::Inprocessings(void) const 
  { return _core[_sID]->Inprocessings(); }

  // Returns the average "Literals Blocks Distance" of all conflict clauses deduced so far.
  double AntomBase::AvgLBD(void) const 
  { return _core[_sID]->AvgLBD(); }

  // Returns the average length of all conflict clauses deduced so far.
  double AntomBase::AvgCCLength(void) const 
  { return _core[_sID]->AvgCCLength(); }
  // Returns the solver's average decision level before backtracking.
  double AntomBase::AvgDL(void) const 
  { return _core[_sID]->AvgDL(); }
  // Returns the average number of decision levels cleared during conflict analysis. 
  double AntomBase::AvgDLclearedCA(void) const 
  { return _core[_sID]->AvgDLclearedCA(); }
  // Returns the average number of variables getting unassigned during conflict analysis. 
  double AntomBase::AvgVarsUnassignedCA(void) const 
  { return _core[_sID]->AvgVarsUnassignedCA(); }

  uint32_t AntomBase::BinaryConstants(void) const 
  { return _preprocessor[_sID]->BinaryConstants(); }
  uint32_t AntomBase::BinaryEquivalences(void) const 
  { return _preprocessor[_sID]->BinaryEquivalences(); }
  uint32_t AntomBase::UplaConstants(void) const 
  { return _preprocessor[_sID]->UplaConstants(); }
  uint32_t AntomBase::UplaEquivalences(void) const 
  { return _preprocessor[_sID]->UplaEquivalences(); }
  uint32_t AntomBase::SatConstants(void) const 
  { return _preprocessor[_sID]->SatConstants(); }
  uint32_t AntomBase::VariableEliminations(void) const 
  { return _preprocessor[_sID]->VariableEliminations(); }
  uint32_t AntomBase::LiteralEliminations(void) const 
  { return _preprocessor[_sID]->LiteralEliminations(); }
  uint32_t AntomBase::BlockedClauses(void) const 
  { return _preprocessor[_sID]->BlockedClauses(); }
  uint32_t AntomBase::HiddenTautologies(void) const
  { return _preprocessor[_sID]->HiddenTautologies(); }
  uint32_t AntomBase::HiddenSubsumptions(void) const
  { return _preprocessor[_sID]->HiddenSubsumptions(); }
  uint32_t AntomBase::MonotoneVariables(void) const
  { return _preprocessor[_sID]->MonotoneVariables(); }
  uint32_t AntomBase::DcVariables(void) const
  { return _preprocessor[_sID]->DcVariables(); }
  uint32_t AntomBase::SubsumedClauses(void) const
  { return _preprocessor[_sID]->SubsumedClauses(); }
  uint32_t AntomBase::SelfsubsumedLiterals(void) const 
  { return _preprocessor[_sID]->SelfsubsumedLiterals(); }
  uint32_t AntomBase::BvaVariables(void) const
  { return _preprocessor[_sID]->BvaVariables(); }
  uint32_t AntomBase::BvaLiterals(void) const
  { return _preprocessor[_sID]->BvaLiterals(); }
  uint32_t AntomBase::VivifySubsumptions(void) const
  { return _preprocessor[_sID]->VivifySubsumptions(); }
  uint32_t AntomBase::VivifyUnits(void) const
  { return _preprocessor[_sID]->VivifyUnits(); }
  uint32_t AntomBase::VivifyDiff(void) const
  { return _preprocessor[_sID]->VivifyDiff(); }
  uint32_t AntomBase::UnitPropagations(void) const
  { return _preprocessor[_sID]->UnitPropagations(); }
  double AntomBase::RuntimeUPLA(void) const
  { return _preprocessor[_sID]->RuntimeUPLA(); }
  double AntomBase::RuntimeSubsumption(void) const
  { return _preprocessor[_sID]->RuntimeSubsumption(); }
  double AntomBase::RuntimeVarElim(void) const
  { return _preprocessor[_sID]->RuntimeVarElim(); }
  double AntomBase::RuntimeBCE(void) const
  { return _preprocessor[_sID]->RuntimeBCE(); }
  double AntomBase::RuntimeBVA(void) const
  { return _preprocessor[_sID]->RuntimeBVA(); }
  double AntomBase::RuntimeHTE(void) const
  { return _preprocessor[_sID]->RuntimeHTE(); }
  double AntomBase::RuntimeVivify(void) const
  { return _preprocessor[_sID]->RuntimeVivify(); }
  // uint32_t runtimeSATConst ( void ) const { return _preprocessor[_sID]->runtimeSATConst(); }
  double AntomBase::RuntimePrepro(void) const
  { return _preprocessor[_sID]->RuntimePrepro(); }

  // Returns a reference to either the satisfying variable assignment (after calling one of the 
  // "solve()" routines) or the set of currently assigned variables (after calling one of the 
  // two "deduceAssumptions()" routines). Example:
  // model[17] =  0 --> x17 = unassigned
  // model[17] = 35 --> x17 = false 
  // model[17] = 34 --> x17 = true
  // In case neither "solve()/maxSolve()" nor "deduceAssumptions()" has been called beforehand, the 
  // vector contains invalid data. 
  const std::vector<uint32_t>& AntomBase::Model(void) const
  { return _core[_sID]->Model(); }

  // Returns a fresh variable index not used so far.
  uint32_t AntomBase::NewVariable(void) 
  {
#ifndef NDEBUG

    // If the threads do not handle the same number of variables, we might have a problem.
    for (uint32_t t = 1; t < _antomSetting->threads; ++t)
      { assert(_core[t]->Variables() == _core[t - 1]->Variables()); }

#endif

    // Get a new variable index.
    uint32_t index(_core[_sID]->Variables() + 1);
	
    // Set the maximum variable index of all threads to "index". 
	
    for (uint32_t t = 0; t < _antomSetting->threads; ++t)
      { _core[t]->SetMaxIndex(index); }

	if (index >= _lastModel.size())
	  {
		_lastModel.resize(index+1,0);
	  }
    
    // Return "index".
    return index; 
  }

  // Sets the maximum variable index to "max". 
  void AntomBase::SetMaxIndex(uint32_t max) 
  { 
	uint32_t maxindex(max);
	for (uint32_t t = 0; t < _antomSetting->threads; ++t) 
	  {
		_core[t]->SetMaxIndex(maxindex); 
	  }
  }

  // Sets the group of variable "var" to "grp". See "setDecisionStrategy()" for more details. 
  void AntomBase::SetVarGroup(uint32_t var, uint32_t grp)
  { 
	for (uint32_t t = 0; t < _antomSetting->threads; ++t) 
	  { _core[t]->SetVarGroup(var, grp); } 
  }

  void AntomBase::UseTernaryClauses(bool val)
  {
	_antomSetting->useTernary = val;
  }

  // Returns whether variable was already deleted in pre-/in-processing
  bool AntomBase::IsDeleted(uint32_t var) const 
  { return _core[_sID]->IsDeleted(var); }

  void AntomBase::SetCPULimit(double t)
  { _antomSetting->cpuLimit = t; }

  void AntomBase::SetMemoryLimit(uint32_t m)
  { _antomSetting->memLimit = m*1024*1024; }

  // Adds a clause to the clause databases of all threads. Returns FALSE if the CNF formula is unsatisfiable, 
  // otherwise TRUE will be returned. Assumes that the solver is on decision level 0 and that "clause" is not 
  // empty. Furthermore, all literals have to be encoded as follows, having variable indices greater 0:
  //  x3 <--> (3 << 1)     = 6
  // -x3 <--> (3 << 1) + 1 = 7
  // All clauses inserted into the clause database using "addClause()" are assumed to belong to 
  // the original CNF formula (having a "Literals Blocks Distance" of 1). 
  // IN THE MULTI-THREADED MODE, "maxSetIndex()" HAS TO BE CALLED BEFOREHAND.
  bool AntomBase::AddClause(std::vector<uint32_t>& clause, uint32_t lbd) 
  {
    // Add the clause to thread with ID 0 to see whether the CNF is already unsatisfiable or not.
    if (!_core[0]->AddClause(clause, lbd))
      { return false; }

    // Clause already satisfied or a tautological one?
    if (clause.empty())
      { return true; }

	// Now, add the (modified) clause to the clause databases of the remaining threads.     
#ifdef PARALLEL
#pragma omp parallel for
    for (uint32_t t = 1; t < _antomSetting->threads; ++t)
      { _result[t] = _core[t]->addClause(clause, lbd); /* assert(_result[t]); */ }

    for (uint32_t t = 1; t < _antomSetting->threads; ++t)
      { if (!_result[t]) { return false; } }
#endif

    // Everything went fine. 
    return true; 
  }

  // Adds a clause to the clause databases of all threads. Returns FALSE if the CNF formula is unsatisfiable, 
  // otherwise TRUE will be returned. Assumes that the solver is on decision level 0 and that "lits != NULL" 
  // and "num > 0" holds. Furthermore, all literals have to be encoded as follows, having variable indices greater 0:
  //  x3 <--> (3 << 1)     = 6
  // -x3 <--> (3 << 1) + 1 = 7
  // All clauses inserted into the clause database using "addClause()" are assumed to belong to 
  // the original CNF formula (having a "Literals Blocks Distance" of 1). 
  // NOTE, THAT THIS VARIANT OF "addClause()" REQUIRES THAT
  // 1) THE CLAUSE TO BE ADDED DOES NOT CONTAIN MULTIPLE COPIES OF THE SAME LITERAL,
  // 2) THE CLAUSE TO BE ADDED IS NOT A TAUTOLOGICAL ONE, AND
  // 3) "maxSetIndex()" HAS BEEN CALLED BEFOREHAND.
  bool AntomBase::AddClause(uint32_t* lits, uint32_t num, uint32_t lbd)
  {
    // Add the current clause to all clauses databases.

#ifdef PARALLEL
#pragma omp parallel for 
    for (uint32_t t = 0; t < _antomSetting->threads; ++t)
      { _result[t] = _core[t]->addClause(lits, num, lbd); }

#ifndef NDEBUG

    // If the results are not all the same, we might have a problem.
    for (uint32_t t = 1; t < _antomSetting->threads; ++t)
      { assert(_result[t] == _result[t - 1]); }

    // Return "_result[0]".
    return _result[0]; 
#endif
#else
    // Add the clause to thread with ID 0 to see whether the CNF is already unsatisfiable or not.
    _result[0] = _core[0]->AddClause(lits, num, lbd);

    return _result[0]; 
#endif
  }

  // Add unit clause, using literal encoding as in "addClause()"
  bool AntomBase::AddUnit(uint32_t lit)
  {
	std::vector<uint32_t > unitclause;
	unitclause.push_back(lit);
	return AddClause(unitclause);
  }

  #if 0
  // Adds encoding of an AND-Gate (output = input1 * input2) to the solver
  bool AntomBase::AddAndGate(uint32_t input1, uint32_t input2, uint32_t output)
  {
	std::vector<uint32_t> inputs {input1, input2};
	return AddAndGate(inputs, output);
  }

  bool AntomBase::AddAndGate(const std::vector<uint32_t>& inputs, uint32_t output)
  {
	_circuit->AddGate(inputs, output, ANDGATE);
	return AddBaseGate(inputs, true, output, false);
  }
  
  // Adds encoding of an NAND-Gate (output = !(input1 * input2)) to the solver
  bool AntomBase::AddNandGate(uint32_t input1, uint32_t input2, uint32_t output)
  {
	std::vector<uint32_t> inputs {input1, input2};
	return AddNandGate(inputs, output);
  }

  bool AntomBase::AddNandGate(const std::vector<uint32_t>& inputs, uint32_t output)
  {
	_circuit->AddGate(inputs, output, NANDGATE);
	return AddBaseGate(inputs, true, output, true);
  }
  
  // Adds encoding of an OR-Gate (output = input1 + input2) to the solver
  bool AntomBase::AddOrGate(uint32_t input1, uint32_t input2, uint32_t output)
  {
	std::vector<uint32_t> inputs {input1, input2};
	return AddOrGate(inputs, output);
  }

  bool AntomBase::AddOrGate(const std::vector<uint32_t>& inputs, uint32_t output)
  {
	_circuit->AddGate(inputs, output, ORGATE);
	return AddBaseGate(inputs, false, output, true);
  }
  
  // Adds encoding of an NOR-Gate (output = !(input1 + input2)) to the solver
  bool AntomBase::AddNorGate(uint32_t input1, uint32_t input2, uint32_t output)
  {
	std::vector<uint32_t> inputs {input1, input2};
	return AddNorGate(inputs, output);
  }

  bool AntomBase::AddNorGate(const std::vector<uint32_t>& inputs, uint32_t output)
  {
	_circuit->AddGate(inputs, output, NORGATE);
	return AddBaseGate(inputs, false, output, false);
  }
  
  bool AntomBase::AddBaseGate(const std::vector<uint32_t>& inputs, bool neginput, uint32_t output, bool negoutput)
  {
	std::vector<uint32_t> clause(inputs.size()+1);
	std::vector<uint32_t> binClause(2);

	for (uint32_t i = 0; i != inputs.size(); ++i)
	  {
		clause[i] = inputs[i]^neginput;
		binClause[0] = inputs[i]^(!neginput);
		binClause[1] = output^(!negoutput);
		if (!AddClause(binClause))
		  { return false; }
	  }
	
	clause[inputs.size()] = output^negoutput;
	if (!AddClause(clause))
	  { return false; }

	return true;
  }
							   
  // Adds encoding of an XOR-Gate (output = input1 XOR input2) to the solver
  bool AntomBase::AddXorGate(uint32_t input1, uint32_t input2, uint32_t output)
  {
	std::vector<uint32_t> inputs {input1, input2};
	_circuit->AddGate(inputs, output, XORGATE);
	return AddBaseXorGate(input1, input2, output);
  }

  // Adds encoding of an XNOR-Gate (output = !(input1 XOR input2)) to the solver
  bool AntomBase::AddXnorGate(uint32_t input1, uint32_t input2, uint32_t output)
  {
	std::vector<uint32_t> inputs {input1, input2};
	_circuit->AddGate(inputs, output, XNORGATE);
	return AddBaseXorGate(input1, input2, output^1);
  }

  bool AntomBase::AddBaseXorGate(uint32_t input1, uint32_t input2, uint32_t output)
  {
	std::vector<uint32_t> clause(3);

	clause[0] = input1^1;
	clause[1] = input2;
	clause[2] = output;
	if (!AddClause(clause))
	  { return false; }

	clause[0] = input1;
	clause[1] = input2^1;
	clause[2] = output;
	if (!AddClause(clause))
	  { return false; }

	clause[0] = input1;
	clause[1] = input2;
	clause[2] = output^1;
	if (!AddClause(clause))
	  { return false; }

	clause[0] = input1^1;
	clause[1] = input2^1;
	clause[2] = output^1;
	if (!AddClause(clause))
	  { return false; }

	return true;
  }
  #endif

  // Collects data for performing a regression analysis afterwards (--> SATzilla-like SAT solving). 
  // NOTE, THAT "getDataRegressionAnalysis()" CAN BE USED IN SINGLE-THREADED MODE ONLY.
  void AntomBase::GetDataRegressionAnalysis(void)
  {
    // Check the number of threads. 
    assert(_antomSetting->threads == 1);

    // Initialization.
    std::vector<uint32_t> assumptions; 

    // Save antom's current status.
    _core[0]->SaveStatus();

    // Now, analyse all configurations.
    for (uint32_t c = 0; c < configs; ++c)
      {
		// Output.
		std::cout << "c checking configuration " << c << "..." << std::endl;

		// Initialize current configuration.	
		_core[0]->SetDecisionStrategy(c_decisionStrategy[c], 0); 

		// Start some SAT solving.
		uint32_t result(_core[0]->Solve(assumptions, 10)); 
	
		// CNF formula not already solved?
		if (result == ANTOM_UNKNOWN)
		  {
			// Get some stats. 
			uint32_t x0(_core[0]->Variables());
			uint32_t x1(_core[0]->Decisions()); 
			uint32_t x2(_core[0]->Bcps()); 
			uint32_t x3(_core[0]->Conflicts());
			uint32_t x4(_core[0]->Restarts()); 
			double x5(_core[0]->AvgLBD()); 
			double x6(_core[0]->AvgCCLength()); 
			double x7(_core[0]->AvgDL()); 
			uint32_t x8(_core[0]->LearntUnitClauses()); 
			uint32_t x9(_core[0]->LearntBinaryClauses()); 
			uint32_t x10(_core[0]->LearntTernaryClauses()); 
			double x11(_core[0]->AvgDLclearedCA()); 
			double x12(_core[0]->AvgVarsUnassignedCA()); 
			uint32_t x13(_core[0]->Clauses()); 
			unsigned long x14(_core[0]->Implications());
			double x15((double) x0 / x13);
			double x16((double) x13 / x0);
			double x17((double) x1 / x3);
			double x18((double) x14 / x1);
			double x19((double) x14 / x3);
			uint32_t x20(_core[0]->Literals()); 
			double x21(_core[0]->Progress()); 
	    
			// Restore status.
			_core[0]->RestoreStatus();

			// Re-initialize the configuration under consideration.
			_core[0]->SetDecisionStrategy(c_decisionStrategy[c], 0); 

			//Solve the CNF with a higher limit wrt. synchronizations.
			result = _core[0]->Solve(assumptions, 12000); 

			// Output.
			if (result != ANTOM_UNKNOWN)
			  {
				std::cout << "c stats, config " << c << ": " << x0 << " " << x1 << " " << x2 << " " << x3 << " " << x4 
						  << " " << x5 << " " << x6 << " " << x7 << " " << x8 << " " << x9 << " " << x10 << " " << x11 
						  << " " << x12  << " " << x13 << " " << x14 << " " << x15 << " " << x16 << " " << x17 << " "
						  << x18 << " " << x19 << " " << x20 << " " << x21 << " " << _core[0]->Synchronizations() << std::endl; 
			  }

			// Reset "_control".
			_control->ResetDone();	    

		  }

		// Restore status. 
		_core[0]->RestoreStatus();
      }
  }

  // Dumps cnf into std::cout
  void AntomBase::DumpCNF(bool printAssignment) const
  {
	std::cout << "p cnf " << _core[_sID]->Variables() << " " << _core[_sID]->Clauses() << std::endl;
	_core[_sID]->DumpCNF(printAssignment);
  }

  /* antom core interface */

  void AntomBase::TrivialAssignment(void) const 
  { _core[_sID]->TrivialAssignment(); }

  // Activates or deactivates "Lazy Hyper Binary Resolution".
  // val = TRUE  --> LHBR enabled (default).
  // val = FALSE --> LHBR disabled.
  void AntomBase::SetLHBR(bool val)
  {
	_antomSetting->lhbr = val;
  }

  // Sets the unit factor of both restart strategies -- Luby & glucose-like -- to "val" (default: 8). 
  // The unit factor directly corresponds to the interval between two restart operations. 
  void AntomBase::SetLuby(uint32_t val)
  {
	_antomSetting->lubyShift = val;
  }

  // Sets the decision strategy of group "group" to mode "val". 
  // Currently, there are four modes that differ wrt. the polarity of a decision variable:
  // 0 (default) --> Use the variable's cached polarity together with antom's "polarity toggling scheme". 
  // 1           --> Use the variable's cached polarity only.
  // 2           --> The polarity will be set to FALSE regardless of the cached value. 
  // 3           --> The polarity will be set to TRUE regardless of the cached value.
  // Furthermore, antom maintains two variable orderings: "group 0" and "group 1". As long as "group 0" is 
  // non-empty, variables from that group will be preferred to serve as decision variables. By default, all 
  // variables belong to "group 1".
  void AntomBase::SetDecisionStrategy(uint32_t val, uint32_t group)
  { 
	assert( val < 4 );
	DecisionStrategy strat(CACHEANDTOGGLE);
	switch (val)
	  {
	  case 0 :
		strat = CACHEANDTOGGLE;
		break;
	  case 1 :
		strat = CACHE;
		break;
	  case 2 :	
		strat = ALWAYSFALSE;
		break;
	  case 3 :	
		strat = ALWAYSTRUE;
		break;
	  default:
		assert(false);
		break;
	  }
   
	for (uint32_t t = 0; t < _antomSetting->threads; ++t) 
	  { 
		_core[t]->SetDecisionStrategy(strat, group); 
	  }
 
  }

  // Set initial polarity of the variable "var" to "pol"
  void AntomBase::SetPolarity(uint32_t var, bool pol)
  {
	for (uint32_t t = 0; t < _antomSetting->threads; ++t) 
	  { 
		_core[t]->SetPolarity(var,pol); 
	  } 
  }

  // Like "setDecisionStrategy()" for a specific variable instead of a group
  void AntomBase::SetDecisionStrategyForVariable(uint32_t val, uint32_t var)
  { 
	assert( val < 4 );
	DecisionStrategy strat(CACHEANDTOGGLE);
	switch (val)
	  {
	  case 0 :
		strat = CACHEANDTOGGLE;
		break;
	  case 1 :
		strat = CACHE;
		break;
	  case 2 :	
		strat = ALWAYSFALSE;
		break;
	  case 3 :	
		strat = ALWAYSTRUE;
		break;
	  default:
		assert(false);
		break;
	  }

	for (uint32_t t = 0; t < _antomSetting->threads; ++t) 
	  { 
		_core[t]->SetDecisionStrategyForVariable(strat, var); 
	  } 
  }

  // Sets the restart strategy to model "val":
  // 0 --> Luby (default).
  // 1 --> Glucose-like.       
  void AntomBase::SetRestartStrategy(RestartStrategy val)
  {
	_antomSetting->restartStrategy = val;
  }

  void AntomBase::SetSimplifyStrategy(SimplifyStrategy val)
  {
	_antomSetting->simplifyStrategy = val;
  }

  void AntomBase::SetSimplifyActivity(SimplifyActivity val)
  {
	_antomSetting->simplifyActivity = val;
  }

  // Sets the decay factor to "val" (default: 1.05). The decay factor is responsible 
  // for controlling how the variable activities evolve during the search process. 
  void AntomBase::SetDecayFactor(double val) 
  {
	_antomSetting->decayFactor = val;
  }

  // Performs unit propagation, taking the current CNF and the specified assumptions into account. Returns 
  // FALSE if a conflict occurs, otherwise the return value is TRUE. NOTE, THAT "deduceAssumptions()" CAN BE 
  // USED IN SINGLE-THREADED MODE ONLY.
  bool AntomBase::DeduceAssumptions(void) 
  { 
	std::vector<uint32_t> assumptions; 
	return DeduceAssumptions(assumptions); 
  }

  bool AntomBase::DeduceAssumptions(const std::vector<uint32_t>& assumptions)
  { 
	assert(_antomSetting->threads == 1); 
	_sID = 0; 
	return _core[0]->DeduceAssumptions(assumptions); 
  }

  std::vector<uint32_t > AntomBase::GetFailedAssumptions(void) const
  {
	return _core[0]->GetFailedAssumptions();
  }

  // Solves the current CNF formula using the most promising configuration of antom. The return values are SAT/UNSAT.
  // NOTE, THAT "solveSATzilla()" CAN BE USED IN SINGLE-THREADED MODE ONLY.
  uint32_t AntomBase::SolveSATzilla(void)
  {
    // Check the number of threads. 
    assert(_antomSetting->threads == 1);

    // Initialization.
    std::vector<uint32_t> assumptions; 
    uint32_t result(ANTOM_UNKNOWN); 
    double estimatedSync(0.0); 
    double bestSync(0.0);
    uint32_t bestConfig(1);
    const double intercept[configs] = {  1.331e+03, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const double p0[configs]        = {  5.751e-03, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const double p1[configs]        = {  1.513e-01, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; 
    const double p2[configs]        = { -1.443e-01, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; 
    const double p3[configs]        = {  6.215e-03, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const double p4[configs]        = { -2.212e+00, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; 
    const double p5[configs]        = {  6.837e+01, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; 
    const double p6[configs]        = { -5.974e+00, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; 
    const double p7[configs]        = {  6.625e-01, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; 
    const double p8[configs]        = { -3.101e+00, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; 
    const double p9[configs]        = {  2.448e-01, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; 
    const double p10[configs]       = {  4.939e-01, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; 
    const double p11[configs]       = {  2.017e+02, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; 
    const double p12[configs]       = {  3.228e-01, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; 
    const double p13[configs]       = { -6.741e-04, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const double p14[configs]       = { -2.868e-04, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const double p15[configs]       = {  7.254e+03, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const double p16[configs]       = {  8.786e+01, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const double p17[configs]       = { -2.389e+02, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const double p18[configs]       = {  3.342e+00, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const double p19[configs]       = { -1.866e-01, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const double p20[configs]       = {          0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const double p21[configs]       = {          0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
 
    // Save antom's current status.
    _core[0]->SaveStatus();

    // Now, analyse all configurations.
    for (uint32_t c = 0; c < configs; ++c)
      {
		// Initialize current configuration.	
		_core[0]->SetDecisionStrategy(c_decisionStrategy[c], 0); 

		// Start some SAT solving.
		result = _core[0]->Solve(assumptions, 10); 
	
		// CNF formula already solved?
		if (result != ANTOM_UNKNOWN)
		  { return result; }
		else
		  {
			// Get some stats. 
			uint32_t x0(_core[0]->Variables());
			uint32_t x1(_core[0]->Decisions()); 
			uint32_t x2(_core[0]->Bcps()); 
			uint32_t x3(_core[0]->Conflicts());
			uint32_t x4(_core[0]->Restarts()); 
			double x5(_core[0]->AvgLBD()); 
			double x6(_core[0]->AvgCCLength()); 
			double x7(_core[0]->AvgDL()); 
			uint32_t x8(_core[0]->LearntUnitClauses()); 
			uint32_t x9(_core[0]->LearntBinaryClauses()); 
			uint32_t x10(_core[0]->LearntTernaryClauses()); 
			double x11(_core[0]->AvgDLclearedCA()); 
			double x12(_core[0]->AvgVarsUnassignedCA()); 
			uint32_t x13(_core[0]->Clauses()); 
			unsigned long x14(_core[0]->Implications());
			double x15((double) x0 / x13);
			double x16((double) x13 / x0);
			double x17((double) x1 / x3);
			double x18((double) x14 / x1);
			double x19((double) x14 / x3);
			uint32_t x20(_core[0]->Literals()); 
			double x21(_core[0]->Progress()); 

			// Estimate the number of synchronizations required to solve the CNF.
			estimatedSync = intercept[0]; 
			estimatedSync += p0[0] * (double) x0; 
			estimatedSync += p1[0] * (double) x1; 
			estimatedSync += p2[0] * (double) x2; 
			estimatedSync += p3[0] * (double) x3; 
			estimatedSync += p4[0] * (double) x4; 
			estimatedSync += p5[0] * x5; 
			estimatedSync += p6[0] * x6; 
			estimatedSync += p7[0] * x7; 
			estimatedSync += p8[0] * (double) x8; 
			estimatedSync += p9[0] * (double) x9; 
			estimatedSync += p10[0] * (double) x10; 
			estimatedSync += p11[0] * x11; 
			estimatedSync += p12[0] * x12; 
			estimatedSync += p13[0] * (double) x13; 
			estimatedSync += p14[0] * (double) x14; 
			estimatedSync += p15[0] * x15; 
			estimatedSync += p16[0] * x16; 
			estimatedSync += p17[0] * x17; 
			estimatedSync += p18[0] * x18; 
			estimatedSync += p19[0] * x19; 
			estimatedSync += p20[0] * x20; 
			estimatedSync += p21[0] * x21; 

			if (estimatedSync < (double) _core[0]->Synchronizations())
			  { estimatedSync = (double) _core[0]->Synchronizations(); }

			// Output. 
			if (_antomSetting->verbosity > 0)
			  { std::cout << "c estimated #synchronizations required by configuration " << c << ": " << estimatedSync << std::endl; }

			// Have we found a more promising configuration?
			if (c == 0 || estimatedSync < bestSync)
			  { bestConfig = c; bestSync = estimatedSync; }	       

			// Restore antom's status saved at the beginning.
			_core[0]->RestoreStatus();
		  }
      }

    // Output.
    if (_antomSetting->verbosity > 0)
      { std::cout << "c most promising configuration: " << bestConfig << std::endl; }

    // Configure antom, using the most promising configuration found before.
    _core[0]->SetDecisionStrategy(c_decisionStrategy[bestConfig], 0); 

    // Solve the CNF formula.
    result = _core[0]->Solve(assumptions); 

    // Output.
    if (_antomSetting->verbosity > 0)
      {
		double error((double) _core[0]->Synchronizations() / bestSync);
		std::cout << "c estimated / real #synchronizations: " << error << std::endl;
      }

    // Return "result".
    return result;
  }
  
  // Simplify Database
  void AntomBase::Simplify(void) 
  {
	for (uint32_t t = 0; t < _antomSetting->threads; ++t)
      { _core[t]->Simplify(false); }
  }

  // Solves the current CNF formula, taking assumptions (if specified) into account. The assumptions have to be encoded in the 
  // same way as the literals of a clause (see "addClause()"). The return values are SAT/UNSAT. In the multi-threaded mode, the 
  // CNF formula gets solved by "_antomSetting->threads" threads running in parallel, according to the well-known algorithm portfolio scheme. 
  uint32_t AntomBase::Solve(void) 
  { 
	std::vector<uint32_t> assumptions; 
	return Solve(assumptions); 
  } 

  uint32_t AntomBase::Solve(const std::vector<uint32_t>& assumptions) 
  {
	++_satSolverCalls;

	/*
	std::cout << __func__ << " with assumptions: ";
	
	for( uint32_t i = 0; i != assumptions.size(); ++i )
	  {
		std::cout << helper::Lit(assumptions[i]) << " ";
	  }
	std::cout << std::endl;
	*/
	
    // Initialization.
    for (uint32_t t = 0; t < _antomSetting->threads; ++t)
      { _result[t] = ANTOM_UNKNOWN; }
    _sID = 0; 

    // (Re-)initialize the "Control" object.
    _control->ResetDone();

#ifdef PARALLEL
#pragma omp parallel 
    {
      // Every thread has its own unique ID.
      uint32_t id(omp_get_thread_num());
	  
      // Solve the CNF formula.
      _result[id] = _core[id]->Solve(assumptions); 
    }

    // Have we solved the CNF formula?
    uint32_t t(0); 
    for (t = 0; t < _antomSetting->threads; ++t)
      { 
		if (_result[t] != ANTOM_UNKNOWN)
		  { _sID = t; return _result[t]; }
      }
#else
	// Solve the CNF formula.
	_result[0] = _core[0]->Solve(assumptions); 

	if (_result[0] != ANTOM_UNKNOWN)
	  { _sID = 0; return _result[0]; }
#endif

    // If we reach this point, we've got a problem.
	// Sven: not in "basic operation timeout mode" :)
    //assert(t < _antomSetting->threads); 

    // Return UNKNOWN.
    return ANTOM_UNKNOWN;
  }

  uint32_t AntomBase::GetExtendedResult(void) const 
  { return _control->GetExtendedResult(); }

  // Stores the current status of all SAT solving cores.
  void AntomBase::SaveStatus(void)
  {
#ifdef PARALLEL
#pragma omp parallel for 
    for (uint32_t t = 0; t < _antomSetting->threads; ++t)
      { _core[t]->SaveStatus(); }
#else
	_core[0]->SaveStatus();
#endif
  }

  // Restores the status of all SAT solving cores saved before by "saveStatus()".
  void AntomBase::RestoreStatus(void)
  {
#ifdef PARALLEL
#pragma omp parallel for 
    for (uint32_t t = 0; t < _antomSetting->threads; ++t)
      { _core[t]->RestoreStatus(); }
#else
	_core[0]->RestoreStatus();
#endif
  }

  // Deletes the status of all SAT solving cores saved before by "saveStatus()". 
  void AntomBase::DeleteStatus(void)
  {
#ifdef PARALLEL
#pragma omp parallel for 
    for (uint32_t t = 0; t < _antomSetting->threads; ++t)
      { _core[t]->DeleteStatus(); }
#else
	_core[0]->DeleteStatus();
#endif
  } 

  // Resets antom and all SAT solving cores
  void AntomBase::Reset(void)
  {
#ifdef PARALLEL
#pragma omp parallel for 
    for (uint32_t t = 0; t < _antomSetting->threads; ++t)
      {
		_core[t]->Reset();
		_preprocessor[t]->Reset();
		_result[t] = ANTOM_UNKNOWN;
	  }
#else
	_core[0]->Reset();
	_preprocessor[0]->Reset();
	_result[0] = ANTOM_UNKNOWN;
#endif

	_control->Reset();

	DataReset();
  }

  void AntomBase::InstanceReset(void)
  {
#ifdef PARALLEL
#pragma omp parallel for 
    for (uint32_t t = 0; t < _antomSetting->threads; ++t)
      { 
		_core[t]->InstanceReset();
		_preprocessor[t]->InstanceReset();
		_result[t] = ANTOM_UNKNOWN;
	  }
#else
	_preprocessor[0]->InstanceReset();
	_core[0]->InstanceReset(); 
	_result[0] = ANTOM_UNKNOWN;
#endif

	_control->InstanceReset();

	DataReset();
  }

  void AntomBase::DataReset(void)
  {
	_globalPropertyTrigger.resize(1,0);
	_stacksize = 0;
	_satSolverCalls = 0;

	_lastModel.clear();

#ifndef NDEBUG
	_debugNetwork.clear();
#endif
  }

  // Writes current clauses into "db"
  void AntomBase::GetClauseDatabase(std::vector< std::vector< uint32_t > >& db) 
  { _core[_sID]->GetClauseDatabase( db ); }

  // Clears all datastructures by deleting every related clause, etc. 
  // with variable indices between "begin" and "end"
  void AntomBase::ClearVariables(uint32_t begin, uint32_t end)
  { 
	for (uint32_t t = 0; t < _antomSetting->threads; ++t)
	  { _core[t]->ClearVariables(begin, end); }
  }

  std::vector< std::pair<std::vector< uint32_t >, uint32_t > > AntomBase::GetConflictClauses(void) const
  {
	return _core[_sID]->GetConflictClauses();
  }

  /* antom preprocessor interface */ 

  /* Interface preprocessor */

  // De-/activates UPLA in preprocessing
  void AntomBase::SetUPLA(bool val)
  {
	_antomSetting->doUpla = val;
  }

  // De-/activates full subsumption check in preprocessing
  void AntomBase::SetSubsumption(bool val)
  {
	_antomSetting->doSubsumption = val;
  }

  // De-/activates variable elimination in preprocessing
  void AntomBase::SetVarElim(bool val)
  {
	_antomSetting->doVarElimination = val;
  }

  // De-/activates blocked clause elimination in preprocessing
  void AntomBase::SetBCE(bool val)
  {
	_antomSetting->doBce = val;
  }

  // De-/activates hidden tautology elimination in preprocessing
  void AntomBase::SetHTE(bool val)
  {
	_antomSetting->doHte = val;
  }

  // De-/activates hidden tautology elimination in preprocessing
  void AntomBase::SetHSE(bool val)
  {
	_antomSetting->doHse = val;
  }

  // De-/activates bounded variable addition in preprocessing
  void AntomBase::SetBVA(bool val)
  {
	_antomSetting->doBva = val;
  }

  // De-/activates bounded variable addition in preprocessing
  void AntomBase::SetTwoLiteralDiffBVA(bool val)
  {
	_antomSetting->bvaTwoLitDiff = val;
  }

  // De-/activates vivification in preprocessing
  void AntomBase::SetVivification(bool val)
  {
	_antomSetting->doVivification = val;
  }

  // Set thresholds for cost increase (default: 0 )
  // Variable elimination is only performed if formula is maximum increased by "increase"
  void AntomBase::SetVarIncrease(int32_t increase)
  {
	_antomSetting->varIncrease = increase;
  }

  // Sets maximum number of preprocessing loops
  void AntomBase::SetMaxLoops(uint32_t val)
  {
	_antomSetting->maxLoops = val;
  }

  // Should variable "var" be a "Don't Touch" variable?
  void AntomBase::SetDontTouch(uint32_t var, bool dt)
  {
	for (uint32_t t = 0; t < _antomSetting->threads; ++t) 
	  { _preprocessor[t]->SetDontTouch(var, dt); } 
  }

  // Return true if "var" is a "Don't Touch" variable
  bool AntomBase::IsDontTouch(uint32_t var)
  {
	return _preprocessor[0]->IsDontTouch(var); 
  }

  // De-/active preprocessing
  // 0 : no prepro
  // 1 : prepro
  // 2 : incremental prepro
  void AntomBase::SetPreprocessing(PreproType val)
  { 
	_antomSetting->doPreprocessing = val;
  }

  // De-/active inprocessing during solving
  void AntomBase::SetInprocessing(bool val)
  {
	_antomSetting->doInprocessing = val;
  }

  // Be verbose
  void AntomBase::SetVerbosity(uint32_t val)
  {
	_antomSetting->verbosity = val;
  }

  // Simplifies the current CNF formula by performing some preprocessing steps.
  // Returns FALSE if the formula is unsatisfiable, otherwise TRUE
  bool AntomBase::Preprocess(PreproType type)
  {
	_antomSetting->doPreprocessing = type;
	// Initialization.
	for (uint32_t t = 0; t < _antomSetting->threads; ++t)
	  { 
		_result[t] = ANTOM_UNKNOWN; 
	  }
	_sID = 0; 

	// (Re-)initialize the "Control" object.
	_control->ResetDone();

#ifdef PARALLEL
#pragma omp parallel 
	{
	  // Every thread has its own unique ID.
	  uint32_t id(omp_get_thread_num());

	  // Solve the CNF formula.
	  _result[id] = _preprocessor[id]->Preprocess(type); 

    // Have we solved the CNF formula?
    uint32_t t(0); 
    for (t = 0; t < _antomSetting->threads; ++t)
      {
		if (_result[t] != ANTO_UNKNOWN)
		  { _sID = t; return _result[t]; } 
	  }
	}
#else
	_result[0] = _preprocessor[0]->Preprocess(type); 

	if (_result[0] != ANTOM_UNKNOWN)
	  { _sID = 0; return _result[0]; } 
#endif

    // Return UNKNOWN.
    return ANTOM_UNKNOWN;
  }

  void AntomBase::EncodeOR(uint32_t output, uint32_t input1, uint32_t input2)
  {

#ifndef NDEBUG
	bool result;
	_debugNetwork.push_back( Gate(output,input1,input2, ORGATE) );
#endif

	/*
	std::cout << "\""<< output << "\" [color=green];" << std::endl;
	std::cout << "\"" << input1 << "\" -> \"" << out << "\";" << std::endl;
	std::cout << "\"" << input2 << "\" -> \"" << out << "\";" << std::endl;
	*/

	//std::cout << __func__ << output << " = " << input1 << " OR " << input2 << std::endl;

	std::vector<uint32_t> clause;
	clause.reserve(3);

	clause.push_back(output << 1); 
	clause.push_back( (input1<<1)^1); 

	if( _antomSetting->incrementalMode > 0 )
	  { clause.push_back( CurrentGlobalAssumption() ); }

#ifndef NDEBUG
    result = AddClause(clause); assert(result); 
#else
	AddClause(clause);
#endif

    clause.clear(); 
    clause.push_back(output << 1); 
	clause.push_back( (input2<<1)^1); 

	if( _antomSetting->incrementalMode > 0 )
	  { clause.push_back( CurrentGlobalAssumption() ); }

#ifndef NDEBUG
    result = AddClause(clause); assert(result); 
#else
	AddClause(clause);
#endif

	clause.clear(); 
	clause.push_back((output << 1) ^ 1); 
	clause.push_back( (input1<<1) ); 
	clause.push_back( (input2<<1) ); 
	
	if( _antomSetting->incrementalMode > 0 )
	  { clause.push_back( CurrentGlobalAssumption() ); }
		
#ifndef NDEBUG
	result = AddClause(clause); assert(result); 
#else
	AddClause(clause);
#endif

  }

  void AntomBase::EncodeAND(uint32_t output, uint32_t input1, uint32_t input2)
  {
#ifndef NDEBUG
	bool result;
	_debugNetwork.push_back( Gate(output,input1,input2, ANDGATE) );
#endif

	/*	
	std::cout << "\""<< output << "\" [color=red];" << std::endl;
	std::cout << "\"" << input1 << "\" -> \"" << out << "\";" << std::endl;
	std::cout << "\"" << input2 << "\" -> \"" << out << "\";" << std::endl;
	*/

	//std::cout << __func__ << output << " = " << input1 << " AND " << input2 << std::endl;

	std::vector<uint32_t> clause;
	clause.reserve(3);

	clause.push_back((output << 1) ^ 1); 
	clause.push_back( (input1<<1) ); 

	if( _antomSetting->incrementalMode > 0 )
	  { clause.push_back( CurrentGlobalAssumption() ); }

#ifndef NDEBUG
	result = AddClause(clause); assert(result); 
#else
	AddClause(clause);
#endif

	clause.clear(); 
	clause.push_back((output << 1) ^ 1); 
	clause.push_back( (input2<<1) ); 

	if( _antomSetting->incrementalMode > 0 )
	  { clause.push_back( CurrentGlobalAssumption() ); }
	
#ifndef NDEBUG
	result = AddClause(clause); assert(result); 
#else
	AddClause(clause);
#endif
	clause.clear(); 
		
    clause.push_back(output << 1); 
	clause.push_back( (input1<<1)^1); 
    clause.push_back( (input2<<1)^1); 

	if( _antomSetting->incrementalMode > 0 )
	  { clause.push_back( CurrentGlobalAssumption() ); }

#ifndef NDEBUG
    result = AddClause(clause); assert(result); 
#else
	AddClause(clause);
#endif
  }

  uint32_t AntomBase::CurrentGlobalAssumption(void) const
  {
	return (_globalPropertyTrigger[_stacksize]<<1);
  }


}
