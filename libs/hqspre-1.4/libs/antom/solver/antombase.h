
/********************************************************************************************
antombase.h -- Copyright (c) 2013-2017, Tobias Schubert, Sven Reimer, Tobias Paxian

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

#ifndef ANTOM_ANTOMBASE_H_
#define ANTOM_ANTOMBASE_H_

// Include standard headers.
#include <iostream>
#include <cassert>
#include <vector>

//#define PARALLEL

#ifdef PARALLEL
#include <omp.h>
#endif

#include "settings.h"

namespace antom
{

  // Some definitions.
#define ANTOM_UNKNOWN                  0
#define ANTOM_SAT                     10
#define ANTOM_UNSAT                   20
#define ANTOM_UNSAT_WITH_ASSUMPTION   30
#define ANTOM_UNKNOWN_WITH_PRE_RESULT 40

  // Some forward declarations.
  class Core;
  class Control;
  class Preprocessor;
  class SolverProxy;
  class Settings;
  //class Circuit;

  struct Gate 
  {
  Gate(uint32_t out, uint32_t i1, uint32_t i2, GateType t) :
	output(out),
	  inputs({i1, i2}),
	  type(t)
	{}

  Gate(uint32_t out, uint32_t i, GateType t) :
	output(out),
	  inputs({i, 0}),
	  type(t)
	{}

  Gate(uint32_t out, const std::vector<uint32_t>& ins, GateType t) :
	output(out),
	  inputs(ins),
	  type(t)
	{}

	~Gate() = default;
	
	uint32_t output;
	std::vector<uint32_t> inputs;
	GateType type;
  };

  // The "Antom" class.
  class AntomBase
  {

  public:

    // Constructor.
    explicit AntomBase();

    // Destructor.
    ~AntomBase(void);

	void SetThreads(uint32_t threads);

    /* Some statistics */

    // Returns the ID of the thread that was able to solve the CNF.
    uint32_t SolvingThread(void) const; 

    // Returns the number of variables for which memory has been reserved.
    uint32_t Variables(void) const;

    // Returns the current number of clauses within the clause database. 
    uint32_t Clauses(void) const; 

    // Returns the current number of literals within the clause database. 
    uint32_t Literals(void) const;

    // Returns the number of decisions made so far.
    uint32_t Decisions(void) const;

    // Returns the number of BCP operations performed so far.
    uint32_t Bcps(void) const;

    // Returns the number of implications found so far. 
    uint64_t Implications(void) const;

    // Returns the number of conflicts encountered so far.
    uint32_t Conflicts(void) const;

    // Returns the number of restarts performed so far.
    uint32_t Restarts(void) const;
	uint32_t BlockedRestarts(void) const;

    // Returns the number of database simplifications performed so far.
    uint32_t Simplifications(void) const;

    // Returns the number of binary clauses deduced due to "Lazy Hyper Binary Resolution".
    uint32_t Lhbr(void) const;

    uint32_t UsedVariables(void) const;
    uint32_t CurrentBinaryClauses(void) const;
    uint32_t CurrentTernaryClauses(void) const;
    uint32_t CurrentNaryClauses(void) const;

    // Returns the number of unit clauses deduced due to conflict analysis.
    uint32_t LearntUnitClauses(void) const;

    // Returns the number of binary clauses deduced due to conflict analysis.
    uint32_t LearntBinaryClauses(void) const;

    // Returns the number of ternary clauses deduced due to conflict analysis.
    uint32_t LearntTernaryClauses(void) const;

	uint32_t MinimizedLiterals(void) const;

    // Returns the number of synchronizations performed so far. 
    uint32_t Synchronizations(void) const;

    // Returns the number of inprocessings steps during solving main routine
    uint32_t Inprocessings(void) const;

    // Returns the average "Literals Blocks Distance" of all conflict clauses deduced so far.
    double AvgLBD(void) const;

    // Returns the average length of all conflict clauses deduced so far.
    double AvgCCLength(void) const;

    // Returns the solver's average decision level before backtracking.
    double AvgDL(void) const;

    // Returns the average number of decision levels cleared during conflict analysis. 
    double AvgDLclearedCA(void) const;

    // Returns the average number of variables getting unassigned during conflict analysis. 
    double AvgVarsUnassignedCA(void) const;

	uint32_t BinaryConstants(void) const;
	uint32_t BinaryEquivalences(void) const;
	uint32_t UplaConstants(void) const;
	uint32_t UplaEquivalences(void) const;
	uint32_t SatConstants(void) const;
	uint32_t VariableEliminations(void) const;
	uint32_t LiteralEliminations(void) const;
	uint32_t BlockedClauses(void) const;
	uint32_t HiddenTautologies(void) const;
	uint32_t HiddenSubsumptions(void) const;
	uint32_t MonotoneVariables(void) const;
	uint32_t DcVariables(void) const;
	uint32_t SubsumedClauses(void) const;
	uint32_t SelfsubsumedLiterals(void) const;
	uint32_t BvaVariables(void) const;
	uint32_t BvaLiterals(void) const;
	uint32_t VivifySubsumptions(void) const;
	uint32_t VivifyUnits(void) const;
	uint32_t VivifyDiff(void) const;
	uint32_t UnitPropagations(void) const;
	double RuntimeUPLA(void) const;
	double RuntimeSubsumption(void) const;
	double RuntimeVarElim(void) const;
	double RuntimeBCE(void) const;
	double RuntimeBVA(void) const;
	double RuntimeHTE(void) const;
	double RuntimeVivify(void) const;
	double RuntimeSATConst(void) const;
	double RuntimePrepro(void) const;

    // Returns a reference to either the satisfying variable assignment (after calling one of the 
    // "solve()" routines) or the set of currently assigned variables (after calling one of the 
    // two "deduceAssumptions()" routines). Example:
    // model[17] =  0 --> x17 = unassigned
    // model[17] = 35 --> x17 = false 
    // model[17] = 34 --> x17 = true
    // In case neither "solve()/maxSolve()" nor "deduceAssumptions()" has been called beforehand, the 
    // vector contains invalid data. 
    const std::vector<uint32_t>& Model(void) const;

	// Creates and returns a fresh variable index not used so far.
    uint32_t NewVariable(void); 

	/* solver options */

    // Activates or deactivates "Lazy Hyper Binary Resolution". 
    // val = TRUE  --> LHBR enabled (default).
    // val = FALSE --> LHBR disabled.
    void SetLHBR(bool val); 

    // Sets the unit factor of both restart strategies -- Luby & glucose-like -- to "val" (default: 8). 
    // The unit factor directly corresponds to the interval between two restart operations. 
    void SetLuby(uint32_t val);

    // Sets the decision strategy of group "group" to mode "val". 
    // Currently, there are four modes that differ wrt. the polarity of a decision variable:
    // 0 (default) --> Use the variable's cached polarity together with antom's "polarity toggling scheme". 
    // 1           --> Use the variable's cached polarity only.
    // 2           --> The polarity will be set to FALSE regardless of the cached value. 
    // 3           --> The polarity will be set to TRUE regardless of the cached value.
    // Furthermore, antom maintains variable orderings: As long the group with the lowest index is non-empty,
	// variables from that group will be preferred to serve as decision variables. By default, all 
    // variables belong to "group 0".
    void SetDecisionStrategy(uint32_t val, uint32_t group);

	// Like "setDecisionStrategy()" for a specific variable instead of a group
    void SetDecisionStrategyForVariable(uint32_t val, uint32_t var);

	// Set initial polarity of the variable "var" to "pol"
	void SetPolarity(uint32_t var, bool pol);

    // Sets the restart strategy to model "val":
    // 0 --> Luby (default).
    // 1 --> Glucose-like.       
    void SetRestartStrategy(RestartStrategy val);

	void SetSimplifyStrategy(SimplifyStrategy val);

	void SetSimplifyActivity(SimplifyActivity val);

    // Sets the decay factor to "val" (default: 1.05). The decay factor is responsible 
    // for controlling how the variable activities evolve during the search process. 
    void SetDecayFactor(double val);

    // Sets the maximum variable index to "max". 
    void SetMaxIndex(uint32_t max);

    // Sets the group of variable "var" to "grp". See "setDecisionStrategy()" for more details. 
    void SetVarGroup(uint32_t var, uint32_t grp);

	void UseTernaryClauses(bool val);

	// Returns whether variable was already deleted in pre-/in-processing
	bool IsDeleted(uint32_t var) const;

	// Defines CPU limit in seconds (default: no CPU limit)
	void SetCPULimit(double t);

	// Defines Memory limit in MB (default: no Memory limit)
	void SetMemoryLimit(uint32_t m);

    // Adds a clause to the clause databases of all threads. Returns FALSE if the CNF formula is unsatisfiable, 
    // otherwise TRUE will be returned. Assumes that the solver is on decision level 0 and that "clause" is not 
    // empty. Furthermore, all literals have to be encoded as follows, having variable indices greater 0:
    //  x3 <--> (3 << 1)     = 6
    // -x3 <--> (3 << 1) + 1 = 7
    // All clauses inserted into the clause database using "addClause()" are assumed to belong to 
    // the original CNF formula (having a "Literals Blocks Distance" of 1). 
    // IN THE MULTI-THREADED MODE, "maxSetIndex()" HAS TO BE CALLED BEFOREHAND.
    bool AddClause(std::vector<uint32_t>& clause, uint32_t lbd = 1);

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
    // 3) "setMaxIndex()" HAS BEEN CALLED BEFOREHAND.
    bool AddClause(uint32_t* lits, uint32_t num, uint32_t lbd = 1);

	// Add unit clause, using literal encoding as in "addClause()"
	bool AddUnit(uint32_t lit);

	// Gate interface of solver
	// In- and Outputs are always literals => supports natively inverted out- and inputs
	// Returns false if instance get unsatisfiable after insertion

	#if 0
	// Adds encoding of an AND-Gate (output = input1 * input2) to the solver
	bool AddAndGate(uint32_t input1, uint32_t input2, uint32_t output);
	bool AddAndGate(const std::vector<uint32_t>& inputs, uint32_t output);
	// Adds encoding of an NAND-Gate (output = !(input1 * input2)) to the solver
	bool AddNandGate(uint32_t input1, uint32_t input2, uint32_t output);
	bool AddNandGate(const std::vector<uint32_t>& inputs, uint32_t output);
	// Adds encoding of an OR-Gate (output = input1 + input2) to the solver
	bool AddOrGate(uint32_t input1, uint32_t input2, uint32_t output);
	bool AddOrGate(const std::vector<uint32_t>& inputs, uint32_t output);
	// Adds encoding of an NOR-Gate (output = !(input1 + input2)) to the solver
	bool AddNorGate(uint32_t input1, uint32_t input2, uint32_t output);
	bool AddNorGate(const std::vector<uint32_t>& inputs, uint32_t output);
	// Adds encoding of an XOR-Gate (output = input1 XOR input2) to the solver
	bool AddXorGate(uint32_t input1, uint32_t input2, uint32_t output);
	// Adds encoding of an XNOR-Gate (output = !(input1 XOR input2)) to the solver
	bool AddXnorGate(uint32_t input1, uint32_t input2, uint32_t output);
	#endif

	/* Solver interface */

    // Collects data for performing a regression analysis afterwards (--> SATzilla-like SAT solving). 
    // NOTE, THAT "getDataRegressionAnalysis()" CAN BE USED IN SINGLE-THREADED MODE ONLY.
    void GetDataRegressionAnalysis(void);

    // Solves the current CNF formula using the most promising configuration of antom. The return values are SAT/UNSAT.
    // NOTE, THAT "solveSATzilla()" CAN BE USED IN SINGLE-THREADED MODE ONLY.
    uint32_t SolveSATzilla(void);

    // Performs unit propagation, taking the current CNF and the specified assumptions into account. Returns 
    // FALSE if a conflict occurs, otherwise the return value is TRUE. NOTE, THAT "deduceAssumptions()" CAN BE 
    // USED IN SINGLE-THREADED MODE ONLY.
    bool DeduceAssumptions(void);
    bool DeduceAssumptions(const std::vector<uint32_t>& assumptions);

	// Returns the vector of failed assumptions leading to UNSAT result
	// Always returns the last failed vector
	// Vector is empty if assumptions are not directly responsible for UNSAT result
	std::vector<uint32_t > GetFailedAssumptions(void) const;

	// Simplify Database
	void Simplify(void);

    // Solves the current CNF formula, taking assumptions (if specified) into account. The assumptions have to be encoded in the 
    // same way as the literals of a clause (see "addClause()"). The return values are SAT/UNSAT. In the multi-threaded mode, the 
    // CNF formula gets solved by "_threads" threads running in parallel, according to the well-known algorithm portfolio scheme. 
    uint32_t Solve(void);
    uint32_t Solve(const std::vector<uint32_t>& assumptions);

	uint32_t GetExtendedResult(void) const;

	// Writes current clauses into "db"
	void GetClauseDatabase(std::vector< std::vector< uint32_t > >& db);

	// Returns learnt conflict clauses and their "LBD"-value
	// Useful for "Internal constraints replication"
	// See O. Strichman "Accelerating Bounded Model Checking of Safefy Properties" for more details
	std::vector< std::pair<std::vector< uint32_t >, uint32_t > > GetConflictClauses(void) const;

    // Stores the current status of all SAT solving cores.
    void SaveStatus(void);

    // Restores the status of all SAT solving cores saved before by "saveStatus()".
    void RestoreStatus(void);
    
    // Deletes the status of all SAT solving cores saved before by "saveStatus()". 
    void DeleteStatus(void); 

    // Resets all SAT solving cores.
    void Reset(void);

	void InstanceReset(void);

	/* Preprocessor interface */

	// De-/activates UPLA in preprocessing (default: TRUE)
	void SetUPLA(bool val);

	// De-/activates full subsumption check in preprocessing  (default: TRUE)
	void SetSubsumption(bool val);

	// De-/activates variable elimination in preprocessing  (default: TRUE)
	void SetVarElim(bool val);

	// De-/activates blocked clause elimination in preprocessing  (default: TRUE)
	void SetBCE(bool val);

	// De-/activates hidden tautology elimination in preprocessing  (default: TRUE)
	void SetHTE(bool val);

	// De-/activates hidden tautology elimination in preprocessing  (default: TRUE)
	void SetHSE(bool val);

	// De-/activates bounded variable addition in preprocessing  (default: TRUE)
	void SetBVA(bool val);

	// De-/activates two literal difference in bounded variable addition
	void SetTwoLiteralDiffBVA(bool val);

	// De-/activates vivification in preprocessing (default: TRUE)
	void SetVivification(bool val);

	// Set thresholds for cost increase factor (default: 0 ) in variable elimination
	void SetVarIncrease(int32_t increase);

	// Sets maximum number of preprocessing loops (default: 1.000.000)
	void SetMaxLoops(uint32_t val);

	// Should variable "var" be a "Don't Touch" variable?
	void SetDontTouch(uint32_t var, bool dt = true);

	// Return true if "var" is a "Don't Touch" variable
	bool IsDontTouch(uint32_t var);

	// De-/actives preprocessing (default: FALSE)
	void SetPreprocessing(PreproType val);

	// De-/actives inprocessing during solving (default: FALSE)
	void SetInprocessing(bool val);

	// De-/active inprocessing during max-antom (default: FALSE)
	void SetMaxInprocessing(bool val);

	// Simplifies the current CNF formula by performing some preprocessing steps.
	// Returns FALSE if the formula is unsatisfiable, otherwise TRUE
	// type = "PREPROCESS" : standard preprocessing
	// type = "INPROCESS" : inprocessing (lightweighted preprocessing)
	// type = "INCREMENTAL" : incremental preprocessing
	// With "e_incremental" the preprocessor only performs preprocessing on the lastly added variables and clauses
	bool Preprocess(PreproType type = PREPROCESS);

	// Be verbose
	void SetVerbosity(uint32_t val);

	// Dumps core data as cnf into std::cout
	void DumpCNF(bool printAssignment = false) const;

	// get trivial Assignemnts of solver core
	void TrivialAssignment(void) const;

	void EncodeOR(uint32_t output, uint32_t input1, uint32_t input2);
	void EncodeAND(uint32_t output, uint32_t input1, uint32_t input2);

  protected:

	// Copy constructor.
    AntomBase(const AntomBase&);

    // Assignment operator.
    AntomBase& operator=(const AntomBase&);

	void DataReset(void);
		
	// Clears all datastructures by deleting every related clause, etc. 
	// with variable indices between "begin" and "end"
	void ClearVariables(uint32_t begin, uint32_t end);

	uint32_t CurrentGlobalAssumption(void) const;

	// Adds a basic gate
	// These method are called by the addGate interfaces
	// Gatetypes by (neginput,negoutput):
	// (false, false) = NOR
	// (false, true) = OR
	// (true, false) = AND
	// (true, true) = NAND
	bool AddBaseGate(const std::vector<uint32_t>& inputs, bool neginput, uint32_t output, bool negoutput);
	bool AddBaseXorGate(uint32_t input1, uint32_t input2, uint32_t output);

	Settings* _antomSetting;

    // A pointer to the "Control" object.
    Control* _control;

    // The SAT solving cores.
    std::vector<Core*> _core; 

    // The SAT solving pre-/inprocessors.
    std::vector<Preprocessor*> _preprocessor;

	// The SAT solving back-end
	//SolverProxy* _solver;

	// The stored circuit
	//Circuit* _circuit;

	// The ID of the thread that was able to solve the CNF.
    uint32_t _sID; 

    // For each core we maintain an uint32_t to represent (temporary) results.
    std::vector<uint32_t> _result; 

	// For incremental maxsat: 

	// Variable index triggering all property clauses in incremental mode
	std::vector< uint32_t > _globalPropertyTrigger;

	// pointer to current valid global assumption
	uint32_t _stacksize;

	uint32_t _satSolverCalls;

	// Store last model 
	std::vector< uint32_t > _lastModel;

#ifndef NDEBUG
	std::vector< Gate > _debugNetwork;
#endif
  };
}

#endif
