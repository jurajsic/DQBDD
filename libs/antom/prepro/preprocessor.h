
/********************************************************************************************
preprocessor.h -- Copyright (c) 2013-2016, Sven Reimer

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

#ifndef ANTOM_PREPROCESSOR_H_
#define ANTOM_PREPROCESSOR_H_

#include <algorithm>
#include <iostream>
#include <cassert>
#include <vector>

#include "antom.h"
#include "helper.h"
#include "clause.h"
#include "control.h"
#include "statistics.h"
#include "watcher.h"
#include "reason.h"
#include "varcandidate.h"
#include "varheap.h"

/* For all routines:
   Add an extra loop, if something changes after unitPropagation, consider the subroutine again
   
   So far:
   At the end of every subroutine unitPropagation is performend, the subroutine does not take the result of the propagation into account
 */

namespace antom {

  class BCE;
  class BVA;
  class HTE;
  class UPLA;
  class ModelRebuilder;

  enum PreproFlag
	{
	  PF_NOTHING = 0,
	  PF_SUBSUMPTION = 1,
	  PF_UPLA = 1<<1,
	  PF_VAR_ELIM = 1<<2,
	  PF_BVA = 1<<3,
	  PF_BCE = 1<<4,
	  PF_HTE = 1<<5,
	  PF_VIVIFY = 1<<6,
	  PF_ALL = (1<<7)-1
	};

  PreproFlag inline operator|(PreproFlag a, PreproFlag b)
	{
	  return PreproFlag(static_cast<int>(a) | static_cast<int>(b));
	}
  PreproFlag inline operator&(PreproFlag a, PreproFlag b)
	{
	  return PreproFlag(static_cast<int>(a) & static_cast<int>(b));
	}
  PreproFlag inline operator~(PreproFlag a)
	{
	  return PreproFlag( ~static_cast<int>(a) & static_cast<int>(PF_ALL));
	}

  class Preprocessor
  {

  public:

	//Constructor.
	explicit Preprocessor ( Core* core );

	//Destructor.
	~Preprocessor ( void );

	// De-/activates UPLA in preprocessing
	void SetUPLA(bool val) 
	{ _setting->doUpla = val; }

	// De-/activates full subsumption check in preprocessing
	void SetSubsumption(bool val)
	{ _setting->doSubsumption = val; }

	// De-/activates variable elimination in preprocessing
	void SetVarElim(bool val) 
	{ _setting->doVarElimination = val; }

	// De-/activates blocked clause elimination in preprocessing
	void SetBCE(bool val)
	{ _setting->doBce = val; }

	// De-/activates hidden tautology elimination in preprocessing
	void SetHTE(bool val)
	{ _setting->doHte = val; }

	// De-/activates hidden subsumption elimination in preprocessing
	void SetHSE(bool val) 
	{ _setting->doHse = val; }

	// De-/activates bounded variable addition in preprocessing
	void SetBVA(bool val)
	{ _setting->doBva = val; }

	// De-/activates two literal difference in bounded variable addition
	void SetTwoLiteralDiffBVA(bool val)
	{ _setting->bvaTwoLitDiff = val; }

	// De-/activates vivification in preprocessing
	void SetVivification(bool val)
	{ _setting->doVivification = val; }

	// Set thresholds for cost increase (default: 0 )
	// Variable elimination is only performed if formula is maximum increased by "increase"
	void SetVarIncrease(int32_t increase)
	{ _setting->varIncrease = increase; }

	// Sets maximum number of preprocessing loops
	void SetMaxLoops(uint32_t val)
	{ _setting->maxLoops = val; }

	// Returns some preprocessor related statistics
	uint32_t BinaryConstants(void) const 
	{ return _statistics.constantVariables; }
	uint32_t BinaryEquivalences(void) const
	{ return _statistics.equivalentVariables; }
	uint32_t UplaConstants(void) const
	{ return _statistics.uplaConstantVariables; }
	uint32_t UplaEquivalences(void) const
	{ return _statistics.uplaEquivalentVariables; }
	uint32_t SatConstants(void) const 
	{ return _statistics.constantVariablesBySAT; }
	uint32_t VariableEliminations(void) const 
	{ return _statistics.resolvedVariables; }
	uint32_t LiteralEliminations(void) const 
	{ return _statistics.resolvedLiterals; }
	uint32_t BlockedClauses(void) const 
	{ return _statistics.blockedClauses; }
	uint32_t HiddenTautologies(void) const 
	{ return _statistics.hiddenTautologies; }
	uint32_t HiddenSubsumptions(void) const 
	{ return _statistics.hiddenSubsumptions; }
	uint32_t MonotoneVariables (void) const 
	{ return _statistics.monotoneVariables; }
	uint32_t DcVariables(void) const 
	{ return _statistics.dontcareVariables; }
	uint32_t SubsumedClauses(void) const 
	{ return _statistics.subsumptions; }
	uint32_t SelfsubsumedLiterals(void) const 
	{ return _statistics.selfSubsumptions; }
	uint32_t BvaVariables(void) const 
	{ return _statistics.bvaVariables; }
	uint32_t BvaLiterals(void) const 
	{ return _statistics.bvaLiterals; }
	uint32_t VivifySubsumptions(void) const 
	{ return _statistics.vivifySubsumptions; }
	uint32_t VivifyUnits(void) const 
	{ return _statistics.vivifyUnits; }
	uint32_t VivifyDiff(void) const 
	{ return _statistics.vivifyDiff; }
	uint32_t UnitPropagations(void) const 
	{ return _statistics.unitPropagations; }
	double RuntimeUPLA (void) const 
	{ return _statistics.runtime_upla; }
	double RuntimeSubsumption(void) const 
	{ return _statistics.runtime_subsumption; }
	double RuntimeVarElim(void) const 
	{ return _statistics.runtime_varElim; }
	double RuntimeBCE(void) const 
	{ return _statistics.runtime_bce; }
	double RuntimeBVA(void) const 
	{ return _statistics.runtime_bva; }
	double RuntimeHTE(void) const 
	{ return _statistics.runtime_hte; }
	double RuntimeVivify(void) const 
	{ return _statistics.runtime_vivify; }
	//	double runtimeSATConst ( void ) const;
	double RuntimePrepro(void) const 
	{ return _statistics.runtime_preprocessing; }

	// Should variable "var" be a "Don't Touch" variable?
	void SetDontTouch(uint32_t var, bool dt = true) 
	{
	  assert( var <= _variables ); 

	  // Update Preprocessor Data structure if don't touch variables are declared
	  if( _donttouch.size() <= _variables )
		{ UpdateDataStructures(); }

	  _donttouch[var] = dt; 
	}

	// Returns true, if "var" is a "Don't Toch" variable
	bool IsDontTouch(uint32_t var)
	{ 
	  assert ( var <= _variables );
	  if( _donttouch.size() < _variables )
		{ return false; }
	  return _donttouch[var];
	}

	bool IsUsed(uint32_t var)
	{
	  return ( !_deleted[var] && !_assignment[var<<1] && !_assignment[(var<<1)^1] );
	}

	// Simplifies the current CNF formula by performing some preprocessing steps.
	// Returns ANTOM_UNSAT if the formula is unsatisfiable, otherwise ANTOM_UNKNOWN
	// type = "PREPROCESS" : standard preprocessing
	// type = "INPROCESS" : inprocessing (lightweighted preprocessing)
	// type = "INCREMENTAL" : incremental preprocessing
	// With e_incremental the preprocessor only performs preprocessing on the lastly added variables and clauses
	uint32_t Preprocess(PreproType type);

	// Some fast preprocessing routines, which can be called often
	// If incremental = TRUE: Just check variables with index >= "_firstPreIndex"
	bool FastPreprocess(bool& didsomething, bool incremental = false);

	void ExtendModel(void);

	void ClearRestoreData(uint32_t begin, uint32_t end);

	// Set variable index at which point the next preprocessing starts
	void SetPreVarIndex(uint32_t val)
	{ _firstPreVarIndex = val; }

	// Set index for internal clause datastructure at which point the next preprocessing starts
	// Assumes that the clause database is not shrinked until next incremental preprocessor call
	void SetPreClauseIndex(void)
	{ _firstPreClauseIndex = static_cast<uint32_t>(_clauseDatabase.size())+1; }

	// DEBUGGING
	void PrintDatabase(bool printOccur = true, bool printWatches = false) const;

	// Resets the Preprocessore. The following status flags/variables remain untouched:
	// * The references of the core
	// * Information about replaced variables and clauses 
	// * Don't touch variables
	// * Overall statistics 
	// * DoPreprocessing and DoInprocessing flags
	void Reset(void);

	void InstanceReset(void);

  protected:

	struct OccurenceSorter
	{
	  explicit OccurenceSorter( std::vector< size_t >& preproOccurSizes ):
		occurSizes(preproOccurSizes)
	  {}
	  bool operator() (uint32_t l1, uint32_t l2) const;
	  std::vector< size_t > occurSizes;
	};

	/* Begin: some helper functions */

    // Updates all data structures depending on the number of variables to be able to handle "_variables" variables.

	void UpdateDataStructures(void);

	// Adds (lit1 + lit2) to binaries. 
	// If (lit1 + lit2) already exists, do nothing#
	// Returns false, if binary already exists
	bool AddBinary(uint32_t lit1, uint32_t lit2, bool learned);

	// Remove the "pos"th binary entry in the binary list of "lit"
	// "pos" and "size" will be updated
	void RemoveBinary(uint32_t lit, uint32_t& pos, uint32_t& size);

	// Remove the binary (lit1 + lit2) from binary list of "lit1"
	void RemoveBinary(uint32_t lit1, uint32_t lit2);

	// Returns "true" if the binary (lit1,lit2) exists
	bool HasBinary(uint32_t lit1, uint32_t lit2) const;

	// Removes the clause "c" from watchlist of "lit"
	void RemoveWatch( uint32_t lit, CRef c);

	// Removes "clause" from occurence list of "lit"
	void RemoveFromOccurenceList(uint32_t lit, CRef cr);

	// Removes every occurence of "clause"
	// Mark "clause" as "to deleted"
	void ClearAllOccurences(CRef cr);

	void ClearClauseDatabase(void);

	void CountOccurences(bool countLearned);
	void CountBinaryOccurences(void);

	// Returns number occurences of "lit" (including binaries)
	uint32_t OccurenceCount(uint32_t lit, bool countLearned = true) const;

	uint32_t GetClauseSize(const Watcher& watcher) const;

	// Removes "literal" from clause
	// Eventually introduces new binary and mark old n-nary for deletion
	// Returns "true" if n-nary clause is deleted, "false" otherwisse
	bool StrengthenClause(CRef cr, uint32_t literal, bool keepoccurence = false);

	// 1. Extract all binaries of watchlist and push them in extra data structure
	// (this can also be done during "addClause", if preprocessing is enabled)
	// 2. Creates occurence lists
	// Assumes that there are no duplicated binaries in "_watches"
	void PreparePreprocessing(PreproType type);

	// Update all watchlists
	void UpdateWatches(bool showstats = true);

	void UpdateOccurenceLists(ClauseAllocator& to);

	bool CopyBinaryList(uint32_t toReplace, uint32_t replace); 
	void CopyOccurenceList(uint32_t toReplace, uint32_t replace); 

	void PreserveModel(uint32_t resvar);
	// Replace every occurence of "toReplace" with "replace"
	// Refresh occurence lists
	// Preserves model for replaced variable
	bool ReplaceVariable(uint32_t toReplace, uint32_t replace);

	// Delete all occurences and clauses of "var" in database
	void DeleteVariable(uint32_t var);

	// Merges two n-nary clauses with common literal "reslit"
	// Store new clause in "newClause"
	void MergeClauses(Clause& clause1, Clause& clause2, uint32_t reslit, std::vector< uint32_t >& newClause) const;

	// Merges a n-nary clauses "c" with a binary (reslit + otherlit) and common literal "reslit"
	// Store new clause in "newClause"
	void MergeClauses(Clause& clause, uint32_t reslit, uint32_t otherlit, std::vector< uint32_t >& newClause) const;

	// Counts the literals of merge of clause "c1" and "c2"
	// Return 0, if result is a tautology
	uint32_t CountMergeClauses(const Clause& clause1, const Clause& clause2, uint32_t reslit) const;

	// Counts the literals of merge of clause "c" and "(reslit + otherlit)"
	// Return 0, if result is a tautology
	uint32_t CountMergeClauses(const Clause& clause, uint32_t reslit, uint32_t otherlit) const;

	// Estimates costs for variable elimination of "var"
	int32_t EstimateCosts(uint32_t var);

	// Estimate costs for variable elimination of "var" with more accurate "countMergeClauses"
	int32_t EstimateCosts2(uint32_t var) const;

	// Estimate costs for variable elimination of "var" with more accurate "countMergeClauses"
	void EstimateCosts3(VarCandidate& varcand) const;

	// Adds a clause to database
	// Do _not_ update watch lists
	// A binary clause is put into the binary datastructure of preprocessor
	// Updates occurence lists for n-nary clauses
	// This method is only used in the preprocessor 
	bool AddClausePrepro(std::vector<uint32_t>& clause, bool updatewatches = false, uint32_t lbd = 1);

	/* End: some helper functions */

	// Propagate new units
	// Delete all satisfied clauses and their occurence	
	// Strengthen Clauses
	// Returns "false" if a contradiction occurs (formula is UNSAT), otherwise true
	bool PropagateUnits(void);

	// Detect trivial monotone/pure literals
	// If incremental = TRUE: Just check variables with index >= "_firstPreIndex"
	void DetectMonotone(bool incremental = false);

	// Search in binary clauses for two cases:
	// 1. [(a + b) * (a+ ~b)] => a is constant, imply a
	// 2. [(a + ~b) * (~a + b )] => a and b are equivalent, replace all occurences of a with b (or b with a)
	// If incremental = TRUE: Just check variables with index >= "_firstPreIndex"
	bool FindBinaryConsAndEquiv(bool& didsomething);

	// Performs full subsumptioncheck over all clauses
	// If incremental = TRUE: Just check variables with index >= "_firstPreIndex"
	bool FullSubsumptionCheck(bool incremental = false);

	template<class T> 
	uint32_t GetSubsumerCandidate(const T& c);

	template<class T>
	  //bool Subsumes(const Clause& c1, const T& c2, uint64_t c2sign);
	  bool Subsumes(const Clause& c1, const T& c2, uint32_t c2sign);

	template<class T>
	  //bool IsSubsumed(const T& clause, uint64_t sign, bool islearned);
	  bool IsSubsumed(const T& clause, uint32_t sign, bool islearned);

	template<class T>
	  //bool IsSubsumedExcept(const T& clause, uint64_t sign, bool islearned, CRef except);
	  bool IsSubsumedExcept(const T& clause, uint32_t sign, bool islearned, CRef except);

	// Mark all n-nary clauses which are subsumed by the binary (i1 + i2)
	// Perform self-subsumption, where one clause is a binary
	void CheckBinSub(uint32_t i1, Watcher& binaryclause);

	//bool CheckAllBinaries(const std::vector< uint32_t >& clause, uint64_t sign);
	bool CheckAllBinaries(const std::vector< uint32_t >& clause, uint32_t sign);

	//bool CheckAllBinariesExcept(const std::vector< uint32_t >& clause, uint64_t sign, uint32_t lit, uint32_t lit2);
	bool CheckAllBinariesExcept(const std::vector< uint32_t >& clause, uint32_t sign, uint32_t lit, uint32_t lit2);

	// Is "clause" subsumed by some clause in database?
	// Assumes that "clause" is sorted
	//bool IsForwardSubsumed(const std::vector< uint32_t >& clause, uint64_t sign);
	bool IsForwardSubsumed(const std::vector< uint32_t >& clause, uint32_t sign);

	// Is "clause" subsumed by some clause in database (Except "except")?
	// Assumes that "clause" is sorted
	//bool IsForwardSubsumedExcept(const std::vector< uint32_t >& clause, uint64_t sign, CRef except);
	bool IsForwardSubsumedExcept(const std::vector< uint32_t >& clause, uint32_t sign, CRef except);

	// Is "clause" subsumed by some clause in database (Except "except")?
	// Assumes that "clause" is sorted
	//bool IsForwardSubsumedExcept(const std::vector< uint32_t >& clause, uint64_t sign, uint32_t lit1, uint32_t lit2);
	bool IsForwardSubsumedExcept(const std::vector< uint32_t >& clause, uint32_t sign, uint32_t lit1, uint32_t lit2);

	// Performs variable elimination
	// If incremental = TRUE: Just check variables with index >= "_firstPreIndex"
	bool VarElimination(bool incremental = false);

	// Do Vivification, see "Vivifying Propositional Clausal Formulae", Piette, Hamadi, Sais
	bool Vivify(void);

	void SetFlag(PreproFlag flag);
	void UnsetFlag(PreproFlag flag);
	bool GetFlag(PreproFlag flag) const;

	// DEBUGGING functions
	void PrintOccurenceList(uint32_t lit) const;
	void PrintWatchList(uint32_t lit) const;
	void PrintCompleteLists(uint32_t var) const;
	void CheckWatchLists(void) const;
	int32_t IsAssigned(uint32_t lit) const;

	friend class Core;
	friend class BVA;
	friend class BCE;
	friend class HTE;
	friend class UPLA;
	friend class ModelRebuilder;

	/* Core references, see also core.hpp for description */
	Core* _core;
	Control* _control;
	Settings* _setting;
	// The bounded variable class
	BCE* _bce;
	BVA* _bva;
	HTE* _hte;
	UPLA* _upla;
	ModelRebuilder* _rebuilder;

	// defines the limit for each method
	int64_t _methodLimit;

	bool& _emptyClause;
	uint32_t& _variables;
	uint32_t& _dsImplIndex;
	uint32_t& _dsEndIndex;
	uint32_t& _decisionLevel;
	std::vector<unsigned char>& _assignment;
	std::vector<uint32_t>& _level;
    std::vector<uint32_t>& _decisionStack;
	std::vector<uint32_t>& _model;	
	std::vector< Reason >& _forcing;
	std::vector< unsigned char >& _deleted;
	std::vector< CRef >& _clauseDatabase;
	std::vector< std::vector< Watcher > >& _binaries;
	ClauseAllocator& _ca;

	/* Additional preprocessing data structures */
	std::vector< std::vector< CRef > > _occur;
	std::vector< size_t > _occurCounts;
	std::vector< size_t > _occurenceCount;
	std::vector< unsigned char > _donttouch;
	std::vector< VarCandidate > _varCandidates;
	
	VarHeap< helper::DescendingOrder, size_t > _binCandidates;

	uint32_t _firstPreVarIndex;
	uint32_t _firstPreClauseIndex;

	/* Some statistics and control values*/
	uint32_t _lastImplIndex;

	Statistics& _statistics;

	PreproFlag _pFlag;
	bool _okay;

  private :
	// Copy constructor.
    Preprocessor (const Preprocessor&) = default;

    // Assignment operator.
    Preprocessor& operator = (const Preprocessor&) = default;
  };
}

#endif
