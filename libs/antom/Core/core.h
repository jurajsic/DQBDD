
/********************************************************************************************
core.h -- Copyright (c) 2013-2017, Tobias Schubert, Sven Reimer

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

#ifndef ANTOM_CORE_H_
#define ANTOM_CORE_H_

// Include standard headers.
#include <sys/resource.h>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <vector>
#include <cstdint>

// Include antom related headers.
#include "clause.h"
#include "varheap.h"
#include "boundedQueue.h"
#include "control.h"
#include "preprocessor.h"

#include "watcher.h"
#include "reason.h"

#define CLEARWATCHES

namespace antom
{
  class SolverState;

  // The "Core" class.
  class Core
  {

	// A struct to be able to store the current state of the SAT solving core. 
	struct SolverState
	{
	  Statistics stats;
	  uint64_t basicOps; 
	  uint32_t variables; 
	  uint32_t decisionLevel;   
	  uint32_t dsEndIndex;
	  uint32_t dsImplIndex;
	  std::vector<DecisionStrategy> modeDS; 
	  std::vector<DecisionStrategy> modeDSForVariables; 
	  double incVarActivity;
	  bool emptyClause;
	  uint32_t startPtr;
	  uint32_t endPtr;
	  std::vector<bool> decVarSign;
	  std::vector<unsigned char> assignment;
	  std::vector<unsigned char> polarity;
	  std::vector<uint32_t> model;
	  std::vector<uint32_t> level; 
	  std::vector<uint32_t> decisionStack;
	  std::vector<uint32_t> dl2ds;   
	  std::vector<uint32_t> varGroup;    
	  std::vector<double> activity;  
	  std::vector<std::vector<Watcher> > watches;
	  std::vector<Reason> forcing;
	  ClauseAllocator ca;
	  std::vector<CRef> clauseDatabase;

	SolverState(void):
	  stats(),
	  basicOps(0), 
		variables(0),
		decisionLevel(0),
		dsEndIndex(0),
		dsImplIndex(0),
		modeDS(),
		modeDSForVariables(),
		incVarActivity(0.0),
		emptyClause(false),
		startPtr(1),
		endPtr(1),
		decVarSign(),
		assignment(0),
		polarity(0),
		model(0),
		level(0),
		decisionStack(0),
		dl2ds(0),
		varGroup(0),
		activity(0),
		watches(),
		forcing(),
		ca(0),
		clauseDatabase()
	  {}
	};

	// A helper struct to compare two clauses wrt. their LBD value & clause length.
	template <typename T>
	  struct LBDOrder
	  {
	  explicit LBDOrder(const std::vector<T>& c) :
		act(c),
		  core(NULL)
		{}
	  
		bool operator()(const uint32_t& x, const uint32_t& y) const
		{
		  return CompareLBD(act[x], act[y]);
		}

		void SetCore(Core* c) { core = c; }

		bool CompareLBD(const CRef& c1, const CRef& c2) const
		{
		  if (core->_ca[c1].Lbd() == core->_ca[c2].Lbd())
			{ return core->_ca[c1].Activity() < core->_ca[c2].Activity(); }
		  return core->_ca[c1].Lbd() > core->_ca[c2].Lbd(); 
		}
		
		const std::vector<T>& act;
		Core* core;
	  };

	template <typename T>
	  struct ActivityOrder
	  {
	  explicit ActivityOrder(const std::vector<T>& c) :
		act(c),
		  core(NULL)
		{}

		bool operator()(const uint32_t& x, const uint32_t& y) const
		{
		  return CompareActivity(act[x], act[y]);
		}

		void SetCore(Core* c) { core = c; }

		bool CompareActivity(const CRef& c1, const CRef& c2) const
		{
		  if (core->_ca[c1].Activity() == core->_ca[c2].Activity())
			{ return core->_ca[c1].size() > core->_ca[c2].size(); }
		  return core->_ca[c1].Activity() < core->_ca[c2].Activity();
		}

		const std::vector<T>& act;
		Core* core;
	  };

  public:

    // Constructor.
  Core(Control* control, Settings* setting, uint32_t id = 0) : 
	_id(id), 
      _control(control),
	  _setting(setting),
	  _ca(1024*1024),
      _emptyClause(false), 
	  _startPtr(1),
	  _endPtr(1),
      _variables(0),
	  _capacity(0),
	  _statistics(),
      _assignment(),
      _level(),
      _activity(),
	  _activityInc(1),
      _polarity(),
      _forcing(),
	  _conflict(),
	  _newUnits(0),
      _watches(),
      _varOrder(),
      _varGroup(), 
	  _noOfVarGroups(1),
      _incVarActivity(1.0),
      _decisionLevel(0),
      _decisionStack(),
	  _stackQueue(),
      _dl2ds(),
      _dsEndIndex(1),
      _dsImplIndex(1),
      _modeDS(), 
	  _modeDSForVariables(),
	  _analyzeStack(),
	  _seenToClear(),
	  _touchedToClear(),
	  _seen(),
	  _touched(),
	  _conflictLiteral(0),
	  _conflictClause(),
	  _failedAssumptions(),
	  _learnedClausesLimit(20000),
      _decVarSign(),
	  _candidates(),
      _basicOps(0), 
      _model(),
      _clauseDatabase(),
      _solverState(),
	  _preprocessor(NULL),
	  _deleted()
    {
      // Consistency check.
      assert(_control != nullptr);
	  assert(_setting != nullptr);

	  // Initialize first varorder level
	  VarHeap<helper::DescendingOrder, double >* varheap = new VarHeap<helper::DescendingOrder, double >( helper::DescendingOrder<double>(_activity) );
	  _varOrder.push_back(varheap);
	  _modeDS.push_back(CACHEANDTOGGLE);
	  _decVarSign.push_back(false);
	  _stackQueue.InitSize(5000); // magic number of watching last 5000 decisions
    }

    // Destructor.
    ~Core(void)
    {
	  // Delete the solver state.
      DeleteStatus();

	  for( uint32_t v = 0; v != _noOfVarGroups; ++v )
		{
		  assert( _varOrder[v] != NULL ); delete _varOrder[v];
		}
    }


	// Set preprocessor
	void SetPreprocessor(Preprocessor* prepro)
	{ _preprocessor = prepro; }
 
    // Returns the number of variables for which memory has been reserved.
    uint32_t VariableCapacity(void) const 
	{ return _capacity; }

    // Returns the maximial variable index which is used.
	// Does not necessary correspond with "variableCapacity()"
    uint32_t Variables(void) const
	{ return _variables; }

    // Returns the current number of clauses within the clause database. 
    uint32_t Clauses(void) const
    {
      uint32_t cl(0);
      for (uint32_t v = 1; v <= _variables; ++v)
		{
		  uint32_t lit(v << 1); 
		  for (uint32_t r = 0; r < 2; ++r)
			{
			  lit ^= 1;
			  const std::vector<Watcher>& watches(_watches[lit]);
			  size_t size(watches.size()); 
			  for( size_t p = 0; p < size; ++p )
				{
				  // Binaries are only stored in watchlists
				  if( watches[p].IsBinary() )
						{ ++cl; }
				}
			}
		}

      return ((cl >> 1) + (uint32_t)_clauseDatabase.size());
    }

    // Returns the current number of literals within the clause database. 
    uint32_t Literals(void) const
    {
      uint32_t li(0); 
      for (uint32_t v = 1; v <= _variables; ++v)
		{
		  uint32_t lit(v << 1); 
		  for (uint32_t r = 0; r < 2; ++r)
			{
			  lit ^= 1;
			  const std::vector<Watcher>& watches(_watches[lit]);
			  size_t size = watches.size(); 
			  for( size_t p = 0; p < size; ++p )
				{
				  // Binary literals are counted twice
				  if( watches[p].IsBinary() )
					{ ++li; }
				}
			}
		}
      for (uint32_t c = 0; c < _clauseDatabase.size(); ++c)
		{
		  li += _ca.GetClause(_clauseDatabase[c])->size();
		}
      return li; 
    }

    // Returns the number of decisions made so far.
    uint32_t Decisions(void) const 
	{ return _statistics.decisions; }

    // Returns the number of BCP operations performed so far.
    uint32_t Bcps(void) const 
	{ return _statistics.bcps; }

    // Returns the number of implications found so far.
    uint64_t Implications(void) const 
	{ return _statistics.implications; }

    // Returns the number of conflicts encountered so far.
    uint32_t Conflicts(void) const 
	{ return _statistics.totalConflicts; }

    // Returns the number of restarts performed so far.
    uint32_t Restarts(void) const 
	{ return _statistics.restarts; }
	
	uint32_t BlockedRestarts(void) const 
	{ return _statistics.blockedRestarts; }

    // Returns the number of database simplifications performed so far.
    uint32_t Simplifications(void) const 
	{ return _statistics.simplifications; }

    // Returns the number of binary clauses deduced due to "Lazy Hyper Binary Resolution".
    uint32_t Lhbr(void) const
	{ return _statistics.lhbr; }

	uint32_t UsedVariables(void) const
	{ return _statistics.usedVariables; }
	uint32_t CurrentBinaryClauses(void) const 
	{ return _statistics.currentBinaryClauses>>1; }
	uint32_t CurrentTernaryClauses(void) const 
	{ return _statistics.currentTernaryClauses; }
	uint32_t CurrentNaryClauses(void) const 
	{ return _statistics.currentNaryClauses; }

    // Returns the number of unit clauses deduced due to conflict analysis.
    uint32_t LearntUnitClauses(void) const 
	{ return _statistics.totalLearntUnitClauses; }

    // Returns the number of binary clauses deduced due to conflict analysis.
    uint32_t LearntBinaryClauses(void) const 
	{ return _statistics.totalLearntBinaryClauses; }

    // Returns the number of ternary clauses deduced due to conflict analysis.
    uint32_t LearntTernaryClauses(void) const
	{ return _statistics.totalLearntTernaryClauses; }

    uint32_t MinimizedLiterals(void) const
	{ return _statistics.minimizedLiterals; }
	
    // Returns the number of synchronizations performed so far. 
    uint32_t Synchronizations(void) const 
	{ return _statistics.synchronizations; }

	// Returns the number of inprocessings steps during solving main routine
	uint32_t Inprocessings(void) const 
	{ return _statistics.inprocessings; }

    // Returns the progress after the search process has been stopped due to reaching the limit wrt. synchronizations.
    double Progress(void) const 
	{ return _statistics.progress; }

    // Returns the average "Literals Blocks Distance" of all conflict clauses deduced so far.
    double AvgLBD(void) const 
	{ return _statistics.AvgLBD(); }
 
    // Returns the average length of all conflict clauses deduced so far.
    double AvgCCLength(void) const 
	{ return _statistics.AvgCCLength(); }

    // Returns the solver's average decision level before backtracking.
    double AvgDL(void) const 
	{ return _statistics.AvgDL(); }

    // Returns the average number of decision levels cleared during conflict analysis. 
    double AvgDLclearedCA(void) const 
	{ return _statistics.AvgDLclearedCA(); }

    // Returns the average number of variables getting unassigned during conflict analysis. 
    double AvgVarsUnassignedCA(void) const 
	{ return _statistics.AvgVarsUnassignedCA(); }

    // Returns a reference to either the satisfying variable assignment (after calling "solve()") or 
    // the set of currently assigned variables (after calling "deduceAssumptions()"). Example:
    // model[17] =  0 --> x17 = unassigned
    // model[17] = 35 --> x17 = false 
    // model[17] = 34 --> x17 = true
    // In case neither "solve()/maxSolve()" nor "deduceAssumptions()" has been called beforehand, the 
    // vector contains invalid data. 
    const std::vector<uint32_t>& Model(void) const 
	{ return _model; }

	void SetModel(const std::vector<uint32_t>& model)
	{
	  assert( _model.size() == model.size() );
	  _model = model;
	}

	// Write trivial assignemnts ( assignemnts on decision level 0 ) into model
	void TrivialAssignment(void)
	{ 
      // Has to be executed on decision level 0.
      assert(_decisionLevel == 0);

      // Clear "_model".
      for (uint32_t v = 1; v <= _variables; ++v)
		{ _model[v] = 0; }

      // Copy all assigned variables to "_model".
      for (uint32_t d = 1; d < _dsEndIndex; ++d)
		{ 
		  _model[_decisionStack[d] >> 1] = _decisionStack[d]; 
		}
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
    void SetDecisionStrategy(DecisionStrategy val, uint32_t group)  
	{ 
	  if( group >= _noOfVarGroups )
		{ UpdateVarorder(group); }
	  _modeDS[group] = val; 
	}

	// Like "setDecisionStrategy()" for a specific variable instead of a group
    void SetDecisionStrategyForVariable(DecisionStrategy val, uint32_t var) 
	{ 
	  assert(var <= _variables); 
	  _modeDSForVariables[var] = val; 
	}

	// Set initial polarity of the variable "var" to "pol"
	void SetPolarity(uint32_t var, bool pol)
	{
	  // Flip "pol", since, FALSE in _polarity sets "var" to TRUE, and vice versa
	  _polarity[var] = !pol;
	}

    // Sets the maximum variable index to "max". 
    void SetMaxIndex(uint32_t max) 
	{
	  assert(max > 0); 
	  if( max < _variables ) 
		{ return; }
	  UpdateDataStructures(max); 
	}

    // Sets the group of variable "var" to "grp". See "setDecisionStrategy()" for more details. 
    void SetVarGroup(uint32_t var, uint32_t grp) 
	{ 
	  assert(var > 0 && var <= _variables ); 

	  if( grp >= _noOfVarGroups )
		{ UpdateVarorder(grp); }
	  _varGroup[var] = grp; 
	}

	// Returns whether variable was already deleted in pre-/in-processing
	bool IsDeleted(uint32_t var) const 
	{ 
	  assert( var <= _capacity); 
	  return _deleted[var]; 
	}

	// Returns whether an empty clause is deduced
	bool EmptyClause(void) const 
	{ return _emptyClause; }

    // Adds a clause to the clause database. Returns FALSE if the CNF formula is unsatisfiable,
    // otherwise TRUE will be returned. Assumes that the solver is on decision level 0 and that 
    // "clause" is not empty. Furthermore, all literals have to be encoded as follows, having 
    // variable indices greater 0:
    //  x3 <--> (3 << 1)     = 6
    // -x3 <--> (3 << 1) + 1 = 7
    // All clauses inserted into the clause database using "addClause()" are assumed to belong to 
    // the original CNF formula (having a "Literals Blocks Distance" of 1). 
    // IN THE MULTI-THREADED MODE, "maxSetIndex()" HAS TO BE CALLED BEFOREHAND.
    bool AddClause(std::vector<uint32_t>& clause, uint32_t lbd = 1)
    {
      // What about the empty clause?
      if (_emptyClause)
		{ assert( _control->GetExtendedResult() == ANTOM_UNSAT ); return false; }

      // Are we really on decision level 0?
      assert(_decisionLevel == 0); 

      // If "clause" is empty, we might have a problem.
      assert(!clause.empty()); 

      // Sort "clause" to speedup the checks below.
      std::sort(clause.begin(),clause.end());

      // "Shifted" variable indices are assumed to be greater 1.
      assert(clause.front() > 1); 

	  // FIXME Does only work with clause sorting!
	  if( (clause.back()>>1) > _variables )
		{ UpdateDataStructures(clause.back()>>1); }


//	  std::cout << __func__ << " ";
//      for (uint32_t c = 0; c < clause.size(); ++c)
//		{
//		  std::cout << helper::Lit(clause[c]) << " [" << IsAssigned(clause[c]) << "] ";
//		}
//	  std::cout << std::endl;


	  // Assume that variable indices was initialized before using "setMaxIndex()"
	  // FIXME Does only work with clause sorting!
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

		  // added variable is either assigned or not deleted
		  assert( !_deleted[l>>1] || _assignment[l] || _assignment[l^1] );

		  // "clause" satisfied by "l"? Do we have a tautological clause?
		  // FIXME Does only work with clause sorting!
		  if (_assignment[l] || (l ^ 1) == lit)
			{ clause.resize(0); return true; }

		  // Do we have to take the current literal into account?
		  // FIXME Does only work with clause sorting!
		  //if (!_assignment[l ^ 1])
		  if (!_assignment[l ^ 1] && l != lit)
			{ clause[size++] = l; lit = l; }
		}

      // Do we have an empty clause? CNF formula unsatisfiable?
      if (size == 0)
		{
		  if (_setting->verbosity > 0)
			{
			  std::cout << "c empty clause " << std::endl;
			}
		  _emptyClause = true; 
		  _control->SetExtendedResult(ANTOM_UNSAT); 
		  return false; 
		}
      
      // Resize "clause" (necessary for the multi-threaded mode to work correctly).
      clause.resize(size); 

      // Do we have a unit clause?
      if (size == 1)
		{
		  // Push the unit literal as a "fake" implication onto the decision stack.
		  AddImplication(clause[0]);
		  
		  // What about the effects of this implication?
		  if (!Deduce())
			{
			  _emptyClause = true;
			  _control->SetExtendedResult(ANTOM_UNSAT);
			  return false; 
			}

		  // Everything went fine.
		  return true; 
		}

      // Do we have a binary clause?
      if (size == 2)
		{
		  // Get both watched literals.
		  uint32_t wl0(clause[0]);
		  uint32_t wl1(clause[1]);
		  
		  // Update "_activity".
		  IncreaseActivity( (wl0 >> 1), false );
		  IncreaseActivity( (wl1 >> 1), false );

		  // choose watch list depending on whether we implicitely deduce binaries first
		  const std::vector<Watcher>& bWatches = _watches[wl0];
		  size_t bSize(bWatches.size()); 
		  for (size_t bW = 0; bW < bSize; ++bW)
			{
			  if (bWatches[bW].IsBinary() && bWatches[bW].GetSecondLit() == wl1)
				{ return true; }
			}

		  _statistics.currentBinaryClauses += 2;
		  _statistics.staticLength += 2;
		  ++_statistics.staticClauses;

		  // Update "_watches".
		  _watches[wl0].push_back( Watcher(wl1, (lbd > 1)) );
		  _watches[wl1].push_back( Watcher(wl0, (lbd > 1)) );

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

	  // Copy all literals of "clause" to "cptr" and update "_activity". 
      for (uint32_t l = 0; l < size; ++l)
		{ 
		  IncreaseActivity( (clause[l] >> 1), false ); 
		}

	  ++_statistics.staticLength += size;
	  ++_statistics.staticClauses;

	  CRef cr = _ca.Alloc(clause, lbd, _activityInc, size );
	  
      // Update "_watches" and "_clauseDatabase"
	  AttachClause(cr);

      // Everything went fine.
      return true;
    }

    // Adds a clause to the clause database. Returns FALSE if the CNF formula is unsatisfiable,
    // otherwise TRUE will be returned. Assumes that the solver is on decision level 0 and that 
    // "lits != NULL" and "num > 0" holds. Furthermore, all literals have to be encoded as follows, 
    // having variable indices greater 0:
    //  x3 <--> (3 << 1)     = 6
    // -x3 <--> (3 << 1) + 1 = 7
    // All clauses inserted into the clause database using "addClause()" are assumed to belong to 
    // the original CNF formula (having a "Literals Blocks Distance" of 1). 
    // NOTE, THAT THIS VARIANT OF "addClause()" REQUIRES THAT
    // 1) THE CLAUSE TO BE ADDED DOES NOT CONTAIN MULTIPLE COPIES OF THE SAME LITERAL,
    // 2) THE CLAUSE TO BE ADDED IS NOT A TAUTOLOGICAL ONE, AND
    // 3) "maxSetIndex()" HAS BEEN CALLED BEFOREHAND.
    bool AddClause (uint32_t* lits, uint32_t num, uint32_t lbd = 1)
    {
 
      // What about the empty clause?
      if (_emptyClause)
		{ assert( _control->GetExtendedResult() == ANTOM_UNSAT ); return false; }
	  
      // Are we really on decision level 0?
      assert(_decisionLevel == 0); 

      // Everything ok wrt. "lits" and "num"?
      assert(lits != NULL && num > 0); 

	  /*
	  std::cout << __func__ << " ";
      for (uint32_t c = 0; c < num; ++c)
		{
		  std::cout << Lit(lits[c]) << " [" << isAssigned(lits[c]) << "] ";
		}
	  std::cout << std::endl;
	  */

      uint32_t lit(0);
      uint32_t size(0); 

	  // sort clause for efficent tautology checks
	  std::sort(lits, lits+num);

      // Check whether the clause to be added is already satisfied. 
      // By the way, eliminate literals evaluating to FALSE.
      for (uint32_t c = 0; c < num; ++c)
		{
		  // Get the next literal.
		  uint32_t l(lits[c]);

		  // Consistency check.
		  assert((l >> 1) <= _variables); 

		  // added variable is either assigned or not deleted
		  assert( !_deleted[l>>1] || _assignment[l] || _assignment[l^1]);

		  // "clause" satisfied by "l"? Do we have a tautological clause?
		  // FIXME Does only work with clause sorting!
		  if (_assignment[l] || (l ^ 1) == lit)
			{ 
			  return true; 
			}

		  // Do we have to take the current literal into account?
		  // FIXME Does only work with clause sorting!
		  if (!_assignment[l ^ 1] && l != lit)
			{ 
			  lits[size++] = l;
			  lit = l; 
			}
		}

	 
      // Do we have an empty clause? CNF formula unsatisfiable?
      if (size == 0)
		{ 
		  _emptyClause = true; 
		  _control->SetExtendedResult(ANTOM_UNSAT); 
		  return false; 
		}
      
      // Do we have a unit clause?
      if (size == 1)
		{
		  // Push the unit literal as a "fake" implication onto the decision stack.
		  AddImplication(lits[0]);
		  // What about the effects of this implication?
		  if (!Deduce())
			{
			  _emptyClause = true; 
			  _control->SetExtendedResult(ANTOM_UNSAT);
			  return false; 
			}

		  // Everything went fine.
		  return true; 
		}

      // Do we have a binary clause?
      if (size == 2)
		{

		  // Get both watched literals.
		  uint32_t wl0(lits[0]);
		  uint32_t wl1(lits[1]);
		  
		  // Update "_activity".
		  IncreaseActivity( (wl0 >> 1), false);
		  IncreaseActivity( (wl1 >> 1), false);

		  // Check whether we have that particular binary clause already.
		  
		  // choose watch list depending on whether we implicitely deduce binaries first
		  const std::vector<Watcher>& bWatches = _watches[wl0];

		  
		  size_t bSize(bWatches.size()); 
		  for (size_t bW = 0; bW < bSize; ++bW)
			{
			  if (bWatches[bW].IsBinary() && bWatches[bW].GetSecondLit() == wl1)
				{ return true; }
			}
		  
		  _statistics.currentBinaryClauses += 2;
		  _statistics.staticLength += 2;
		  ++_statistics.staticClauses;

		  // Update "_watches".
		  _watches[wl0].push_back( Watcher(wl1, (lbd > 1)) );
		  _watches[wl1].push_back( Watcher(wl0, (lbd > 1)) );

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
      
      // Consistency check.
      assert(size <= num); 

      // Update "_activity". 
      for (uint32_t l = 0; l < size; ++l)
		{ 
		  IncreaseActivity( (lits[l] >> 1), false);
		}

	  CRef cr(_ca.Alloc(lits, lbd, _activityInc, size));

	  ++_statistics.staticLength += size;
	  ++_statistics.staticClauses;
  
      // Update "_watches" and "_clauseDatabase"
	  AttachClause(cr);

      // Everything went fine.
      return true;
    }

    // Performs unit propagation, taking the current CNF and the specified assumptions into 
    // account. Returns FALSE if a conflict occurs, otherwise the return value is TRUE. 
    bool DeduceAssumptions(const std::vector<uint32_t>& assumptions)
    {
      // What about the empty clause?
      if (_emptyClause)
		{ assert( _control->GetExtendedResult() == ANTOM_UNSAT); return false; }

      // If there are no variables, we don't have to perform unit propagation.
      if (_variables == 0) 
		{ return true; }

      // "deduceAssumptions()" has to be executed on decision level 0.
      assert(_decisionLevel == 0);

      // Clear "_model".
	  std::fill(_model.begin(), _model.end(), 0 );

	  // Update "_varOrder".
	  for( uint32_t v = 0; v != _noOfVarGroups; ++v )
		{ 
		  _varOrder[v]->clear();
		  _varOrder[v]->resize(_variables); 
		}

      // Push all assumptions onto the decision stack and perform unit propagation.
      if (!SetAssumptions(assumptions))
		{ _control->SetDone(); return false; }

      // Copy all assigned variables to "_model".
      for (uint32_t d = 1; d < _dsEndIndex; ++d)
		{ 
		  _model[_decisionStack[d] >> 1] = _decisionStack[d]; 
		}

      // Backtrack to decision level 0.
      if (_decisionLevel != 0)
		{ 
		  Backtrack(0); 
		}

	  // Everything went fine.
	  return true; 
	}

	// Returns learnt conflict clauses and their "LBD"-value containing clauses with index <= "lastindex"
	// This method assumes that all learnt clauses can be perceived
	std::vector< std::pair<std::vector< uint32_t >, uint32_t > > GetConflictClauses(void) const
	{
	  assert( _decisionLevel == 0 );

	  std::vector< std::pair<std::vector< uint32_t >, uint32_t > > result;

	  std::vector< uint32_t > candidate;
	  
	  for( uint32_t i = 0; i != _clauseDatabase.size(); ++i )
		{
		  // Consider learnt clauses
		  if( _ca[_clauseDatabase[i]].IsLearned() )
			{
			  candidate.clear();
			  const Clause& clause(_ca[_clauseDatabase[i]]);
			  uint32_t length(clause.size());
			  uint32_t pos(0);
			  for( ; pos < length; ++pos )
				{
				  assert( ( clause[pos] >> 1 ) <= _variables );
				  assert( !_deleted[clause[pos]>>1]);
				  // Skip satisfied conflict clauses
				  if( _assignment[clause[pos]] )
					{ break; }
				  // Skip unsatisfied literals
				  else if( !_assignment[clause[pos]^1] )
					{ 
					  candidate.push_back(clause[pos]); 
					}
				}

			  if(pos == length )
				{ result.push_back(std::make_pair(candidate, clause.Lbd() ) ); }
			}
		}

	  candidate.resize(2);
	  // Now add binary clauses
	  for( uint32_t v = 1; v <= _variables; ++v )
		{
		  for( uint32_t literal = (v<<1); literal != ((v<<1)+2); ++literal )
			{
			  const std::vector< Watcher >& watches = _watches[literal];
			  for( uint32_t i = 0; i != watches.size(); ++i )
				{
				  if( watches[i].IsLearnedBinary() && ( (watches[i].GetSecondLit()>>1) > v ) )
					{	
					  assert( !_deleted[watches[i].GetSecondLit()>>1] );
					  candidate[0] = literal;
					  candidate[1] = watches[i].GetSecondLit();
					  result.push_back( std::make_pair( candidate, 2 ) );
					}
				}
			}
		}
	  return result;
	}

	// Solves the current CNF formula, taking the specified assumptions into account. The assumptions have to be encoded 
	// in the same way as the literals of a clause (see "addClause()"). The return values are SAT/UNSAT/UNKNOWN. In case 
	// "limit" is not equal to 0, the search process will be stopped after "limit" synchronizations. 
	// NOTE, THAT THE CURRENT VERSION OF "solve()" ASSUMES THAT ALL THREADS HAVE THE SAME SET OF ASSUMPTIONS AND THAT THE 
	// MULTI-THREADED MODE FOLLOWS AN ALGORITHM PORTFOLIO APPROACH. 
	uint32_t Solve(const std::vector<uint32_t>& assumptions, uint64_t limit = 0)
	{
	  // restart factors
	  uint32_t nxtRestart(_statistics.totalConflicts + 256);
	  uint32_t lubyNumber(0);

	  // minisat restart parameters
	  double learntinc(1.1);
	  double adjustconflicts(100.0);
	  uint32_t adjustcounts(100);
	  double adjustinc(1.5);

	  if( _setting->simplifyStrategy == MINISAT )
		{
		  double learntfactor(1.0/3.0);
		  _learnedClausesLimit = static_cast<uint32_t>(Clauses()*learntfactor); 
		}

	  uint32_t nextinpro(0);
	  if( _setting->restartStrategy == LUBY && ( (_setting->lubyShift>>1) < 12 ) )
		{
		  nextinpro = (1<<(12-(_setting->lubyShift>>1)))-1;
		}
	  else if ( _setting->restartStrategy == GLUCOSE && ( (_setting->lubyShift) < 13 ))
		{
		  nextinpro = (1<<(13-(_setting->lubyShift)))-1;
		}

	  uint32_t rootLevel(0);
	  _startPtr = 1; 
	  _endPtr = 1; 
	  uint32_t result(ANTOM_UNKNOWN);
	  bool limitSpecified(limit != 0);

	  // needed for timeout specification
	  uint32_t lastConflictCount(0);

	  // Initilaize start time
	  struct rusage resources;

	  getrusage(RUSAGE_SELF, &resources);
	  double timeS  = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
	  _control->SetStartTime(timeS);

	  // If there are no variables, the CNF formula is satisfiable by definition.
	  if (_variables == 0) 
		{
		  result = ANTOM_SAT; 
		  _control->SetExtendedResult(ANTOM_SAT);
		  _control->SetDone(); 
		  return Synchronize(result);
		}

	  // What about the empty clause?
	  if (_emptyClause)
		{ 
		  result = ANTOM_UNSAT; 
		  assert( _control->GetExtendedResult() == ANTOM_UNSAT); 
		  _control->SetDone(); 
		  return Synchronize(result);
		}

	  _control->SetExtendedResult(ANTOM_UNKNOWN);

	  // Update "_varOrder".
	  for( uint32_t v = 0; v != _noOfVarGroups; ++v )
		{ 
		  _varOrder[v]->clear();
		  _varOrder[v]->resize(_variables); 
		  _decVarSign[v] = false;
		}

	  // The SAT solving process has to be started on decision level 0.
	  assert(_decisionLevel == 0);
	  assert(_dsEndIndex == _dsImplIndex);

	  // Push all assumptions onto the decision stack and perform unit propagation.
	  if (!SetAssumptions(assumptions))
		{ 
		  result = ANTOM_UNSAT; 
		  _control->SetDone(); 
		  return Synchronize(result);
		}

    
	  // Initialization.
	  rootLevel = _decisionLevel;
	        
	  if(limitSpecified)
		{ limit += _statistics.synchronizations; }

	  _statistics.usedVariables = 0;
	  // Initialize "_varOrder".
	  for (uint32_t v = _variables; v > 0; --v)
		{
		  // Empty watch list? => Variable is not used
		  if ( !_deleted[v] && !_assignment[v << 1] && !_assignment[(v << 1) ^ 1])
			{ 
			  _varOrder[_varGroup[v]]->insert(v); 
			  ++_statistics.usedVariables;
			}
		}

	  // Re-initialize some of our status variables.
	  _basicOps      = 0; 
	  _statistics.localLBD      = 0; 

	  // Consistency check.
	  assert(_dsEndIndex == _dsImplIndex); 

	  // The main SAT solving loop.
	  while ( Decide() )
		{
		  // Perform "Boolean Constraint Propagation".
		  while ( !Deduce() )
			{
			  // Unresolvable conflict?
			  if (_decisionLevel <= rootLevel)
				{ 
				  result = ANTOM_UNSAT;
				  // Unsatisfiable CNF formula? If the decision level is equal to 0, we either don't 
				  // have assumptions or all assumptions are unit clauses (or at least forced by unit 
				  // clauses). Anyway, the assumptions don't play a role and the CNF formula ist not 
				  // only unsatisfiable under assumptions, but unsatisfiable at all.
				  if (_decisionLevel == 0)
					{ 
					  _emptyClause = true; 
					  _control->SetExtendedResult(ANTOM_UNSAT); 
					  _control->SetDone(); 
					  return Synchronize(result);
					}

				  // Backtrack to decision level 0.
				  Backtrack(0); 

				  // At this point, the CNF formula is unsatisfiable under assumptions.
				  assert( _control->GetExtendedResult() != ANTOM_UNSAT );

				  _control->SetExtendedResult(ANTOM_UNSAT_WITH_ASSUMPTION); 
				  _control->SetDone(); 
				  return Synchronize(result);
				}

			  _stackQueue.Push(_dsEndIndex);

			  // block the next restart
			  if (_setting->restartBlocking && (_dsEndIndex > (1.2 * _stackQueue.GetAvg())))
				{
				  ++_statistics.blockedRestarts;
				  _statistics.localLBD = 0;
				  _statistics.localConflicts = 0;
				}

			  // Analyze the current conflict and backtrack.
			  Analyze();
      
			  // What about synchronizing with other threads?
			  if (_basicOps >= 6000000)
				{
				  // Increment "_synchronizations".
				  ++_statistics.synchronizations; 
		  
				  // Reset "_basicOps".
				  _basicOps = 0; 

				  // Initialization.
				  _endPtr = _dsEndIndex; 
				  if (_decisionLevel != 0)
					{ _endPtr = _dl2ds[1]; }

				  // Consistency check.
				  assert(_startPtr <= _endPtr); 

				  // Export all assignments made on decision level 0 to the "Control" object.
				  _control->ExportUnitClauses(&_decisionStack[0], _startPtr, _endPtr, _id); 

				  if( Synchronize(result) != ANTOM_UNKNOWN)
					{
					  return result;
					}
				}

			  if( _setting->simplifyStrategy == MINISAT )
				{
				  if( --adjustcounts == 0 )
					{
					  adjustconflicts *= adjustinc;
					  adjustcounts = (int)adjustconflicts;
					  _learnedClausesLimit = static_cast<uint32_t>(_learnedClausesLimit*learntinc);
					  EstimateProgress();
					}
				}

			  // Should we check for another restart?
			  if (_statistics.totalConflicts >= nxtRestart)
				{
				  //std::cout << "tC: " << _statistics.totalConflicts << " nxtre: " << nxtRestart << std::endl;
				  bool didrestart(false);

				  // Do we use a glucose-like restart strategy?
				  if (_setting->restartStrategy == GLUCOSE)
					{
					  // Initialization. 
					  double global((double) _statistics.globalLBD / _statistics.totalConflicts);
					  //double local((double) _statistics.localLBD / (1 << _setting->lubyShift));
					  double local((double) _statistics.localLBD / _statistics.localConflicts); 

					  // What about performing a restart?
					  if ((0.8 * local) > global) 
						{
						  // Increment "_restarts".
						  ++_statistics.restarts; 
						  didrestart = true;
			  
						  // Backtrack to decision level 0.
						  if (_decisionLevel > 0)
							{ Backtrack(0); }
			  
						  // Flip "_decVarSign" (in case the decision strategy has been set to mode 0).
						  for( uint32_t v = 0; v != _noOfVarGroups; ++v )
							{
							  if (_modeDS[v] == CACHEANDTOGGLE)
								{ _decVarSign[v] = _decVarSign[v] ^ true; }
							}
						}
	      
					  // Update "nxtRestart".
					  nxtRestart += 1 << _setting->lubyShift; 
		      
					  // Reset "_localAvgLBD".
					  _statistics.localLBD = 0;
					  _statistics.localConflicts = 0;
					}
				  else if (_setting->restartStrategy == LUBY)
					{
					  // Increment both, "_restarts" and "lubyNumber".
					  ++_statistics.restarts; 
					  ++lubyNumber; 
					  didrestart = true;
	      
					  // Backtrack to decision level 0.
					  if (_decisionLevel > 0)
						{ Backtrack(0); }
		      
					  // Update "nxtRestart".
					  nxtRestart += helper::Luby(lubyNumber) << _setting->lubyShift; 
		      
					  // Flip "_decVarSign" (in case the decision strategy has been set to mode 0).
					  for( uint32_t v = 0; v != _noOfVarGroups; ++v )
						{
						  if (_modeDS[v] == CACHEANDTOGGLE)
							{ _decVarSign[v] = _decVarSign[v] ^ true; }
						}
					} 

				  // Perform inprocessing... wo only perfrom inprocessing when we are on decision level 0
				  if (didrestart)
					{
					  if( _setting->doInprocessing && ( (_statistics.restarts&nextinpro) == nextinpro ) )
						{
						  assert( _decisionLevel == 0 );
						  if( _preprocessor->Preprocess(INPROCESS) == ANTOM_UNSAT )
							{ 
							  _emptyClause = true; 
							  result = ANTOM_UNSAT; 
							  _control->SetExtendedResult(ANTOM_UNSAT); 
							  _control->SetDone(); 
							  return Synchronize(result);
							}
						}
					}
				}
			}
		
		  // after deduction...

		  // What about simplifying the clause database?
		  if ( _statistics.learnedClauses >= _learnedClausesLimit ) 
			{
			  // Simplify the clause database, in particular perform conflict clause deletion.
			  Simplify(true);

			  if( _setting->simplifyStrategy == ANTOM )
				{
				  // Update "learntClausesLimit".
				  _learnedClausesLimit = static_cast<uint32_t>(_learnedClausesLimit*1.1);
				}
			}
		
		  // What about the time limit?
		  if ( (_statistics.totalConflicts-lastConflictCount) > 512 )
			{
			  if( _control->ReachedLimits() )
				{
				  assert( _control->GetExtendedResult() != ANTOM_UNSAT);
				  result = ANTOM_UNKNOWN; 
				  _control->SetExtendedResult(ANTOM_UNKNOWN); 
				  _control->SetDone(); 
				  return Synchronize(result);
				}
			  lastConflictCount = _statistics.totalConflicts;
			}

		  // Do we have to stop the search process?
		  if (limitSpecified && _statistics.synchronizations > limit)
			{
			  // Determine the progress.
			  EstimateProgress();

			  // Backtrack to decision level 0.
			  if (_decisionLevel != 0)
				{ Backtrack(0); }
	      
			  // Return unknown.
			  assert( _control->GetExtendedResult() != ANTOM_UNSAT );
			  result = ANTOM_UNKNOWN; 
			  _control->SetExtendedResult(ANTOM_UNKNOWN);
			  _control->SetDone(); 
			  return Synchronize(result);
			}
	  
		  // What about the decision level?
		  if (_decisionLevel < rootLevel)
			{
			  // Push all assumptions onto the decision stack and perform unit propagation.
			  if (!SetAssumptions(assumptions))
				{ 
				  result = ANTOM_UNSAT;
				  _control->SetDone(); 
				  return Synchronize(result);
				}
	      	      
			  // Update "rootLevel".
			  rootLevel = _decisionLevel;
			}
		}

	  // At this point, there's no variable left to serve as a decision variable. 
	  // So, we have found a satisfying variable assignment.
#ifndef NDEBUG
	  for( uint32_t v = 0; v != _noOfVarGroups; ++v )
		{ assert(_varOrder[v]->empty()); }
#endif
      
	  // Clear "_model".
	  std::fill(_model.begin(), _model.end(), 0);
	  // Store the satisfying variable assignment in "_model".
	  for (uint32_t d = 1; d < _dsEndIndex; ++d)
		{
		  _model[_decisionStack[d] >> 1] = _decisionStack[d];
		}


	  /*
	  for (uint32_t v = 1; v <= _variables; ++v)
		{

		  if( _model[v] == 0 && !_deleted[v])
			{
			  // Push fake assignment to model
			  _model[v] = (_model[v]<<1);
			}
		}
	  */

	  // Add assignment of replaced variables and resolved literals to model
	  if( _setting->doPreprocessing == PREPROCESS || _setting->doInprocessing )
		{ _preprocessor->ExtendModel(); }
      
	  // Backtrack to decision level 0.
	  if (_decisionLevel != 0)
		{ Backtrack(0); }
      
	  // Store the result and terminate execution.
	  assert(_control->GetExtendedResult() != ANTOM_UNSAT);
	  result = ANTOM_SAT;
	  _control->SetExtendedResult(ANTOM_SAT);
	  _control->SetDone(); 
	  return Synchronize(result);
	}

	uint32_t Synchronize(uint32_t& currentResult)
	{
#ifdef PARALLEL		  
#pragma omp barrier
#endif
	  // CNF formula solved?
	  if (_control->Done())
		{ 		 		    
		  // Backtrack to decision level 0 (performed by threads that get aborted only).
		  if (_decisionLevel > 0)
			{ Backtrack(0); }

		  assert(currentResult == _control->GetExtendedResult() || _control->GetExtendedResult() == ANTOM_UNSAT_WITH_ASSUMPTION);
		  return currentResult; 
		}

	  // Initialization.
	  uint32_t* unitClauses(NULL); 

	  // Get a pointer to the set of unit clauses found by other threads.
	  unitClauses = _control->ImportUnitClauses(_id);
	  assert(unitClauses != NULL); 

	  // Update "_startPtr".
	  _startPtr = _endPtr; 

#ifdef PARALLEL		  
#pragma omp barrier 
#endif

	  // Import unit clauses found by other threads.
	  uint32_t pos(0); 
	  while (unitClauses[pos] != 0)
		{
		  // Get the next unit clause.
		  uint32_t lit(unitClauses[pos]); 
		  // Consistency check.
		  assert((lit>>1) <= _variables);
		  assert(lit > 1); 
		  assert( _decisionLevel == 0 );
		      
		  // "lit" incorrectly assigned on decision level 0?
		  if (_assignment[lit ^ 1] && _level[lit >> 1] == 0)
			{
			  // At this point the CNF formula is unsatisfiable.
			  _emptyClause = true; 
			  _control->SetExtendedResult(ANTOM_UNSAT);
			  _control->SetDone(); 
			  // Backtrack to decision level 0 (performed by threads that get aborted only).
			  if (_decisionLevel > 0)
				{ Backtrack(0); }

			  assert(currentResult == _control->GetExtendedResult());
			  return currentResult; 
			}

		  // "lit" currently unassigned or assigned on a decision level greater than 0?
		  if ((!_assignment[lit] && !_assignment[lit ^ 1]) || _level[lit >> 1] > 0)
			{ 
			  // Push "lit" as a "fake" implication onto the decision stack.
			  AddImplication(lit);
		      
			  // What about the effects of this implication?
			  if (!Deduce())
				{ 
				  // At this point the CNF formula is unsatisfiable.
				  _emptyClause = true; 
				  _control->SetExtendedResult(ANTOM_UNSAT);
				  _control->SetDone(); 
				  // Backtrack to decision level 0 (performed by threads that get aborted only).
				  if (_decisionLevel > 0)
					{ Backtrack(0); }
				  assert(currentResult == _control->GetExtendedResult());
				  return currentResult;
				}
			}
	 
		  // Increment "pos".
		  ++pos;
		}
	  assert(currentResult != ANTOM_SAT && currentResult != ANTOM_UNSAT);
	  return ANTOM_UNKNOWN;
	}

	void ClearVariables(uint32_t begin, uint32_t end)
	{
	  assert( begin <= end );

	  _statistics.ClearClauseStatistics();

	  // Consistency checks.
	  assert(_dsEndIndex == _dsImplIndex); 
	  assert(_decisionLevel == 0 );

	  size_t csize(_clauseDatabase.size());
	  // Mark all n-nary clauses with (begin <= variables <= end) as "to delete"
	  for( size_t c = 0; c != csize; ++c )
		{
		  Clause& clause( _ca[_clauseDatabase[c]] );

		  uint32_t size(clause.size());
		  for(uint32_t pos = 0; pos != size; ++pos )
			{
			  if( ( (clause[pos]>>1) >= begin ) && ( (clause[pos]>>1) <= end ) )
				{
			  
				  // The clause should not be a forcing clause
				  assert( _forcing[clause[0] >> 1].NoReason() );
				  assert( _forcing[clause[1] >> 1].NoReason() );
				  _ca[_clauseDatabase[c]].MarkForDeletion();
				  break;
				}
			}
		}

	  // Delete all clauses with (begin <= variables <= end) from watchlists 
	  for( uint32_t v = 1; v <= _variables; ++v )
		{
		  if( v >= begin && v <= end )
			{
			  // Clear some status flags for "v"
			  _deleted[v] = false;
			  _preprocessor->_donttouch[v] = false;

			  // Reset some status values for "v"
			  _activity[v] = 0.0;
			  _modeDSForVariables[v] = CACHEANDTOGGLE;
			  _varGroup[v] = 0;
			  _level[v] = 0;
			  // Clear watchlists of "v"
#ifdef CLEARWATCHES
			  _watches[v<<1].clear();
			  _watches[(v<<1)^1].clear();
#else
			  std::vector<Watcher>().swap(_watches[v<<1]);
			  std::vector<Watcher>().swap(_watches[(v<<1)^1]);
#endif
			}

		  for( uint32_t literal = v<<1; literal < ((v<<1)+2); ++literal )
			{
			  // now clear watch structure
			  std::vector< Watcher >& watches( _watches[literal] );
			  size_t size = watches.size();
			  for( size_t i = 0; i != size; ++i )
				{
				  if( watches[i].IsBinary() )
					{
					  if( ( (watches[i].GetSecondLit()>>1) >= begin ) && ( ( watches[i].GetSecondLit()>>1) <= end ) )
						{ 
						  watches[i--] = watches[--size];
						  watches.pop_back();
						}
					  else
						{
						  ++_statistics.currentBinaryClauses;
						  if( watches[i].IsLearnedBinary() )
							{ 
							  ++_statistics.learnedBinary; 
							  ++_statistics.learnedLength;
							}
						  else
							{ 
							  ++_statistics.staticClauses; 
							  ++_statistics.staticLength;
							}
						}
					}
				  else
					{
					  // Delete all watches of n-nary clauses... will be refilled later
					  watches[i--] = watches[--size];
					  watches.pop_back();
					}
				}
			}
		}

	  _statistics.learnedBinary >>= 1; 
	  _statistics.staticClauses >>= 1; 

	  // Now delete marked n-nary clauses and rebuild watches
	  for( uint32_t c = 0; c != csize; ++c )
		{
		  Clause& clause(_ca[_clauseDatabase[c]]);

		  if( clause.ToDelete() )
			{		
			  _clauseDatabase[c--] = _clauseDatabase[--csize];
			  _clauseDatabase.pop_back();
			}
		  else
			{
			  // Update "_watches".
			  _watches[clause[0]].push_back(Watcher(_clauseDatabase[c], clause[1]));
			  _watches[clause[1]].push_back(Watcher(_clauseDatabase[c], clause[0]));
			  
			  // Update "_learntClauses" if necessary.
			  if (clause.IsLearned())
				{ 
				  ++_statistics.learnedClauses; 
				  _statistics.learnedLength += clause.size();
				}
			  else
				{ 
				  ++_statistics.staticClauses; 
				  _statistics.staticLength += clause.size();
				}
			}
		}

	  _preprocessor->ClearRestoreData(begin, end);

	  // We may have put a deleted variable as fake implication in the decision stack
	  uint32_t i(1);
	  uint32_t j(1);
	  for( ; i < _dsEndIndex; ++i )
		{
		  uint32_t var( _decisionStack[i]>>1 );
		  if( ( var < begin ) || ( var > end ) )
			{
			  _decisionStack[j] = _decisionStack[i];
			  ++j;
			}
		  else
			{
			  _assignment[_decisionStack[i]] = false;
			  assert( !_assignment[_decisionStack[i]^1] );
			  assert( _forcing[var].NoReason());
			}
		}
	  _dsEndIndex = j;
	  _dsImplIndex = j;

	  for( uint32_t v = 1; v <= _variables; ++v )
		{
		  if( !(_watches[v<<1].empty() && _watches[(v<<1)^1].empty()) )
			{ ++_statistics.usedVariables; }
		}

	  // update _variables
	  if( end == _variables )
		{ _variables = begin-1; }
	}

	// Simplifies the clause database by removing
	// -- clauses satisfied due to an assignment made on decision level 0,
	// -- literals evaluating to false due to an assignment made on decision level 0,
	// -- conflict clauses with a high "Literals Blocks Distance" (requires "extended" set to TRUE). 
	void Simplify(bool extended)
	{
	  //std::cout << __func__ << std::endl;

	  // Increment "_simplifications".
	  ++_statistics.simplifications;

	  uint32_t tmpusedvars    = _statistics.usedVariables;
	  uint32_t tmpbinary      = _statistics.currentBinaryClauses>>1;
	  uint32_t tmpternary     = _statistics.currentTernaryClauses;
	  uint32_t tmpnary        = _statistics.currentNaryClauses;

	  _statistics.ClearClauseStatistics();

	  // Consistency check.
	  assert(_dsEndIndex == _dsImplIndex); 

	  // Clear "_forcing" for all implications forced on decision level 0.
	  // By the way, remove all variables assigned on decision level 0 from "_varOrder".
	  uint32_t p(1); 
	  while (p < _dsEndIndex && _level[_decisionStack[p] >> 1] == 0)
		{ 
		  uint32_t var(_decisionStack[p] >> 1); 
		  _forcing[var].ClearReason(); 
		  if (_varOrder[_varGroup[var]]->inHeap(var))
			{ _varOrder[_varGroup[var]]->remove(var); }
		  ++p;
		}

	  // Update "_watches".
	  for (uint32_t v = 1; v <= _variables; ++v)
		{
		  // Initialization.
		  uint32_t pLit(v << 1);
		  uint32_t nLit((v << 1) ^ 1);

		  // Variable "v" assigned on decision level 0?
		  if (_deleted[v] || ( (_assignment[pLit] || _assignment[nLit]) && _level[v] == 0) )
			{
			  // Clear "_watches[pLit]/_watches[nLit]".
#ifdef CLEARWATCHES
			  _watches[pLit].clear();
			  _watches[nLit].clear();
#else
			  std::vector<Watcher>().swap(_watches[pLit]);
			  std::vector<Watcher>().swap(_watches[nLit]);
#endif
			}
		  else
			{
			  ++_statistics.usedVariables;
			  // Shrink "_watches[pLit]/_watches[nLit]".
			  for (uint32_t r = 0; r < 2; ++r)
				{
				  // Flip "pLit" (in order to check both polarities). 
				  pLit ^= 1;

				  std::vector<Watcher>& watches(_watches[pLit]); 
				  size_t wSize = watches.size();

				  // Increment "_basicOps".
				  _basicOps += wSize; 

				  // Shrink "_watches" by removing binary clauses satisfied on decision level 0. 
				  // Furthermore, all elements corresponding to non-binary clauses are removed, too. 
				  for (size_t w = 0; w < wSize; ++w)
					{
					  
						if ( !watches[w].IsBinary() || (_level[watches[w].GetSecondLit()>>1] == 0 && _assignment[watches[w].GetSecondLit()]) )
						{ 
						  watches[w--] = watches[--wSize]; 
						  watches.pop_back(); 
						}
					}

				  // Sort "watches".
				  std::sort(watches.begin(), watches.end(), WatchedSorter()); 

				  // Initialization.
				  uint32_t lit(0);
				  uint32_t size(0);

				  // Remove duplicates from "watches".
				  for (size_t w = 0; w < wSize; ++w)
					{
					  assert(watches[w].IsBinary()); 
					  
					  if (watches[w].GetSecondLit() != lit)
						{ 
						  watches[size++] = watches[w]; 
						  lit = watches[w].GetSecondLit(); 

						  ++_statistics.currentBinaryClauses;
						  if( watches[w].IsLearnedBinary() )
							{ 
							  ++_statistics.learnedBinary; 
							  ++_statistics.learnedLength;
							}
						  else
							{ 
							  ++_statistics.staticClauses; 
							  ++_statistics.staticLength;
							}
						}
					}

				  // Shrink "watches".
				  watches.resize(size); 
				}
			}
		}

	  // Initialization.
	  size_t size(_clauseDatabase.size());

	  assert(_candidates.empty());
  
	  // Increment "_basicOps".
	  _basicOps += size;

	  bool betterLbd = false;

	  // Update "candidates".
	  for (size_t c = 0; c < size; ++c)
		{
		  // Get a pointer to the next clause.
		  Clause& lits(_ca[_clauseDatabase[c]]);
		  betterLbd = false;

		  // "clause" currently not forcing an implication and not marked as "to be deleted"?
		  if ( !_forcing[lits[0] >> 1].ForcedBy(_clauseDatabase[c]) && 
			   !_forcing[lits[1] >> 1].ForcedBy(_clauseDatabase[c]) &&
			   !lits.ToDelete() )
			{	

			  // Initialization.
			  uint32_t p(0);
			  uint32_t q(0);
	      
			  // Check whether "clause" is satisfied by an assignment made on decision level 0 or 
			  // if it contains literals evaluating to false due to an assignment made on decision level 0.
			  uint32_t size(lits.size());
			  for( ; p != size; ++p )
				{
				  // Current literal unassigned or assigned on a decision level greater 0?
				  if ((!_assignment[lits[p]] && !_assignment[lits[p] ^ 1]) || _level[lits[p] >> 1] > 0)
					{
					  // At this point, we have to keep the current literal.
					  lits[q++] = lits[p]; 
					}
				  else
					{
					  // Consistency check.
					  assert((_assignment[lits[p]] || _assignment[lits[p] ^ 1]) && _level[lits[p] >> 1] == 0); 

					  // "clause" satisfied by "lits[p]"?
					  if (_assignment[lits[p]])
						{ break; }

					}
				}

			  // Have we reached the end of the clause?
			  if (p == size)
				{
				  // What about the clause length?
				  if (q > 2)
					{
					  // Update the clause length.
					  lits.SetLength(q);

					  _ca.FreeLiterals(p-q);

					  // Update the "Literals Blocks Distance" if necessary.
					  if (lits.Lbd() > q)
						{
						  betterLbd = true;
						  lits.SetLBD(q);
						}
					}
				  else
					{
					  // At this point we have found a new binary clause to be stored separately.
					  assert(q == 2); 

					  // Check whether we already have that particular binary clause.		      
					  std::vector<Watcher>& watches = _watches[lits[0]];

					  size_t wSize(watches.size());
					  size_t w(0); 
					  for (; w < wSize; ++w)
						{
						  assert(watches[w].IsBinary()); 
						  if (watches[w].GetSecondLit() == lits[1])
							{ break; }
						}
		      
					  // Update "_watches".
					  if (w == wSize)
						{
						  _watches[lits[0]].push_back(Watcher(lits[1], lits.IsLearned()));
						  _watches[lits[1]].push_back(Watcher(lits[0], lits.IsLearned()));
						  _statistics.currentBinaryClauses += 2;
						  if( lits.IsLearned() )
							{
							  ++_statistics.learnedBinary; 
							  ++_statistics.learnedLength;							  
							}
						  else 
							{
							  ++_statistics.staticClauses; 
							  ++_statistics.staticLength;
							}
						}

					  // Mark the clause as "to be deleted".
					  lits.MarkForDeletion();
					}
				}
			  else
				{

				  // Current clause satisfied by an assignment made on decision level 0, 
				  // so let's mark it as "to be deleted".
				  lits.MarkForDeletion();
				}
		  
			  // "clause" not part of the original CNF, LBD value greater than 2, and "extended" set to TRUE?

			  if (extended)
				{
				  if (_setting->simplifyActivity == SIMP_LBD && lits.Lbd() > 2)
					{
					  if (!betterLbd)
						{
						  _candidates.push_back(_ca.GetClause(_clauseDatabase[c]));
						}
					}
				  // Do not delete ternary clauses if used, since activites are not updated
				  else if (_setting->simplifyActivity == SIMP_CONF && lits.Lbd() > 1 && (!_setting->useTernary || q > 3) )
					{
					  _candidates.push_back(_ca.GetClause(_clauseDatabase[c]));
					}
				}
			}
		}

	  _statistics.learnedBinary >>= 1;
	  _statistics.staticClauses >>= 1;

	  // Sort "_watches" again, since we might have some new binary clauses (mandatory for "saveStatus()" & "restoreStatus()".

	  // Sort "candidates" by decreasing "Literals Blocks Distance" values.
	  if (_setting->simplifyActivity == SIMP_LBD)
		{
		  std::sort(_candidates.begin(), _candidates.end(), antom::CompareLBD);
		  assert(_candidates.empty() || _candidates[0]->Lbd() >= _candidates.back()->Lbd());
		}
	  else 
		{
		  assert(_setting->simplifyActivity == SIMP_CONF);
		  std::sort(_candidates.begin(), _candidates.end(), antom::CompareActivity);
		  assert(_candidates.empty() || _candidates[0]->Activity() <= _candidates.back()->Activity());
		}
    
	  // Resize "candidates".
	  _candidates.resize(_candidates.size() >> 1); 

	  size_t unlearnedclauses = _candidates.size();

	  // Mark all candidates as "to be deleted" by setting 
	  // the corresponding "Literals Blocks Distance" values to 0.
	  
	  while (!_candidates.empty())
		{ 
		  assert(_candidates.back()->IsLearned());
		  _candidates.back()->MarkForDeletion(); 
		  _candidates.pop_back(); 
		}
	  	  
	  // Remove all clauses marked as "to be deleted" from 
	  // the clause database and update "_watches". 
	  for (uint32_t c = 0; c < size; ++c)
		{
		  Clause& clause(_ca[_clauseDatabase[c]]);
		  assert(clause.Activity()>0.0);
		  // Clause marked for deletion?
		  if (clause.ToDelete())
			{	      
			  _ca.Free(_clauseDatabase[c]);
			  --size;
			  _clauseDatabase[c] = _clauseDatabase[size];
			  _clauseDatabase.pop_back();
			  --c;
			}
		  else
			{
			  // Update "_watches".
			  AttachClause(_clauseDatabase[c], false);

			  if (clause.size() == 3 )
				{
				  ++_statistics.currentTernaryClauses; 
				}
			  else
				{
				  ++_statistics.currentNaryClauses;
				}
		  
			  // Update "_learntClauses" if necessary.
			  if (clause.IsLearned())
				{ 
				  ++_statistics.learnedClauses;
				  _statistics.learnedLength += clause.size();
				}
			  else
				{
				  ++_statistics.staticClauses;
				  _statistics.staticLength += clause.size();
				}
			}
		}

	  for (uint32_t v = 1; v <= _variables; ++v)
		{
		  uint32_t lit(v << 1); 
		  for (uint32_t r = 0; r < 2; ++r)
			{ 
			  lit ^= 1;
			  std::sort(_watches[lit].begin(), _watches[lit].end(), WatchedSorter()); 
			}
		}


	  uint32_t deletedclauses = 
		(tmpbinary-(_statistics.currentBinaryClauses>>1)) + 
		(tmpternary-_statistics.currentTernaryClauses) +
		(tmpnary-_statistics.currentNaryClauses);

	  if ( _setting->verbosity > 2)
		{
		  std::cout << "c-------------------------" << std::endl;
		  std::cout << "c simplify #" << _statistics.simplifications << ", " << deletedclauses << " clauses deleted, " << unlearnedclauses << " by unlearning " << std::endl;
		  if ( _setting->verbosity > 3 )
			{
			  std::cout << "c #used Variables..: " << _statistics.usedVariables << " (" << tmpusedvars << ")" << std::endl
						<< "c #binary clauses..: " << (_statistics.currentBinaryClauses>>1) << " (" << tmpbinary << ")" << std::endl
						<< "c #ternary clauses.: " << _statistics.currentTernaryClauses << " (" << tmpternary << ")" << std::endl
						<< "c #nary clauses....: " << _statistics.currentNaryClauses << " (" << tmpnary << ")" << std::endl
						<< "c #total clauses...: " << ((_statistics.currentBinaryClauses>>1)+_statistics.currentTernaryClauses+_statistics.currentNaryClauses) << " (" << (tmpbinary+tmpternary+tmpnary) << ")" << std::endl
					<< "c-------------------------" << std::endl;
			}
		}

	  _newUnits = 0;
	  CheckGarbage();
	}

	bool CheckMaxIndex(uint32_t lastIndex) const;

	void PrintDatabase(bool printWatches = true) const;

	// Writes current clauses into "db"
	void GetClauseDatabase(std::vector< std::vector< uint32_t > >& db) const
	{
	  std::vector< uint32_t > clause(2);

	  for ( uint32_t v = 0; v != _variables; ++v )
		{
		  const std::vector<Watcher>& watches = _watches[v<<1];
		  size_t size(watches.size());

		  for( size_t i = 0; i != size; ++i )
			{
			  if( !watches[i].IsBinary() )
				{ continue; }

			  uint32_t lit = watches[i].GetSecondLit();
			  if( v > ( lit>>1 ) )
				{
				  clause[0] = (v<<1);
				  clause[1] = lit;
				  db.push_back(clause);
				}
			}

		  const std::vector<Watcher>& watchesneg = _watches[(v<<1)^1];
		  size = watchesneg.size();
	
		  for( size_t i = 0; i != size; ++i )
			{
			  if( !watchesneg[i].IsBinary() )
				{ continue; }

			  uint32_t lit = watchesneg[i].GetSecondLit();
		
			  if( v > ( lit>>1 ) )
				{
				  clause[0] = (v<<1)^1;
				  clause[1] = lit;
				  db.push_back(clause);
				}
			}
		}
	  
	  for( uint32_t d = 0; d != _clauseDatabase.size(); ++d )
		{
		  const Clause& lits(_ca[_clauseDatabase[d]]);
		  assert( !lits.ToDelete() );

		  uint32_t size(lits.size());
		  clause.clear();
		  for( uint32_t pos = 0; pos < size; ++pos )
			{
			  clause.push_back( lits[pos] ); 
			}
		  db.push_back(clause);
		}
	}

	// Dumps cnf into std::cout
	void DumpCNF(bool printAssignment = false) const;

	// Saves the current status of the SAT solving core. 
	// The following variables/vectors are not stored:
	// -- "_control",
	// -- "_id",
	// -- "_progress", 
	// -- "_binaryConflictingClause". 
	void SaveStatus(void);

	// Deletes the status of the SAT solving core saved before by "saveStatus()". 
	void DeleteStatus(void);

	// Restores the status of the SAT solving core saved before by "saveStatus()".
	void RestoreStatus(void);

	// Resets the SAT solving core. The following status flags/variables remain untouched:
	// -- The SAT solving threads unique ID number: "_id".
	// -- The pointer to the "Control" object: "_control".
	// -- The number of variables for which memory has been reserved: "_variables".
	void Reset(void);

	void InstanceReset(void);

	friend class AntomBase;
	friend class Antom;
	friend class UPLA;
	friend class BCE;
	friend class BVA;
  protected:

	// Updates all data structures depending on the number of variables to be able to handle "var" variables.
	void UpdateDataStructures(uint32_t var)
	{
	  //std::cout << __FUNCTION__ << " " << var << std::endl;
	  if( var <= _capacity )
		{ 
		  _variables = var;
		  return; 
		}

	  // Update "_level".
      _level.resize(var + 1, 0); 

      // Update "_activity".
      _activity.resize(var + 1, 0.0);

      // Update "_polarity".
      _polarity.resize(var + 1, _setting->initialPolarity);

      // Update "_forcing".
      _forcing.resize(var + 1); 

      // Update "_model".
      _model.resize(var + 1, 0);

      // Update "_dl2ds".
      _dl2ds.resize(var + 1, 1); 

	  // Update "_modeDSforvariables"
	  _modeDSForVariables.resize(var+1,CACHEANDTOGGLE);

      // Update "_varGroup".
      _varGroup.resize(var + 1, 0); 

      // Update "_deleted".
      _deleted.resize(var+1, false);

      // Update "_decisionStack".
      _decisionStack.resize(var + 1, 0); 

	  _conflictClause.resize(var+1, 0);

	  _seen.resize(var+1, false);
	  _touched.resize(var+1, false);

	  // Update "_varOrder".
	  for( uint32_t v = 0; v != _noOfVarGroups; ++v )
		{ 
		  _varOrder[v]->resize(_variables); 
		}
	  
      // Initialization.
      uint32_t max((var << 1) + 2);

      // Update "_assignment".
      _assignment.resize(max, false);
	  
      // Update "_watches".
      _watches.resize(max);

      // Update "_variables" and "_capacity".
      _variables = var; 
	  _capacity = var;
    }

	// Update variable order data structures for "group" var order groups
	void UpdateVarorder(uint32_t group)
	{
	  assert( group >= _noOfVarGroups);
	  _varOrder.resize(group+1);
	  _modeDS.resize(group+1,CACHEANDTOGGLE);
	  _decVarSign.resize(group+1, false);

	  
	  for( uint32_t v = _noOfVarGroups; v <= group; ++v )
		{ 
		  VarHeap<helper::DescendingOrder, double >* varheap = new VarHeap<helper::DescendingOrder, double >( helper::DescendingOrder<double>(_activity) );
		  _varOrder[v] = varheap;
		}	  
	  
	  _noOfVarGroups = group+1;
	}

    // Helper function for pushing assumptions onto the decision stack.
    // In case of a conflict, FALSE will be returned, otherwise the 
    // return value is TRUE. 
    bool SetAssumptions(const std::vector<uint32_t>& assumptions)
    {
	  _failedAssumptions.clear();
      // What about the empty clause?
      if (_emptyClause)
		{ assert( _control->GetExtendedResult() == ANTOM_UNSAT); return false; }

      // Initialization.
      size_t aSize(assumptions.size());

      // Set all assignments specified by "assumptions".
      for (size_t a = 0; a < aSize; ++a)
		{
		  // Get the next assignment.
		  uint32_t lit(assumptions[a]);

		  // Variable indices are assumed to be greater 0.
		  assert(lit > 1); 

		  // Increment "_basicOps".
		  ++_basicOps; 

		  if ((lit >> 1) > _variables )
			{
			  if (_setting->verbosity > 1 )
				{
				  std::cout << "c WARNING: declaration of assumption literal " << helper::Lit(lit) << " which exceeds max index " << _variables << std::endl;
				}
			  SetMaxIndex(lit>>1);
			  // Update "_varOrder".
			  for( uint32_t v = 0; v != _noOfVarGroups; ++v )
				{ 
				  _varOrder[v]->resize(_variables); 
				}
			  continue;
			}

		  // "lit" already incorrectly assigned?
		  if (_assignment[lit ^ 1])
			{
			  // Collect failed assumptions 
			  uint32_t i(0);
			  do
				{
				  _failedAssumptions.push_back(assumptions[i]);
				}
			  while( assumptions[i++] != lit );

			  _control->SetExtendedResult(ANTOM_UNSAT_WITH_ASSUMPTION);
			  if (_decisionLevel != 0)
				{ Backtrack(0); }
			  return false; 
			}

		  // "lit" currently unassigned?
		  if (!_assignment[lit])
			{ 
			  // Push "lit" as a "fake decision" onto the decision stack.
			  AddDecision(lit); 
	    
			  // What about the effects of this implication?
			  if ( !Deduce() )
				{
				  // Collect failed assumptions 
				  uint32_t i(0);
				  do
					{
					  _failedAssumptions.push_back(assumptions[i]);
					}
				  while( assumptions[i++] != lit );

				  _control->SetExtendedResult(ANTOM_UNSAT_WITH_ASSUMPTION);
				  Backtrack(0); 
				  return false; 
				}
			}
		}

      // Everything went fine.
      return true;
    }

	std::vector<uint32_t> GetFailedAssumptions(void) const
	  {
		return _failedAssumptions;
	  }
	

    // Adds a decision to the decision stack.
    void AddDecision(uint32_t lit)
    {
	  //std::cout << __FUNCTION__ << " " << helper::Lit(lit) << " on dl: " << (_decisionLevel+1) << " activity: " << _activity[lit>>1] << std::endl;
	
      // Increment "_decisionLevel".
      ++_decisionLevel; 

      // The variable corresponding to "lit" has to be undefined.
      assert(!_assignment[lit] && !_assignment[lit ^ 1]); 
    
      // Update "_assignment".
      _assignment[lit] = true;

	  // Update "_level".
      _level[lit >> 1] = _decisionLevel;

      // Update "_forcing".
      _forcing[lit >> 1].ClearReason();

      // Update "_dl2ds".
      _dl2ds[_decisionLevel] = _dsEndIndex;  

      // Push "lit" onto the decision stack.    
      _decisionStack[_dsEndIndex++] = lit; 
    }


    // Adds an unconditionally forced implication to the decision stack.
    void AddImplication(uint32_t lit)
    {
	  //std::cout << __FUNCTION__ << " " << helper::Lit(lit) << " on dl: " << _decisionLevel << std::endl;

	  // Update "_implications".
	  ++_statistics.implications;

	  // The variable corresponding to "lit" has to be undefined.
	  assert(!_assignment[lit] && !_assignment[lit ^ 1]);

	  // Update "_assignment".
	  _assignment[lit] = true;

	  // Update "_level".
	  _level[lit >> 1] = _decisionLevel;

	  // Update "_forcing".
	  _forcing[lit >> 1] = Reason();

	  // Push "lit" onto the decision stack.
	  _decisionStack[_dsEndIndex++] = lit;
    }

    // Adds an implication forced by a non-binary clause to the decision stack.
    void AddImplication(uint32_t lit, const CRef& reason)
    {
	  //std::cout << __FUNCTION__ << " " << helper::Lit(lit) << " on dl: " << _decisionLevel << " forced by: ";
	  //_ca[reason].Print();
		  
      // Update "_implications".
      ++_statistics.implications; 
	
      // The variable corresponding to "lit" has to be undefined.
      assert(!_assignment[lit] && !_assignment[lit ^ 1]); 

      // Update "_assignment".
      _assignment[lit] = true;

      // Update "_level".
      _level[lit >> 1] = _decisionLevel;
      
      // Update "_forcing".
      _forcing[lit >> 1] = Reason(reason, false);

      // Push "lit" onto the decision stack.    
      _decisionStack[_dsEndIndex++] = lit; 
    }

    // Adds an implication forced by a binary clause to the decision stack.
	void AddBinaryImplication(uint32_t lit, uint32_t forcingLit)
    {
	  //std::cout << __FUNCTION__ << " " << helper::Lit(lit) << " on dl: " << _decisionLevel << " forced by: " << helper::Lit(forcingLit) << std::endl;

	
	  // Update "_implications".
	  ++_statistics.implications; 
	
      // The variable corresponding to "lit" has to be undefined.
      assert(!_assignment[lit] && !_assignment[lit ^ 1]); 

      // Update "_assignment".
      _assignment[lit] = true;

	  assert( !_assignment[lit^1] );
   
      // Update "_level".
      _level[lit >> 1] = _decisionLevel;
      
      // Update "_forcing".
      _forcing[lit >> 1] = Reason(forcingLit, true);

      // Push "lit" onto the decision stack.    
      _decisionStack[_dsEndIndex++] = lit; 
    }

    // Adds an implication forced by a binary clause to the decision stack.
    void AddTernaryImplication(uint32_t lit, uint32_t forcingLit1, uint32_t forcingLit2)
    {    
	  //std::cout << __FUNCTION__ << " " << helper::Lit(lit) << " forcing: " << helper::Lit(forcingLit1) << " " << helper::Lit(forcingLit2) << " on dl: " << _decisionLevel << std::endl;

      // Update "_implications".
      ++_statistics.implications; 

      // The variable corresponding to "lit" has to be undefined.
      assert(!_assignment[lit] && !_assignment[lit ^ 1]); 
    
      // Update "_assignment".
      _assignment[lit] = true;

      // Update "_level".
      _level[lit >> 1] = _decisionLevel;
      
      // Update "_forcing".
      _forcing[lit >> 1] = Reason(forcingLit1, forcingLit2);

      // Push "lit" onto the decision stack.    
      _decisionStack[_dsEndIndex++] = lit; 
    }

    // Selects the next decision variable and adds it to the decision stack using "addDecision()".
    // Returns TRUE as long as decision variables (still unassigned variables) are available.
    bool Decide(void)
    {
	  // Do we have another variable within "_varOrder[v]"?
	  for( uint32_t v = 0; v != _noOfVarGroups; ++v )
		{
		  while (!_varOrder[v]->empty())
			{
			  // Increment "_basicOps".
			  ++_basicOps; 

			  // Get the literal corresponding to the variable with the maximum activity
			  // and take the variable's last polarity into account. If "_decVarSign" 
			  // is set to TRUE, the cached truth value gets inverted.
			  uint32_t var(_varOrder[v]->top());
			  uint32_t lit(((var << 1) ^ _polarity[var]) ^ _decVarSign[v]);
			  // "lit" currently undefined?
			  if (!_assignment[lit] && !_assignment[lit ^ 1])
				{
				  // Increment "_decisions".
				  ++_statistics.decisions;

				  // Modify "lit" 
				  if (_modeDS[v] == ALWAYSFALSE)
					{ 
					  lit = (var << 1) ^ 1; 
					}
				  else
					{ 
					  if (_modeDS[v] == ALWAYSTRUE)
						{ 
						  lit = var << 1; 
						}
					}

				  // Do the same for variable specific values
				  if (_modeDSForVariables[var] == ALWAYSFALSE)
					{ 
					  lit = (var << 1) ^ 1; 
					}
				  else if (_modeDSForVariables[var] == ALWAYSTRUE)
					{
					  lit = var << 1; 
					}
		
				  // Push "lit" as a decision onto the decision stack.
				  AddDecision(lit); 
	      
				  // Everything went fine.
				  return true;
				}
			}
		}
 
      // At this point, there's no variable left to serve as a decision variable.
      return false;
    }

	bool DeduceBinary(Watcher& watcher, uint32_t& lit)
	{
	  uint32_t otherlit(watcher.GetSecondLit());
	  // At this point, we have either an implication or a conflict.
	  if (_assignment[otherlit ^ 1])
		{
		  // conflict case
		  _conflictLiteral = otherlit;
		  _conflict = Reason(lit, true);
		  return false;  
		}
	  else if ( !_assignment[otherlit] )
		{
		  // Push the second literal as an implication forced
		  // by a binary clause onto the decision stack.
		  AddBinaryImplication(otherlit, lit);
		}
	  return true;
	}


	bool DeduceTernary(Watcher& watcher, uint32_t& lit)
	{
	  uint32_t secondlit(watcher.GetSecondLit());
	  if( _assignment[secondlit])
		{ return true; }

	  uint32_t thirdlit(watcher.GetThirdLit());
	  if( _assignment[thirdlit])
		{ 
		  // Set second lit to assigned lit
		  watcher.SetSecondLit(thirdlit);
		  watcher.SetThirdLit(secondlit);
		  return true;
		}

	  // secondlit unassigned, thirdlit assigned to false
	  // -> propagte secondlit
	  if( !_assignment[secondlit^1] && _assignment[thirdlit^1] )
		{
		  AddTernaryImplication(secondlit, lit, thirdlit); 
		}

	  // thirdlit unassigned, secondlit assigned to false
	  // -> propagte thirdlit
	  else if( _assignment[secondlit^1] && !_assignment[thirdlit^1] )
		{
		  AddTernaryImplication(thirdlit, lit, secondlit); 

		  // Set second lit to assigned lit
		  watcher.SetSecondLit(thirdlit);
		  watcher.SetThirdLit(secondlit);
		}
	  // both lits assigned to false
	  // -> conflict
	  else if( _assignment[secondlit^1] && _assignment[thirdlit^1] )
		{

		  // conflict case
		  assert( ( _level[secondlit>>1] == _level[lit>>1] ) || ( _level[thirdlit>>1] == _level[lit>>1] ) );
		  // mark secondlit as conflict
		  _conflictLiteral = secondlit;
		  // Set lit and thirdlit as reason for conflict
		  _conflict = Reason(lit, thirdlit);
		  return false;  
		}

	  // add this point secondlit and thirdlit are unassigned
	  // -> continue with bcp
	  return true;
	}

	bool DeduceNary(uint32_t& wpos, uint32_t& size, uint32_t& lit)
	{
	  Watcher& watch(_watches[lit][wpos]);

	  // Get the literals of the next clause to be checked.
	  Clause& clause(_ca[watch.GetClause()]);

	  // Ensure that watchlit is at second position of clause
	  // and clause[0] contains the potentially "open" literal
	  // important for conflict analysis
	  if (clause[0] == lit)
		{
		  clause[0] = clause[1];
		  clause[1] = lit;
		}
	  assert(clause[1] == lit);
	  
	  // Current clause satisfied by "wl"?
	  if (_assignment[clause[0]])
		{
		  // Update the blocking literal.
		  watch.SetBlockingLiteral(clause[0]);
		}
	  else
		{
		  // Initialization ("pos" is set to 2, since there is no need 
		  // to check the current watched literals again).
		  uint32_t clausesize(_ca[watch.GetClause()].size());
			  
		  // Try to find a new watched literal within the remainder of "cLits".
		  // This loop assumes that "_assignment[0]" and "_assignment[1]" 
		  // are both set to FALSE.

		  uint32_t pos(2);
		  //assert((clause[pos]>>1) <= _variables);
		  //assert((clausesize+watch.GetClause()) <= _ca.GetSize());
		  while( pos < clausesize && _assignment[clause[pos]^1] )
			{
			  ++pos;
			}
							  
		  if(pos != clausesize)
			{
			  // Found new watch
			  // Get the literal checked last.
			  uint32_t newlit(clause[pos]); 
			      
			  // swap new watch literal with old one:
			  clause[1]   = newlit;
			  clause[pos] = lit;

			  // Update "_watches" for both "lit" and "newlit".
			  _watches[newlit].push_back(watch);
			  --wpos; --size;
			  return true;
			}
		
		  // At this point, we either have an implication or a conflict.
		      	      
		  // Conflict?
		  if (_assignment[clause[0] ^ 1])
			{
			  _conflict = Reason(watch.GetClause(), false);
			  /* std::cout << " conflict of literal: " << helper::Lit(clause[0]) << std::endl << "clause: "; */
			  /* clause.Print(); */
			  /* std::cout << "reason: "; */
			  /* _conflict.Print(_ca); */
			  _conflictLiteral = clause[0];
			  return false;
			}
		  else
			{
			  // Push "wl" as an implication onto the decision stack.
			  AddImplication(clause[0], watch.GetClause());
			}
		}
	  return true;
	}

	// Performs "Boolean Constraint Propagation". In case of a conflict, the reason
	// will be updated and "true" returned, otherwise "false" will be returned
    bool Deduce()  
	{
	  // Increment "_bcps".
      ++_statistics.bcps;

	  _conflict.ClearReason();

	  uint32_t decision(_decisionStack[_dl2ds[_decisionLevel]] ^ 1);
	  bool lhbr(_setting->lhbr);
      if (_dsImplIndex != _dl2ds[_decisionLevel])
		{ lhbr = false; } 

	  bool okay(true);

      // Process all assignments not checked so far. 
      while (_conflict.NoReason() && (_dsImplIndex != _dsEndIndex) )
		{
		  // Get the next assignment.
		  uint32_t wlit(_decisionStack[_dsImplIndex++] ^ 1);

		  // "wlit" has to be assigned, since it's part of the decision stack.
		  assert(_assignment[wlit] || _assignment[wlit ^ 1]);

		  // Initialization.
		  std::vector<Watcher>& watches(_watches[wlit]);

		  uint32_t size((uint32_t)watches.size()); 

		  // Increment "_basicOps".
		  _basicOps += size; 
 
		  uint32_t v(0);
		  uint32_t w(0);
		  
		  // Process all clauses in which "wlit" currently serves as a watched literal.
		  for (; w < size; ++w)
			{
			  watches[w] = watches[v++];

			  if( _assignment[watches[w].GetSecondLit()] )
				{ 
				  continue; 
				}

			  if( watches[w].IsBinary() )
				{
				  if ( !DeduceBinary(watches[w], wlit) )
					{ 
					  assert(_conflict.IsBinary());
					  ++w;
					  okay = false;
					  break;
					}
				  else
					{ 
					  if ( lhbr && (decision != wlit) )
						{
						  // Increment "_lhbr".
						  ++_statistics.lhbr; 
						  _statistics.currentBinaryClauses += 2;
						  _statistics.learnedLength += 2;
						  _statistics.globalLBD += 2; 
						  _statistics.localLBD += 2;
						  ++_statistics.localConflicts;
						  ++_statistics.learnedBinary;
						  
						  uint32_t firstlit(watches[w].GetSecondLit());
						  // Consistency check.
						  assert(firstlit != wlit); 

						  IncreaseActivity(firstlit>>1);
						  IncreaseActivity(decision>>1);

						  // Create a new binary clause by updating "_watches".
						  _watches[firstlit].push_back(Watcher(decision, true));
						  _watches[decision].push_back(Watcher(firstlit, true));
						}
					  assert(_conflict.NoReason()); 
					  continue; 
					}
				}
	  
			  if( watches[w].IsTernary() )
				{

				  lhbr=false;
				  if( !DeduceTernary(watches[w], wlit) )
					{ 
					  assert(_conflict.IsTernary());
					  ++w;
					  okay = false;
					  break;
					}
				  else
					{ 
					  assert(_conflict.NoReason()); 
					  continue; 
					}
				}
			  if( watches[w].IsClause() )
				{
				  lhbr=false;
				  if ( !DeduceNary(w, size, wlit) )
					{ 
					  assert(_conflict.IsClause());
					  ++w;
					  okay = false;
					  break;
					}
				  else
					{ 
					  assert(_conflict.NoReason()); 
					  continue; 
					}
				}
			}
		  
		  // Finally, copy remaining watches
		  // Note, we never remove binaries and ternaries from watchlist
		  // -> We only need to copy remaining watches after "deduceNary"
		  while ( w != size )
			{
			  watches[w++] = watches[v++];
			}
		  watches.resize(size);
		}

	  return okay;
	}

    // Performs conflict analysis according to the 1UIP scheme and backtracks afterwards.
    void Analyze(void)
    {
      // Increment "_conflicts".
      ++_statistics.totalConflicts;
	  ++_statistics.localConflicts;

      // Everything OK wrt. "_conflict"?
      assert(_conflict.HasReason());

      // Update "_incVarActivity".
      _incVarActivity *= _setting->decayFactor;
 
      // Initialization.
      uint32_t elements(0);
      uint32_t pos(_dsEndIndex - 1); 
      uint32_t uip(0);
	  uint32_t size(1);

	  _seenToClear.clear();
	  _touchedToClear.clear();

	  // Set Reason
	  //std::cout << "init reason: conf " << helper::Lit(_conflictLiteral) << " clause: ";
	  //_conflict.Print(_ca);
	  ReasonComplete cc(_conflict, _conflictLiteral, _ca);

      // Perform conflict analysis according to the 1UIP scheme.
      do
		{
		  // Increment "_basicOps".
		  ++_basicOps;

		  // Analyze the literals of the current conflicting clause.
		  uint32_t reasonsize(cc.size());
		  for( uint32_t cpos = 0; cpos < reasonsize; ++cpos )
			{
			  // Get the index of "clit".
			  uint32_t index(cc[cpos] >> 1);

			  // "clit" not checked so far and not assigned on decision level 0? 
			  // Checking the decision level on which a particular variable has been 
			  // assigned, requires that assumptions are stored on decision levels 
			  // greater 0. Otherwise we run into problems in the incremental mode. 
			  if (!_seen[index] && _level[index] > 0)
				{
				  // Update "seen".
				  _seen[index] = true;
				  _seenToClear.push_back(index);

				   // Increase the activity of the variable corresponding to "clit".
				  IncreaseActivity(index);

				  // Has "clit" been assigned on the current decision level?
				  if (_level[index] == _decisionLevel)
					{ 
					  ++elements; 
					}
				  else
					{ 
					  _conflictClause[size++] = cc[cpos];
					}
				}
			}

		  // Determine the next clause to be processed. 
		  while (!_seen[_decisionStack[pos] >> 1]) 
			{ --pos; }

		  uip = _decisionStack[pos];

		  // Update the status variables.

		  Reason& r(_forcing[uip>>1]);

		  if (r.IsClause())
			{
			  IncreaseActivity(_ca[r.GetClause()]);
			}

		  cc = ReasonComplete(r, uip, _ca);
		  // Either we reached the last element or there has to be a reason for "uip"
		  assert(elements == 1 || !(_forcing[uip >> 1].NoReason()));
		  --pos;
		  --elements;

		}
      while (elements > 0);

      // Flip the sign of the UIP.
      uip = uip ^ 1;

	  _conflictClause[0] = uip;
  
	  // Update lbd value and backtrack level
      uint32_t lv(0);
      uint32_t at(1);
      uint32_t lbd(1);

	  // Perform a simple conflict clause minimization step. See also "Towards Understanding 
	  // and Harnessing the Potential of Clause Learning" by Beame, Kautz, and Sabharwal.

	  // abstraction level allowing to skip unnecessary operations within "LitRedundant"
	  uint32_t abstractLevel = 0;
	  if (_setting->ccMinimization == DEEP)
		{
		  for (uint32_t l = 1; l < size; ++l)
			{
			  abstractLevel |= (1<<(_level[_conflictClause[l]>>1]&31));
			}
		}

	  for (uint32_t l = 1; l < size; ++l)
		{
		  // Get the next literal of "conflictClause".
		  uint32_t lit(_conflictClause[l]); 
		  uint32_t index(lit>>1);

		  // Do we talk about an implication?
		  if ( (_forcing[index].HasReason() ) )
			{
			  bool redundant = false;
			  if (_setting->ccMinimization == BASIC)
				{
				  ReasonComplete rc(_forcing[index], lit^1, _ca);
				  uint32_t rcpos = 1;
				  uint32_t rcsize = rc.size();
				  for( ; rcpos != rcsize; ++rcpos )
					{
					  uint32_t curIndex( rc[rcpos]>>1 );
					  if ( !_seen[curIndex] && _level[curIndex] > 0 )
						{ break; }
					}
				  redundant = (rcpos == rcsize);
				}
			  else if(_setting->ccMinimization == DEEP)
				{
				  redundant = LitRedundant(lit, abstractLevel);
				}
			  
			  // Is it safe to remove "lit"?
			  if( redundant )
				{
				  ++_statistics.minimizedLiterals;
				  --size; 
				  _conflictClause[l] = _conflictClause[size];
				  --l; 
				  continue;
				}
			}
		}

	  //clear seen
	  for(uint32_t toClear : _seenToClear)
		{
		  _seen[toClear] = false;
		}

	  // Determine lbd and forcing lit
	  for (uint32_t l = 1; l < size; ++l)
		{
		  uint32_t index = _conflictClause[l]>>1;
		  // Do we have to update the backtrack level?
		  if (_level[index] > lv)
			{ 
			  lv = _level[index]; 
			  at = l;
			}
		  
		  // Do we have to update the "Literals Blocks Distance"?
		  if (!_touched[_level[index]])
			{ 
			  _touched[_level[index]] = true;
			  _touchedToClear.push_back(_level[index]);
			  ++lbd; 
			}
		}

	  for(uint32_t toClear : _touchedToClear)
		{
		  _touched[toClear] = false;
		}
 
	  // Do we have an unit clause?
      if (size == 1)
		{
		  // Update some of our statistics related variables. 
		  ++_statistics.globalLBD; 
		  ++_statistics.localLBD; 
		  ++_statistics.learnedLength;
		  ++_statistics.totalLearntUnitClauses;
		  _statistics.DLclearedCA += _decisionLevel; 
		  _statistics.VarsUnassignedCA += _dsEndIndex - _dl2ds[1]; 

		  // Backtrack to decision level 0.
		  Backtrack(0);

		  //std::cout << "new unit: " << Lit(uip) << std::endl;

		  // Add the UIP as a "fake" implication to the decision stack.
		  AddImplication(uip);
		  ++_newUnits;
		  return;
		}
      // Do we have a binary clause?
      else if (size == 2)
		{	  
		  // Update some of our statistics related variables. 
		  _statistics.currentBinaryClauses += 2;
		  _statistics.globalLBD += 2; 
		  _statistics.localLBD += 2; 
		  _statistics.learnedLength +=2;
		  _statistics.totalLearnedLength +=2;
		  _statistics.DLclearedCA += _decisionLevel - _level[_conflictClause[1] >> 1]; 
		  _statistics.VarsUnassignedCA += _dsEndIndex - _dl2ds[_level[_conflictClause[1] >> 1] + 1]; 
		  ++_statistics.learnedBinary;
		  ++_statistics.totalLearntBinaryClauses; 

		  _watches[_conflictClause[1]].push_back( Watcher(uip, true) );
		  _watches[uip].push_back( Watcher(_conflictClause[1], true) );

		  // Backtrack to the decision level on which "conflictClause[1]" has been assigned.
		  Backtrack(_level[_conflictClause[1] >> 1]); 

		  // Push the UIP as an implication forced by a binary clause onto the decision stack.
		  AddBinaryImplication(uip, _conflictClause[1]);

		  return; 
		}
	  else if ( size == 3 )
		{
		  ++_statistics.currentTernaryClauses;
		  ++_statistics.totalLearntTernaryClauses;
		}
	  else
		{
		  ++_statistics.currentNaryClauses;
		}

      // Update "_learnedClauses".
      ++_statistics.learnedClauses;

      // Swap the literal on the second position of "cptr" for the one on 
      // position "at" to have the "correct" pair of watched literals.
	  std::swap(_conflictClause[1],_conflictClause[at]);

      // Update statistics
      _statistics.globalLBD += lbd; 
      _statistics.localLBD += lbd; 
      _statistics.learnedLength += size; 
      _statistics.totalLearnedLength += size; 
      _statistics.DLclearedCA += _decisionLevel - lv; 
	  _statistics.VarsUnassignedCA += _dsEndIndex - _dl2ds[lv + 1]; 

      // Backtrack to decision level "lv".
      Backtrack(lv);

      // Create a new clause.
	  CRef cr(_ca.Alloc(_conflictClause, lbd, _activityInc, size));

      // Add a new clause to the clause database.
	  AttachClause(cr);

      // Add the UIP as an implication to the decision stack.
      AddImplication(uip, cr);
    }

	// checks if lit is redundant wrt the current conflict clause
	// recursively checks whether all literals on implication chain of "lit"  are already included in the current conflicht clause "_conflichtClause"
	// -> we can remove "lit" in this case
	bool LitRedundant(uint32_t lit, uint32_t abstractLevel)
	{
	  _analyzeStack.clear();
	  _analyzeStack.push_back(lit);
	  size_t currentStacksize = _seenToClear.size();

	  while ( !_analyzeStack.empty() )
		{
		  uint32_t curLit = _analyzeStack.back();
		  _analyzeStack.pop_back();

		  assert(_forcing[curLit>>1].HasReason());
		  
		  ReasonComplete rc(_forcing[curLit>>1], curLit^1, _ca);

		  for(size_t rcpos = 1; rcpos != rc.size(); ++rcpos )
			{
			  uint32_t rcIndex( rc[rcpos]>>1 );
			  if ( !_seen[rcIndex] && _level[rcIndex] > 0 )
				{
				  uint32_t alevel =	(1<<(_level[rcIndex]&31));
					if (_forcing[rcIndex].HasReason() && ((alevel & abstractLevel) != 0 ))
					{
					  _seen[rcIndex] = true;
					  _seenToClear.push_back(rcIndex);
					  _analyzeStack.push_back(rc[rcpos]);
					}
				  else
					{
					  for (size_t j = currentStacksize; j < _seenToClear.size(); ++j)
						{
						  _seen[_seenToClear[j]] = false;
						}
					  _seenToClear.resize(currentStacksize);
					  return false;
					}
				}
			}
		}
	  return true;
	}

		
    // Backtracks to decision level "bL". 
    void Backtrack(uint32_t bL)
    {
      // "bL" has to be less than the current decision level.
      assert(bL < _decisionLevel);

      // Update "_avgDL".
      _statistics.DL += _decisionLevel; 

      // Initialization.
      uint32_t stopper(_dl2ds[bL + 1]);

      // Increment "_basicOps".
      _basicOps += _dsEndIndex - stopper; 
   
      // Undo all variable assignments on decision levels higher than "bL".
      // This loop assumes that we have a dummy assignment at the first 
      // position of the decision stack, which has been assigned on decision level 0. 
      do
		{
		  // Decrement "_dsEndIndex".
		  --_dsEndIndex;

		  // Get the next assignment.
		  uint32_t lit(_decisionStack[_dsEndIndex]); 

		  // Cache the polarity of the current assignment.
		  _polarity[lit >> 1] = (lit & 1); 
	  
		  // Undo the current assignment.
		  _assignment[lit] = false;

		  // Clear forcing clause (needed for incremental mode)
		  _forcing[lit>>1].ClearReason();
	  
		  // Update "_varOrder".
		  if (!_varOrder[_varGroup[lit >> 1]]->inHeap(lit >> 1))
			{
			  _varOrder[_varGroup[lit >> 1]]->insert(lit >> 1); 
			}
		}
      while (stopper != _dsEndIndex); 

      // Consistency check.
      assert(stopper == _dsEndIndex); 

      // Update "_decisionLevel".
      _decisionLevel = bL; 
   
      // Update "_dsImplIndex".
      _dsImplIndex = _dsEndIndex;
    }


    // Increases the activity of variable "var".
    void IncreaseActivity(uint32_t var, bool checkheap = true)
    {
      // "var" has to be less or equal "_variables".
      assert(var <= _variables); 

      // Increase "var's" activity.
      _activity[var] += _incVarActivity; 

      // Update the position of "var" within the "_varOrder".
      if (checkheap && _varOrder[_varGroup[var]]->inHeap(var))
		{ _varOrder[_varGroup[var]]->update(var); } 

      // Do we have to "normalize" the variables' activities?
      if (_activity[var] > 1e100)
		{
		  // Increment "_basicOps".
		  _basicOps += _variables; 

		  // Divide all activities by 1e100.
		  for (uint32_t v = 1; v <= _variables; ++v)
			{ _activity[v] *= 1e-100; }
	  
		  // Update "_incVarActivity".
		  _incVarActivity *= 1e-100;
		}
    }

	void IncreaseActivity(Clause& clause)
	{
	  clause.IncreaseActivity(_activityInc);
	  if ( clause.Activity() > 1e20 )
		{
		  // Rescale:
		  for (uint32_t i = 0; i < _clauseDatabase.size(); ++i)
			{
			  _ca[_clauseDatabase[i]].ScaleActivity();
			}
		  //_activityInc *= static_cast<float>(1e-20);
		  _activityInc *= 1e-20;  
		}
	}

	void EstimateProgress(void) 
	{
	  _statistics.progress = 0.0; 
	  double div(1.0); 
	  for (uint32_t l = 0; l <= _decisionLevel; ++l)
		{
		  uint32_t end(_dsEndIndex); 
		  if (l < _decisionLevel)
			{ end = _dl2ds[l + 1]; }
		  if (l != 0)
			{ div *= _variables + 1; }
		  _statistics.progress += static_cast<double>(end - _dl2ds[l]) / div; 
		}
	}
	      
	// Lifting methods
	std::vector< uint32_t > AnalyzeLifting(uint32_t* cc );

	// lifts a solution, given by the variables in "assumptions"
	// liftingmodi:
	// 0001 : conflict driven
	// 0010 : brute force driven
	// 0011 : combine 1+2
	// 01xx : add all assumptions at once (only for conflict driven approach)
	// note: needs a negated property
	// the formula with negated property must be UNSAT, with non-negated property SAT
	// Returns lifted solution
	std::vector< uint32_t > SolveLifting(std::vector<uint32_t>& assumptions, uint32_t liftingmode, uint32_t sortmode);

	// Returns:
	// "-1" if variable is assigned to false
	// " 1" if variable is assigned to true
	// " 0" if variable is unassigned
	int32_t IsAssigned(uint32_t lit) const;

	void PrintClause(Clause* clause, bool assignment) const;
	bool CheckClauses() const;

	void PrintBinaryList(uint32_t lit) const;
	void PrintWatchList(uint32_t lit) const;

	// Build watchers 
	void AttachClause(CRef cr, bool addToDatabase = true)
	{
	  //std::cout << __func__ << " cref: " << cr << std::endl;
	  Clause& clause(_ca[cr]);

	  // Special storage for ternary clauses 
	  if( _setting->useTernary && clause.size() == 3 )
		{

		  _watches[clause[0]].push_back(Watcher(clause[1], clause[2],TERNARY));
		  _watches[clause[1]].push_back(Watcher(clause[0], clause[2],TERNARY));
		  _watches[clause[2]].push_back(Watcher(clause[0], clause[1],TERNARY));
		}
	  else
		{
		  _watches[clause[0]].push_back(Watcher(cr,clause[1],NNARY));
		  _watches[clause[1]].push_back(Watcher(cr,clause[0],NNARY));
		}

	  if( addToDatabase )
		{
		  // Finally, add the new clause to the clause database.
		  _clauseDatabase.push_back(cr);
		}
	}

	// Garbage collection
	
	// relocate all clauses to a new allocator "to"
	void RelocAll(ClauseAllocator& to)
	{
	  // All watchers:
	  //
	  for (uint32_t v = 1; v <= _variables; ++v)
		{
		  for (uint32_t s = 0; s < 2; ++s)
			{
			  uint32_t literal = (v<<1)+s;
			  //std::cout << "relocating watchers of " << Lit(literal) << std::endl;
			  std::vector< Watcher >& watches(_watches[literal]);
			  size_t wsize(watches.size());
			  for (size_t j = 0; j < wsize; ++j)
				{
				  if( watches[j].IsClause() )
					{
					  CRef ref(watches[j].GetClause());
					  _ca.Reloc(ref, to);
					  watches[j].SetClause(ref);
					}
				}
			}
        }

	  // All reasons:
	  //
	  for (uint32_t i = 1; i < _dsEndIndex; ++i)
		{
		  uint32_t index(_decisionStack[i]>>1);
		  Reason& reason(_forcing[index]);
		  if ( reason.IsClause() )
			{
			  CRef ref(reason.GetClause());
			  _ca.Reloc(ref, to);
			  reason.SetClause(ref);
			}
		}
	  
	  // All clauses:
	  //
	  size_t size(_clauseDatabase.size());
	  for (uint32_t i = 0; i < size; ++i)
	  {
        _ca.Reloc(_clauseDatabase[i], to);
	  }
	}
	
	void CheckGarbage(void)
	{
	  //std::cout << __FUNCTION__ << " wasted: " << _ca.Wasted() << " size: " << _ca.size() << std::endl;
	  // collect garbage if 20% of the memory is wasted
	  if( _ca.Wasted() > (_ca.size() * 0.2) )
		{
		  // Initialize the next region to a size corresponding to the estimated utilization degree. This
		  // is not precise but should avoid some unnecessary reallocations for the new region:
		  ClauseAllocator to(_ca.size() - _ca.Wasted());

		  // reallocate all clauses to new allocator
		  RelocAll(to);
		  if (_setting->verbosity >= 2)
			{
			  std::cout << "c Garbage collection: " << ((_ca.size()*ClauseAllocator::Unit_Size)/(1024*1024)) << " MB => " << ((to.size()*ClauseAllocator::Unit_Size)/(1024*1024)) << " MB\n";
			}

		  // move "to" the core allocator "_ca" 
		  to.MoveTo(_ca);
		}
	}

    // Copy constructor.
    Core (const Core&) = default;

    // Assignment operator.
    Core& operator = (const Core&) = default;

    // Each core has its own unique ID.
    uint32_t _id; 

    // Each core has a pointer to the "Control" object.
    Control* _control;

	Settings* _setting;

	// The clause allocator
    ClauseAllocator _ca;

    // A flag indicating whether we have been able to deduce the empty clause.
    // In particular important for the incremental mode to distinguish between
    // "unsatisfiable" and "unsatisfiable under assumptions". 
    bool _emptyClause;

	// Start and end pointer for unit clause synchronization
	uint32_t _startPtr;
	uint32_t _endPtr;

    // The current number of variables.
    uint32_t _variables; 

	// The reserved size of the variable related data structures
	uint32_t _capacity;

	Statistics _statistics;

    // The (partial) variable assignment.
    std::vector<unsigned char> _assignment;

    // The decision level on which a particular literal has been assigned.
    std::vector<uint32_t> _level; 

    // For each variable we store its activity.
    std::vector<double> _activity;  

	// Clause activity increment value;
	float _activityInc;

    // For each variable we cache its "last" polarity.
    std::vector<unsigned char> _polarity;

    // For each variable we store whether it is an implication or a decision.
    std::vector<Reason> _forcing;
	// Temporal storage for last conflict reason
	Reason _conflict;
	uint32_t _newUnits;

    // For each literal we store a vector to represent in which binary & non-binary
    // clauses a particular literal currently serves as a watched literal.
    std::vector<std::vector< Watcher > > _watches;

    // Some variable activity related variables.
	std::vector< VarHeap<helper::DescendingOrder, double >* > _varOrder; 
    std::vector<uint32_t> _varGroup; 
	uint32_t _noOfVarGroups;
    double _incVarActivity;

    // Some decision stack related variables.
    uint32_t _decisionLevel; 
    std::vector<uint32_t> _decisionStack;
	BoundedQueue<size_t> _stackQueue;
    std::vector<uint32_t> _dl2ds;
    uint32_t _dsEndIndex;
    uint32_t _dsImplIndex;

    // Some decision strategy related variables. 
	std::vector< DecisionStrategy > _modeDS; 
	std::vector< DecisionStrategy > _modeDSForVariables; 

    // Some conflict analysis related variables.
	std::vector<uint32_t> _analyzeStack;
	std::vector<uint32_t> _seenToClear;
	std::vector<uint32_t> _touchedToClear;
	std::vector<unsigned char> _seen;
	std::vector<unsigned char> _touched;
	uint32_t _conflictLiteral; 
	std::vector<uint32_t> _conflictClause;
	std::vector<uint32_t> _failedAssumptions;
	uint32_t _learnedClausesLimit;

    // Some restart related variables.
	std::vector< unsigned char > _decVarSign;
	std::vector<Clause*> _candidates;
	   
    // Some more status variables.
    uint64_t _basicOps; 

    // The satisfying variable assignment (in case the CNF formula is 
    // satisfiable and has been solved beforehand; otherwise "_model" 
    // might contain invalid data).
    std::vector<uint32_t> _model;
	
    // The clause database.
    std::vector<CRef> _clauseDatabase;

	SolverState _solverState; 

	friend class Preprocessor;
	Preprocessor* _preprocessor;

	// Deleted variables
	std::vector< unsigned char > _deleted;
  };
}

#endif
