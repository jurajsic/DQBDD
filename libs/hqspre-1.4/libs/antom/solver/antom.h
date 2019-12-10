
/********************************************************************************************
antom.h -- Copyright (c) 2014-2016, Sven Reimer

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

#ifndef ANTOM_ANTOM_H_
#define ANTOM_ANTOM_H_

// Include standard headers.
#include "antombase.h"
#include "settings.h"

namespace antom
{
  enum SolverType
  {
	ANTOMSOLVER,
	CRYPTOMINISATSOLVER
  };
	
  class SoftClause;
  class Sorter;
  class Cascade;
  class MultipleCascade;
  struct TimeVariables;

  // The "MaxAntom" class.
  class Antom : public AntomBase
  {
  public:

	// Constructor. 
    Antom(void);
   
	// Destructor.
	~Antom(void); 

	void SetControl(void);

	/* MaxSAT interface */

	// Returns the last used index. Does not necessary meet with "variables()"
	uint32_t LastIndex(void) const;

	// Add a clause to the soft clause database
    bool AddSoftClause(std::vector<uint32_t>& clause, uint64_t weight = 1);
    // Add weighted SoftClauses to _softClauseVector
    bool AddWeightedSoftClause(std::vector<uint32_t>& clause, uint64_t weight = 1);
    // Add whole _softclause Vector to clause database
    bool AddSoftClauseVector();

	// Set predefined bounds
	void SetMaxBounds(uint32_t low, uint32_t high);
	void SetLowerBound(uint32_t low);
	void SetHigherBound(uint32_t high);

	// Set a target for optimal solution.
	// Declares minimal number of unsatisfied clauses
	// If target is reached, stop maxSAT-computation. 
	void SetOptTarget(int32_t target);
	// States whether the target should be reached as precise as possible
	void SetPreciseTarget(bool val);

	// Search modes for maxSAT
	void SetSearchMode(SearchMode val);

	// Partial modes for maxSAT
	void SetPartialMode(PartialMode val);

	// Solve partial maxSAT Problem
	uint32_t MaxSolve(int64_t& optimum);
	uint32_t MaxSolve(const std::vector< uint32_t >& externalAssumptions, int64_t& optimum);

    // Solve wmSAT with warners
    uint32_t MaxSolveWarnersWeightedPartial(int64_t& optimum);
    uint32_t MaxSolveWarnersWeightedPartial(const std::vector< uint32_t >& externalAssumptions, int64_t& optimum);

//        void genWarners0(std::vector<long long int>& weights, std::vector<Lit>& blockings,
//                long long int max,long long int k, int comp, Solver& S,
//                std::vector<Lit>& lits, std::vector<Lit>& linkingVar)
        void GenWarners0(   int comp,
                            std::vector<uint32_t>& lits,
                            std::vector<uint32_t>& linkingVar);
        void GenWarners(    std::vector< SoftClause* > softCl,
                            uint64_t max,
                            int logk,
                            int comp,
                            uint32_t zero,
                            std::vector<uint32_t>& lits,
                            std::vector<uint32_t>& linkingVar);
        void GenWarnersHalf(    uint32_t& a,
                                uint32_t& b,
                                uint32_t& carry,
                                uint32_t& sum,
                                int comp,
                                std::vector<uint32_t>& lits);
        void GenWarnersFull(    uint32_t& a,
                                uint32_t& b,
                                uint32_t& c,
                                uint32_t& carry,
                                uint32_t& sum,
                                int comp,
                                std::vector<uint32_t>& lits);


	// Solve weighted partial maxSAT Problem
    uint32_t MaxSolveWeightedPartial(int64_t& optimum);
    uint32_t MaxSolveWeightedPartial(const std::vector< uint32_t >& externalAssumptions, int64_t& optimum);

	/* Some maxSAT related interface functions */
	void SetMaxWidth(uint32_t width);
	void SetSplittedWidth(uint32_t width);
	void SetMaxParts(uint32_t parts);
	void SetSortSoftClauses(uint32_t val);
	void CalcSplitWidth(void);
	void SetSkip(bool val);
	void SetCSC(uint32_t val);
	void SetRelaxationLits(bool val);
	void SetIncompleteMode(bool val);
	void SetHorizontalWidth(uint32_t val);
	void SetGridMode(uint32_t val);
	void SetNetworktype(SorterType val);
	void SetMaxInprocessing(bool val);

	// val = 0 --> switch off incremental mode for maxSAT (default)
	// val = 1 --> use maxSAT part of SAT incrementally (soft clauses are deleted after each call).
	// val = 2 --> use maxSAT part of SAT incrementally (soft clauses are kept after each call).
	void SetIncrementalMode(uint32_t val);

	// De-/activates constant detection with SAT in preprocessing
	void SetSatConst(uint32_t val);
 
	// Resets MaxSAT related data
	void Reset(void);
	void InstanceReset(void);

	void SetStrategyforSoftVariables(bool pos);

	void PrintStatistics() const;
	
    // Weighted MaxSat Stuff
    StructureInfo AnalyzeandConvertStructure();
    void SetEncode01Mode(bool val, bool val2);
    void SetSumOfSoftWeights(uint64_t val);
    void SetBaseMode(uint32_t val);
    void SetPartitionStrategy(PartitionStrategy val);
    void SetHeuristic(uint32_t val);
    void SetPercentOff(uint32_t val);
    void SetPercentOffreinsert(bool val);
    void SetEqualWeight(uint32_t val);
    void SetEncodeStrategy(EncodeStrategy val);
    void SetCreateGraphFile(std::string val);
    void SetFeatureTest(bool val);
    void SetCascadeDivider(uint32_t val);
    void SetMaxBucketSize(uint32_t val);
    void SetNOfCascades(uint32_t val);
    void SetSolveTareCascadeOnlyByTares(bool val);
    void SetMultipleCascadeDivideStrategy(MultipleCascadeDivideStrategy val);
    void SetSepHiWeight(bool val);
    void SetOnlyByTares(bool val);

    void SetCorrectInputFormat(bool val);
	void SetDecStratMode(uint32_t val);
	void SetBypassesMode(uint32_t val);
    void SetMoreThanTwoWeights(bool val);
    void SetTopWeight(uint64_t val);
    void SetMinWeight(uint64_t val);
    void SetMaxWeight(uint64_t val);
    void SetSolveAtFirst(bool val);
    void SetWeightPlusOne(bool val);
    void SetInterimResult(InterimResult val);

    friend class Cascade;
    friend class MultipleCascade;
    friend class Bucket;
    friend class Sorter;
	friend class Counter;
	friend class Parser;
	
  private:

	// Copy constructor.
    Antom (const Antom&) = default;

    // Assignment operator.
    Antom& operator = (const Antom&) = default;

	void DataReset(void);

    // Solves the current (partial) MaxSAT formula. Returns SAT/UNSAT and modifies "optimum", representing the minimum number of 
    // unsatisfied clauses.
    uint32_t MaxSolveLocal(int64_t& optimum);
    uint32_t MaxSolveLocal(const std::vector<uint32_t>& externalAssumptions, int64_t& optimum);

	// Check for constants in the candidates list by applying the SAT solver
	// Attention! Use with care, process might be costly
	// If quickcheck is set, only assumptions are deduced, instead of solving the instance
	// -> less powerful, but much faster
	bool FindConstantsWithSat(std::vector< uint32_t >& candidates, bool quickcheck);
	
	uint32_t SatBasedSearch(Sorter* sorter, uint32_t& pos);
	uint32_t UnsatBasedSearch(Sorter* sorter, uint32_t& pos);
	uint32_t BinaryBasedSearch(Sorter* sorter, uint32_t& pos);

	uint32_t GetPreciseTarget(Sorter* sorter, uint32_t& pos, int32_t gap);

	bool AddBound(uint32_t lit);

	uint32_t SetLowerPartialParts(std::vector< uint32_t >& localassumptions) const;
	int32_t UpdateSorterOptimum(void);

    uint64_t CountSatisfiedSoftClauses(Sorter* sorter, const std::vector<uint32_t>& model);
    uint64_t CountSatisfiedSoftClauses(std::vector<SoftClause *> softclauses, const std::vector<uint32_t> &model);

	//void mergePartTrigger ( Sorter& trigger1, Sorter& trigger2 );
	void CheckAllConflictingSoftclauses(void);

	// Store already found assignments of relaxation lits to "assumptions"
	void CollectRelaxationLits(std::vector< uint32_t >& assumptions);

	Sorter* CreateNextSorter(uint32_t depth = 0, uint32_t targetdepth = 1);
	Sorter* MergeNextSorter(void);
	Sorter* GetNextSorter(void);

	// Deactives soft clasues, deletes sorter clauses and variables
    void InvalidateSoftClauses(void);

	// Cout the model of the buckets and which position has the first set bit.
	void DumpBucketModel(const std::vector<uint32_t>& model);

    // Called from Cascade and Bucket - they only calculate their own optimum.
    uint64_t CalculateOverallOptimum(uint64_t SatWeight, bool countAgain);

    // Add to each TopBucket _base-1 many Tare variables, exept the last one.
    void CheckAllWeightedConflictingSoftclauses(void);

	void SetTareBounds(void);

	void ExitTimeout(void);

	void SetDecisionStrategiesForMaxSAT(void);

#ifndef NDEBUG
	void CheckGates(void) const;
#endif

    std::vector< SoftClause* > _softClauses;

	std::vector< std::vector< Sorter* > > _sorterTree;

	std::vector< uint32_t > _bestModel;

	uint32_t _maxSorterDepth;

	uint32_t _maxsatResult;

    // needed for updating optimum value from Weighted Mode - Cascade and Bucket
    int64_t* _optimum;

	bool _lastpartialsorter;
	// some statistics
	uint32_t _horizontalbypasses;
	uint32_t _verticalbypasses;
	uint32_t _comparator;
	uint32_t _skipped;

	uint32_t _triggervar;
	// last index of original problem (without sorter network variables)
	uint32_t _lastIndex;

	// Bounds for maximization mode
	uint32_t _low;
	uint32_t _high;

	//Weighted MaxSAT stuff
    TimeVariables* _timeVariables;
    Cascade* _mainCascade;
    MultipleCascade* _mainMultipleCascade;
    void FillBuckets();
    bool EncodeBuckets();
    int64_t _greatestCommonDivisor;
    enum StructureInformation
    {
        NOERROR,
        CONVERTTOMAXSATFORMULA,
        ISSATFORMULA,
        DIVIDEWEIGHTSBYDIVISOR,
        ISMAXSATFORMULA,
        HASONLYTWOWEIGHTS,      // not yet used
        NOTOPWEIGHT,            // not yet used
        SOFTWEIGHTSTOOBIG       // not yet used
    } _formulaStructure;
    StructureInfo AnalyzeStructure();
    //StructureInformation AnalyzeStructure();
    void ConvertFormulaToMaxSAT(uint64_t maxWeight);
    void DivideAllSoftClausesByFactor(uint64_t factor);
    uint64_t GreatestCommonDivisor(uint64_t a, uint64_t b);

    uint32_t _currentBucketForTare;
	uint64_t _satWeight;
    bool _moreThanTwoWeights;
    uint64_t _topWeight;
    uint64_t _minWeight;
    uint64_t _maxWeight;
    uint64_t _sumOfSoftWeights;
	
    // if last solver call returned ANTOMUNKNOWN - this variable is set to true
    bool _resultUnknown;

    // maybe later on as local vars again - just for testing!!
    uint32_t _clausesBefore;
    uint32_t _binaryClausesBefore;
    uint32_t _ternaryClausesBefore;
    uint32_t _addedClauses;
    uint32_t _addedBinaryClauses;
    uint32_t _addedTernaryClauses;

    //for coding clause estimation in standard mode
    double _binaryTopClauses;
    uint64_t _ternaryTopClauses;
    uint64_t _binaryBottomClauses;
    uint64_t _ternaryBottomClauses;

	// partial status variables
    uint64_t _highestBucketMultiplicator;
    uint64_t _currentBucketMultiplicator;
	// the value calculated highestBucketMultiplicator * bucketInfo[1];
	uint64_t _estimatedSatWeight;
	uint64_t _diffToSatWeight;
	std::vector<uint32_t> _collectedAssumptions;
  };
}

#endif
