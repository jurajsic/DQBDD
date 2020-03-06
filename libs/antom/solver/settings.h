/********************************************************************************************
settings.hpp -- Copyright (c) 2017, Sven Reimer

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

#ifndef ANTOM_SETTINGS_H_
#define ANTOM_SETTINGS_H_

#include <stdint.h>
#include <string>
#include <iostream>
#include <cassert>

namespace antom 
{
  
  enum Application
  {
	SAT,
	MAXSAT,
	WEIGHTEDMAXSAT,
  };
  
  enum SimplifyStrategy 
  {
	ANTOM,
	MINISAT,
  };

  enum SimplifyActivity
  {
	SIMP_LBD,
	SIMP_CONF,
  };

  enum CCMinimization
  {
	NOMIN,
	BASIC,
	DEEP,
  };

  enum PreproType 
  {
	NOPREPRO,
	PREPROCESS,
	INPROCESS,
	INCREMENTAL,
  };

  enum RestartStrategy 
  {
	LUBY,
	GLUCOSE,
  };

  enum SearchMode 
  {
	SATBASED,
	UNSATBASED,
	BINARYBASED,
  };

  enum PartialMode 
  {
	NONE,
	DEPTHFIRST,
	BREADTHFIRST,
  };

  enum EncodeStrategy
  {
    ENCODEALL,
    ENCODEONLYIFNEEDED
  };

  enum PartitionStrategy
  {
    NOPARTITION,
    GROUPBYWEIGHTADDATLAST,
    GROUPBYWEIGHT,
    GROUPBYBIGGESTREPEATINGENTRY,
  };

  enum MultipleCascadeSolvingState
  {
    SINGLECASCADE,
    MAINCASCADE,
    TARECASCADE,
    TARETARES
  };

  enum MultipleCascadeDivideStrategy
  {
	SOLVEINNORMALCASCADEMODE,
	SORTEDNUMBEROFSOFTCLAUSES,
	RANDOMNUMBEROFSOFTCLAUSES,
	SOFTCLAUSESINORDER,
	SORTEDGOODDIVIDEPOINTS,
  };

  enum InterimResult
  {
	NOINTERIMRESULT,
	CUTATTOP,
	CUTATBOTTOM,
	CUTBOTH,
  };

  enum StructureInfo
  {
    ISSAT,
    CONVERTTOMAXSAT,
    ISMAXSAT,
    ISWEIGHTEDMAXSAT,
    DIVIDEWEIGHTSBYDIVISOR,
  };

  enum GateType 
  {
	ANDGATE,
	NANDGATE,
	ORGATE,
	NORGATE,
	XORGATE,
	XNORGATE,
	HALFANDGATE,
	HALFORGATE,
	INPUT,
  };

  std::string inline GateTypeToString(GateType type)
	{
	  switch(type)
		{
		case ANDGATE :
		  return "AND";
		case NANDGATE :
		  return "NAND";
		case ORGATE :
		  return "OR";
		case NORGATE :
		  return "NOR";
		case XORGATE :
		  return "XOR";
		case XNORGATE :
		  return "XNOR";
		case HALFANDGATE :
		  return "HALFAND";
		case HALFORGATE :
		  return "HALFOR";
		case INPUT :
		  return "INPUTAND";
		default:
		  std::cerr << "ERROR! Unknown gate type" << std::endl;
		  assert(false);
		  return "";
		}
	  return "";
	}


  enum DecisionStrategy 
  {
	CACHEANDTOGGLE,
	CACHE,
	ALWAYSFALSE,
	ALWAYSTRUE,
  };

  enum SorterType
  {
	BITONIC,
	ODDEVEN,
	TOTALIZER,
    WARNERS,
  };
  
  // Container for all solver settings
  struct Settings
  {
	
  Settings()
	{
	  Reset(); //initialize with default settings
	}
	
	Settings(const Settings&)            = default;
	Settings(Settings&&)                 = default;
	~Settings()                          = default;
	Settings& operator=(const Settings&) = default;
	Settings& operator=(Settings&&)      = default;
	
	void SetInstance(Settings* Settings);
	void Reset(void);
	void ResetGeneral(void);
	void ResetCore(void);
	void ResetPrepro(void);
	void ResetAntom(void);

	void Print() const;

	// General Settings
	uint32_t threads;
	uint32_t verbosity;
	double cpuLimit;
	double memLimit;
	Application application;
	
	// Core Settings
	double decayFactor;
	uint32_t lubyShift;
	RestartStrategy restartStrategy;
	bool restartBlocking;
	CCMinimization ccMinimization;
	bool initialPolarity;
	SimplifyStrategy simplifyStrategy;
	SimplifyActivity simplifyActivity;
	DecisionStrategy decisionStrategy; 
	bool useTernary;
	bool lhbr;
	PreproType doPreprocessing;
	bool doInprocessing;
	
	// Prepro Settings
	uint32_t maxLoops;
	bool doSubsumption;
	bool doUpla;
	bool doVarElimination;
	bool doBce;
	bool doHte;
	bool doHse;
	bool doBva;
	bool bvaTwoLitDiff;
	bool doVivification;
	int32_t varIncrease;
	uint32_t satconst;

	// Antombase Settings

	// MaxSAT Settings
	// Incremental MaxSAT mode
	// 0 - no incremental mode (default)
	// 1 - "BMC"ish: soft clauses are deleted after each call
	// 2 - soft clauses are kept after each call
	uint32_t incrementalMode;
	SearchMode searchMode;
	PartialMode partialMode;
	uint32_t decStrat;
	uint32_t splittedWidth;
	uint32_t maxWidth;
	uint32_t maxParts;
	uint32_t horizontalWidth;
	uint32_t bypassGrid;
	uint32_t sortSoftClauses;
	uint32_t setCSC;
	bool doSkipping;
	bool setRelaxationLits;
	bool incompleteMode;
	bool maxInprocess;
	uint32_t satConst;
	SorterType networkType;
	bool preciseTarget;
	int32_t targetOpt;

	// Weighted MaxSAT Settings
	bool analyze;
	bool encode01;
	bool lastPos1;
	uint32_t base;
	uint32_t groupHeuristic;
    uint32_t percentOff;
    bool percentOffReinsert;
    uint32_t equalWeight;
    PartitionStrategy partitionStrategy;
    bool solveAtFirst;
    EncodeStrategy encodeStrategy;
    std::string createGraphFile;
    bool onlyByTares;

    // multiple cascade mode
    MultipleCascadeDivideStrategy mcDivideStrategy;
    uint32_t cascadeDivider;
    uint32_t maxBucketSize;
    uint32_t nOfCasc;
    bool tareCascadeOnlyByTares;
    bool sepHiWeight;
    bool weightPlusOne;

	// interim results
    InterimResult interimResult;

    bool featureTest;
  };
}

#endif
