/********************************************************************************************
parser.cpp -- Copyright (c) 2013-2017, Tobias Schubert, Sven Reimer, Tobias Paxian

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

#include "parser.h"
#include "antom.h"
//#include "solverproxy.h"

#include <iostream>
#include <sstream>
#include <string>
#include <limits>

namespace antom 
{
  Parser::Parser(void) :
	antomSettings(nullptr),
	storeResult(false),
	source(),
	filename(""),
	resultFile(""),
	mode(0),
	verify(false),
	assumptions(),
	instance(1),
	doIncrementalSolve(false),
	doSolve(false),
	doReset(false)
  {
  }

  uint32_t Parser::ParseCommandline(antom::Antom& myAntom, int argc, char** argv)
  {
	antomSettings = myAntom._antomSetting;
	// Get the additional command line parameters.
	for (uint32_t c = 1; c < ((uint32_t) argc - 1); ++c)
	  {
		// Initialization.
		std::string argument(argv[c]); 
		bool matched(false);      

		// Basic options
		// Result file
		if (argument.compare(0, 14, "--result-file=") == 0)
		  { std::stringstream ss(argument.substr(14, argument.length())); ss >> resultFile; storeResult = true; matched = true; continue; }
		// Operating mode
		else if (argument == "--mode=0")
		  { mode = 0; matched = true; continue;}
		else if (argument == "--mode=1")
		  { mode = 1; matched = true; continue; }
		else if (argument == "--mode=10")
		  { mode = 10; matched = true; continue; }
		else if (argument == "--mode=20")
		  { mode = 20; matched = true; continue; }
		else if (argument == "--mode=30")
		  { mode = 30; matched = true; continue; }
		else if (argument == "--mode=31")
		  { mode = 31; matched = true; continue; }
                else if (argument == "--mode=99")
		  { mode = 99; matched = true; continue; }
		// Number of threads
		else if (argument == "--threads=1")
		  { antomSettings->threads = 1; matched = true; continue; }
		else if (argument == "--threads=2")
		  { antomSettings->threads = 2; matched = true; continue; }
		else if (argument == "--threads=3")
		  { antomSettings->threads = 3; matched = true; continue; }
		else if (argument == "--threads=4")
		  { antomSettings->threads = 4; matched = true; continue; }
		else if (argument == "--threads=5")
		  { antomSettings->threads = 5; matched = true; continue; }
		else if (argument == "--threads=6")
		  { antomSettings->threads = 6; matched = true; continue; }
		else if (argument == "--threads=7")
		  { antomSettings->threads = 7; matched = true; continue; }
		else if (argument == "--threads=8")
		  { antomSettings->threads = 8; matched = true; continue; }
		else if (argument == "--verbose" || argument == "--v")
		  { ++antomSettings->verbosity; matched = true; continue; }
		else if (argument == "--verbose=0")
		  { antomSettings->verbosity = 0; matched = true; continue; }
		// Verify model 
		else if (argument == "--verify=true")
		  { verify = true; matched = true; continue; }
		else if (argument == "--verify=false")
		  { verify = false; matched = true; continue; }
		// Backend solver
		else if (argument == "--solver=antom")
		  { solver = ANTOMSOLVER; matched = true; continue; }
		else if (argument == "--solver=cryptominisat")
		  { solver = CRYPTOMINISATSOLVER; matched = true; continue; }
		// CPU Time limit
		else if (argument.compare(0,11,"--cpuLimit=") == 0)
		  { std::stringstream sCPU(argument.substr(11,argument.length())); sCPU >> antomSettings->cpuLimit; matched = true; continue; }
		// Memory Limit
		else if (argument.compare(0,11,"--memLimit=") == 0)
		  { std::stringstream sMem(argument.substr(11,argument.length())); sMem >> antomSettings->memLimit; matched = true; continue; }
		// Core Options
		// Restart strategy
		else if (argument == "--restart=luby")
		  { antomSettings->restartStrategy = LUBY; matched = true; continue; }
		else if (argument == "--restart=glucose")
		  { antomSettings->restartStrategy = GLUCOSE; matched = true; continue; }
		else if (argument == "--restartBlock=true")
		  { antomSettings->restartBlocking = true; matched = true; continue; }
		else if (argument == "--restartBlock=false")
		  { antomSettings->restartBlocking = false; matched = true; continue; }
		// What about the unit factor?
		else if (argument.compare(0, 13, "--unitFactor=") == 0)
		  { std::stringstream ss(argument.substr(13, argument.length())); ss >> antomSettings->lubyShift; matched = true; continue; }      // What about the decay factor?
		else if (argument.compare(0, 14, "--decayFactor=") == 0)
		  { std::stringstream ss(argument.substr(14, argument.length())); ss >> antomSettings->decayFactor; matched = true; continue; }
		// Simplifaction strategy
		else if (argument == "--simplify=antom")
		  { antomSettings->simplifyStrategy = ANTOM; matched = true; continue; }
		else if (argument == "--simplify=minisat")
		  { antomSettings->simplifyStrategy = MINISAT; matched = true; continue; }
		// Conflict clause minimization
		else if (argument == "--ccmin=none")
		  { antomSettings->ccMinimization = NOMIN; matched = true; continue; }
		else if (argument == "--ccmin=basic")
		  { antomSettings->ccMinimization = BASIC; matched = true; continue; }
		else if (argument == "--ccmin=deep")
		  { antomSettings->ccMinimization = DEEP; matched = true; continue; }
		// Base for simplifaction activity 
		else if (argument == "--activity=lbd")
		  { antomSettings->simplifyActivity = SIMP_LBD; matched = true; continue; }
		else if (argument == "--activity=conf")
		  { antomSettings->simplifyActivity = SIMP_CONF; matched = true; continue; }
		// Decision strategy?
		else if (argument == "--decision=0")
		  { antomSettings->decisionStrategy = CACHEANDTOGGLE; matched = true; continue; }
		else if (argument == "--decision=1")
		  { antomSettings->decisionStrategy = CACHE; matched = true; continue; }
		else if (argument == "--decision=2")
		  { antomSettings->decisionStrategy = ALWAYSFALSE; matched = true; continue; }
		else if (argument == "--decision=3")
		  { antomSettings->decisionStrategy = ALWAYSTRUE; matched = true; continue; }
		// Initial polarity of variables
		else if (argument == "--initPol=false")
		  { antomSettings->initialPolarity = false; matched = true; continue; }
		else if (argument == "--initPol=true")
		  { antomSettings->initialPolarity = true; matched = true; continue; }
		// Special treatment for ternary clauses?
		else if (argument == "--ternary=true")
		  { antomSettings->useTernary = true; matched = true; continue; }
		else if (argument == "--ternary=false")
		  { antomSettings->useTernary = false; matched = true; continue; }
		// Lazy Hyper Binary Resolution
		else if (argument == "--lhbr=true") 
		  { antomSettings->lhbr = true; matched = true; continue; }
		else if (argument == "--lhbr=false")
		  { antomSettings->lhbr = false; matched = true; continue; }
		// Preprocessor
		// Perfoming preprocessing
		else if (argument == "--prepro=false")
		  { antomSettings->doPreprocessing = NOPREPRO; matched = true; continue; }
		else if (argument == "--prepro=true")
		  { antomSettings->doPreprocessing = PREPROCESS; matched = true; continue; }
		else if (argument == "--prepro=incremental")
		  { antomSettings->doPreprocessing = INCREMENTAL; matched = true; continue; }
		// Perfoming inprocessing
		else if (argument == "--inpro=true")
		  { antomSettings->doInprocessing = true; matched = true; continue; }
		else if (argument == "--inpro=false")
		  { antomSettings->doInprocessing = false; matched = true; continue; }
		// inprocessing during MaxSAT
		else if (argument == "--maxInpro=true")
		  { antomSettings->maxInprocess = true; matched = true; continue; }
		else if (argument == "--maxInpro=false")
		  { antomSettings->maxInprocess = false; matched = true; continue; }
		// Maximum preprocessor loops
		else if (argument.compare(0, 11, "--maxLoops=") == 0)
		  { std::stringstream ss(argument.substr(11, argument.length())); ss >> antomSettings->maxLoops; matched = true; continue; }
		// UPLA
		else if (argument == "--upla=true")
		  { antomSettings->doUpla = true; matched = true; continue; }
		else if (argument == "--upla=false")
		  { antomSettings->doUpla = false; matched = true; continue; }
		// Subsumption checks
		else if (argument == "--subsumption=true")
		  { antomSettings->doSubsumption = true; matched = true; continue; }
		else if (argument == "--subsumption=false")
		  { antomSettings->doSubsumption = false; matched = true; continue; }
		// Variable elimination
		else if (argument == "--varElim=true")
		  { antomSettings->doVarElimination = true; matched = true; continue; }
		else if (argument == "--varElim=false")
		  { antomSettings->doVarElimination = false; matched = true; continue; }
		// Blocked Clause Elimination
		else if (argument == "--bce=true")
		  { antomSettings->doBce = true; matched = true; continue; }
		else if (argument == "--bce=false")
		  { antomSettings->doBce = false; matched = true; continue; }
		// Hidden Tautology Elimination
		else if (argument == "--hte=true")
		  { antomSettings->doHte = true; matched = true; continue; }
		else if (argument == "--hte=false")
		  { antomSettings->doHte = false; matched = true; continue; }
		// Hidden Subsumption Elimination
		else if (argument == "--hse=true")
		  { antomSettings->doHse = true; matched = true; continue; }
		else if (argument == "--hse=false")
		  { antomSettings->doHse = false; matched = true; continue; }
		// Bounded variable addition
		else if (argument == "--bva=true")
		  { antomSettings->doBva = true; matched = true; continue; }
		else if (argument == "--bva=false")
		  { antomSettings->doBva = false; matched = true; continue; }
		else if (argument == "--2litdiff=true")
		  { antomSettings->bvaTwoLitDiff = true; matched = true; continue; }
		else if (argument == "--2litdiff=false")
		  { antomSettings->bvaTwoLitDiff = false; matched = true; continue; }
		// Vivification
		else if (argument == "--vivify=true")
		  { antomSettings->doVivification = true; matched = true; continue; }
		else if (argument == "--vivify=false")
		  { antomSettings->doVivification = false; matched = true; continue; }
		// Constant checking with SAT
		else if (argument == "--satconst=full")
		  { antomSettings->satconst = 2; matched = true; continue; }
		else if (argument == "--satconst=true")
		  { antomSettings->satconst = 1; matched = true; continue; }
		else if (argument == "--satconst=false")
		  { antomSettings->satconst = 0; matched = true; continue; }
		// MaxSAT options
		// Applying incomplete mode
		else if (argument == "--incomplete=false")
		  { antomSettings->incompleteMode = false; matched = true; continue; }
		else if (argument == "--incomplete=true")
		  { antomSettings->incompleteMode = true; matched = true; continue; }
		// Search mode
		else if (argument == "--search=0")
		  { antomSettings->searchMode = UNSATBASED; matched = true; continue; }
		else if (argument == "--search=1")
		  { antomSettings->searchMode = SATBASED; matched = true; continue; }
		else if (argument == "--search=2")
		  { antomSettings->searchMode = BINARYBASED; matched = true; continue; }
		// Network type
		else if (argument == "--network=0")
		  { antomSettings->networkType = BITONIC; matched = true; continue; }
		else if (argument == "--network=1")
		  { antomSettings->networkType = ODDEVEN; matched = true; continue; }
		else if (argument == "--network=3")
		  { antomSettings->networkType = TOTALIZER; matched = true; continue; }
                else if (argument == "--network=4")
                  { antomSettings->networkType = WARNERS; matched = true; continue; }
		// Decision strategies for MaxSAT
		else if (argument == "--decstrat=0")
		  { antomSettings->decStrat = 0; matched = true; continue; }
		else if (argument == "--decstrat=1")
		  { antomSettings->decStrat = 1; matched = true; continue; }
		else if (argument == "--decstrat=2")
		  { antomSettings->decStrat = 2; matched = true; continue; }
		// Encode 01?
		else if (argument == "--encode01=false")
          { antomSettings->encode01 = false; matched = true; continue; }
		else if (argument == "--encode01=true")
          { antomSettings->encode01 = true; matched = true; continue; }
		else if (argument == "--encode01=lastPos1")
          { antomSettings->encode01 = false; antomSettings->lastPos1 = true; matched = true; continue; }
		else if (argument == "--sortsoft=0")
		  { antomSettings->sortSoftClauses = 0; matched = true; continue; }
		else if (argument == "--sortsoft=1")
		  { antomSettings->sortSoftClauses = 1; matched = true; continue; }
		else if (argument == "--sortsoft=2")
		  { antomSettings->sortSoftClauses = 2; matched = true; continue; }
		else if (argument == "--sortsoft=3")
		  { antomSettings->sortSoftClauses = 3; matched = true; continue; }
		// Bypasses
		else if (argument == "--gridMode=0")
		  { antomSettings->bypassGrid = 0; matched = true; continue; }
		else if (argument == "--gridMode=1")
		  { antomSettings->bypassGrid = 1; matched = true; continue; }
		else if (argument == "--gridMode=2")
		  { antomSettings->bypassGrid = 2; matched = true; continue; }
		else if (argument == "--gridMode=3")
		  { antomSettings->bypassGrid = 3; matched = true; continue; }
		// Width for horizontal bypasses
		else if (argument.compare(0,12,"--bypassWidth=") == 0)
		  { std::stringstream sWidth(argument.substr(12,argument.length())); sWidth >> antomSettings->horizontalWidth; matched = true; continue; }
		// Find conflicting soft clauses
		else if (argument == "--csc=false")
		  { antomSettings->setCSC = 0; matched = true; continue; }
		else if (argument == "--csc=true")
		  { antomSettings->setCSC = 1; matched = true; continue; }
		else if (argument == "--csc=all")
		  { antomSettings->setCSC = 2; matched = true; continue; }
		// Skip comparators
		else if (argument == "--skip=true")
		  { antomSettings->doSkipping = true; matched = true; continue; }
		else if (argument == "--skip=false")
		  { antomSettings->doSkipping = false; matched = true; continue; }
		// Partial MaxSAT Mode
		// Type of partial mode
		else if (argument == "--partial=0")
		  { antomSettings->partialMode = NONE; matched = true; continue; }
		else if (argument == "--partial=1")
		  { antomSettings->partialMode = DEPTHFIRST; matched = true; continue; }
		else if (argument == "--partial=2")
		  { antomSettings->partialMode = BREADTHFIRST; matched = true; continue; }
		// Fixed splitting width for partial mode
		else if (argument.compare(0,13,"--splitWidth=") == 0)
		  { std::stringstream sWidth(argument.substr(13,argument.length())); sWidth >> antomSettings->splittedWidth; matched = true; continue; }
		// Maximal splitting width for partial mode
		else if (argument.compare(0,11,"--maxWidth=") == 0)
		  { std::stringstream sWidth(argument.substr(11,argument.length())); sWidth >> antomSettings->maxWidth; matched = true; continue; }
		// Maximial number of subparts for partial mode
		else if (argument.compare(0,11,"--maxParts=") == 0)
		  { std::stringstream sParts(argument.substr(11,argument.length())); sParts >> antomSettings->maxParts; matched = true; continue; }
		// Use fixed relaxation lits
		else if (argument == "--relax=true")
		  { antomSettings->setRelaxationLits = true; matched = true; continue; }
		else if (argument == "--relax=false")
		  { antomSettings->setRelaxationLits = false; matched = true; continue; }
		// Define target optimum
		else if (argument.compare(0, 12, "--targetOpt=") == 0)
		  { std::stringstream ss(argument.substr(12, argument.length())); ss >> antomSettings->targetOpt; matched = true; continue; }
		// Find closed solution to target value
		else if (argument == "--preciseTarget=true")
		  { antomSettings->preciseTarget = true; matched = true; continue; }
		else if (argument == "--preciseTarget=false")
		  { antomSettings->preciseTarget = false; matched = true; continue; }
		// Base for weighted MaxSAT
		else if (argument.compare(0, 7, "--base=") == 0)
		  { std::stringstream ss(argument.substr(7, argument.length())); ss >> antomSettings->base; matched = true; continue; }
		        // How to partition the weight to reuse in different buckets
        else if (argument == "--partitionStrategy=0")
          { antomSettings->partitionStrategy = NOPARTITION; matched = true; continue; }
		else if (argument == "--partitionStrategy=1")
          { antomSettings->partitionStrategy = GROUPBYWEIGHTADDATLAST; matched = true; continue; }
		else if (argument == "--partitionStrategy=2")
          { antomSettings->partitionStrategy = GROUPBYWEIGHT; matched = true; continue; }
        else if (argument == "--partitionStrategy=3")
          { antomSettings->partitionStrategy = GROUPBYBIGGESTREPEATINGENTRY; matched = true; continue; }
        // If partitionStrategy=GROUPBYBIGGESTREPEATINGENTRY
        // a heuristic has to be chosen to divide the formula.
        else if (argument == "--heuristic=0")
          { antomSettings->groupHeuristic = 0; matched = true; continue; }
        else if (argument == "--heuristic=1")
          { antomSettings->groupHeuristic = 1; matched = true; continue; }
        else if (argument == "--heuristic=2")
          { antomSettings->groupHeuristic = 2; matched = true; continue; }
        else if (argument == "--heuristic=3")
          { antomSettings->groupHeuristic = 3; matched = true; continue; }
        else if (argument == "--heuristic=4")
          { antomSettings->groupHeuristic = 4; matched = true; continue; }
        else if (argument == "--heuristic=5")
          { antomSettings->groupHeuristic = 5; matched = true; continue; }
        else if (argument == "--heuristic=6")
          { antomSettings->groupHeuristic = 6; matched = true; continue; }
        else if (argument == "--heuristic=7")
          { antomSettings->groupHeuristic = 7; matched = true; continue; }
        else if (argument == "--heuristic=8")
          { antomSettings->groupHeuristic = 8; matched = true; continue; }
        // If partitionStrategy=GROUPBYBIGGESTREPEATINGENTRY
        // max difference in sizes to merge SofClauseNodes
        else if (argument.compare(0, 13, "--percentOff=") == 0)
          { std::stringstream ss(argument.substr(13, argument.length())); ss >> antomSettings->percentOff; matched = true;
            antomSettings->percentOff = (antomSettings->percentOff > 100) ? 100 : antomSettings->percentOff; continue; }
        else if (argument == "--percentOffReinsert=true")
          { antomSettings->percentOffReinsert = true; matched = true; continue; }
        else if (argument == "--percentOffReinsert=false")
          { antomSettings->percentOffReinsert = false; matched = true; continue; }
        // If partitionStrategy=GROUPBYBIGGESTREPEATINGENTRY
        // check every n rounds if there are equal weights and merge them together.
        else if (argument.compare(0, 14, "--equalWeight=") == 0)
          { std::stringstream ss(argument.substr(14, argument.length())); ss >> antomSettings->equalWeight; matched = true; continue; }
        // Analyze Structure to recognize type of formula and common divisor!
        // Additionally it converts the formula to that type or divides it by the divisor.
        else if (argument == "--analyze=true")
          { antomSettings->analyze = true; matched = true; continue; }
        else if (argument == "--analyze=false")
          { antomSettings->analyze = false; matched = true; continue; }
        else if (argument == "--solveAtFirst=true")
          { antomSettings->solveAtFirst = true; matched = true; continue; }
        else if (argument == "--solveAtFirst=false")
          { antomSettings->solveAtFirst = false; matched = true; continue; }
        else if (argument == "--encodeStrategy=0")
          { antomSettings->encodeStrategy = ENCODEALL; matched = true; continue; }
        else if (argument == "--encodeStrategy=1")
          { antomSettings->encodeStrategy = ENCODEONLYIFNEEDED; matched = true; continue; }
        else if (argument.compare(0,14,"--createGraph=") == 0)
        { std::stringstream ss(argument.substr(14, argument.length())); ss >> antomSettings->createGraphFile; matched = true;
		  std::string file = (argv[(uint32_t) argc - 1]);
		  unsigned long pos = file.find_last_of("/");
		  file = file.substr(pos + 1);
		  antomSettings->createGraphFile = antomSettings->createGraphFile.append(file);
		  continue;
		}
        // Divide - for multiple cascade mode
        else if (argument.compare(0, 10, "--cascDiv=") == 0)
          { std::stringstream ss(argument.substr(10, argument.length())); ss >> antomSettings->cascadeDivider; matched = true; continue; }
        else if (argument.compare(0, 16, "--maxBucketSize=") == 0)
          { std::stringstream ss(argument.substr(16, argument.length())); ss >> antomSettings->maxBucketSize; matched = true; continue; }
        else if (argument.compare(0, 10, "--nOfCasc=") == 0)
          { std::stringstream ss(argument.substr(10, argument.length())); ss >> antomSettings->nOfCasc; matched = true; continue; }
        else if (argument == "--onlyByTares=true")
          { antomSettings->onlyByTares = true; matched = true; continue; }
        else if (argument == "--onlyByTares=false")
          { antomSettings->onlyByTares = false; matched = true; continue; }
        else if (argument == "--tcOnlyByTares=true")
          { antomSettings->tareCascadeOnlyByTares = true; matched = true; continue; }
        else if (argument == "--tcOnlyByTares=false")
          { antomSettings->tareCascadeOnlyByTares = false; matched = true; continue; }
        else if (argument == "--mcDivideStrategy=0")
          { antomSettings->mcDivideStrategy = SOLVEINNORMALCASCADEMODE; matched = true; continue; }
        else if (argument == "--mcDivideStrategy=1")
          { antomSettings->mcDivideStrategy = SORTEDNUMBEROFSOFTCLAUSES; matched = true; continue; }
        else if (argument == "--mcDivideStrategy=2")
          { antomSettings->mcDivideStrategy = RANDOMNUMBEROFSOFTCLAUSES; matched = true; continue; }
        else if (argument == "--mcDivideStrategy=3")
          { antomSettings->mcDivideStrategy = SOFTCLAUSESINORDER; matched = true; continue; }
        else if (argument == "--mcDivideStrategy=4")
          { antomSettings->mcDivideStrategy = SORTEDGOODDIVIDEPOINTS; matched = true; continue; }
        else if (argument == "--interimResult=0")
          { antomSettings->interimResult = NOINTERIMRESULT; matched = true; continue; }
        else if (argument == "--interimResult=1")
          { antomSettings->interimResult = CUTATTOP; matched = true; continue; }
        else if (argument == "--interimResult=2")
          { antomSettings->interimResult = CUTATBOTTOM; matched = true; continue; }
        else if (argument == "--interimResult=3")
          { antomSettings->interimResult = CUTBOTH; matched = true; continue; }
        else if (argument == "--sepHiWeight=true")
          { antomSettings->sepHiWeight = true; matched = true; continue; }
        else if (argument == "--sepHiWeight=false")
          { antomSettings->sepHiWeight = false; matched = true; continue; }
        else if (argument == "--weightPlus1=true")
          { antomSettings->weightPlusOne = true; matched = true; continue; }
        else if (argument == "--weightPlus1=false")
          { antomSettings->weightPlusOne = false; matched = true; continue; }
        else if (argument == "--featureTest=true")
          { antomSettings->featureTest = true; matched = true; continue; }
        else if (argument == "--featureTest=false")
          { antomSettings->featureTest = false; matched = true; continue; }

  
		// Unknown option?
		if (!matched)
		  {
			// Output. 
			std::cout << "c Unknown option: " << argv[c] << std::endl;
			PrintUsage();
			return ANTOM_UNKNOWN;
		  }
	  }

	filename = argv[(uint32_t) argc - 1];

	// Check inconsistencies
	if( antomSettings->incompleteMode&& (antomSettings->partialMode == 0 ) )
	  {
		std::cout << "c Incomplete mode only supported in combination with a partial mode" << std::endl;
		return ANTOM_UNKNOWN;
	  }
	if( !antomSettings->encode01 && (antomSettings->networkType != TOTALIZER ) )
	  {
		std::cout << "c encoding just 0's only suppted in combination with totalizer network" << std::endl;
		//return ANTOM_UNKNOWN;
	  }
	
	return ANTOM_SAT;
  }

#if 0
  void Parser::SetSettings(antom::Antom& myAntom) const
  {
	// Set antom's various parameters.
	// TODO: this is obsolete
	myAntom.SetMemoryLimit(antomSettings->memLimit);
	myAntom.SetVerbosity(antomSettings->verbosity);
	myAntom.SetCPULimit(antomSettings->cpuLimit);

	myAntom.SetRestartStrategy(antomSettings->restartStrategy); 
	myAntom.SetLuby(antomSettings->lubyShift); 
	myAntom.SetDecayFactor(antomSettings->decayFactor); 
	myAntom.SetSimplifyStrategy(antomSettings->simplifyStrategy);
	myAntom.SetSimplifyActivity(antomSettings->simplifyActivity);
	myAntom.SetDecisionStrategy(antomSettings->decisionStrategy, 0); 
	myAntom.UseTernaryClauses(antomSettings->useTernary);
	myAntom.SetLHBR(antomSettings->lhbr); 

	myAntom.SetPreprocessing(antomSettings->doPreprocessing);
	myAntom.SetInprocessing(antomSettings->doInprocessing);
	myAntom.SetMaxInprocessing(antomSettings->maxInprocess);
	myAntom.SetMaxLoops(antomSettings->maxLoops);
	myAntom.SetUPLA(antomSettings->doUpla);
	myAntom.SetSubsumption(antomSettings->doSubsumption);
	myAntom.SetVarElim(antomSettings->doVarElimination);
	myAntom.SetBCE(antomSettings->doBce);
	myAntom.SetHTE(antomSettings->doHte);
	myAntom.SetHSE(antomSettings->doHse);
	myAntom.SetBVA(antomSettings->doBva);
	myAntom.SetTwoLiteralDiffBVA(antomSettings->bvaTwoLitDiff);
	myAntom.SetVivification(antomSettings->doVivification);
	myAntom.SetSatConst(antomSettings->satconst);

	myAntom.SetIncompleteMode(antomSettings->incompleteMode);
	myAntom.SetSearchMode(antomSettings->searchMode);
	myAntom.SetNetworktype(antomSettings->networkType);
	myAntom.SetDecStratMode(antomSettings->decStrat);
    myAntom.SetEncode01Mode(antomSettings->encode01, antomSettings->lastPos1);
	myAntom.SetSortSoftClauses(antomSettings->sortSoftClauses);
	myAntom.SetGridMode(antomSettings->bypassGrid);
	myAntom.SetHorizontalWidth(antomSettings->horizontalWidth);
	myAntom.SetCSC(antomSettings->setCSC);
	myAntom.SetSkip(antomSettings->doSkipping);
	myAntom.SetOptTarget(antomSettings->targetOpt);
	myAntom.SetPreciseTarget(antomSettings->preciseTarget);

	myAntom.SetPartialMode(antomSettings->partialMode);
	myAntom.SetSplittedWidth(antomSettings->splittedWidth);
	myAntom.SetMaxWidth(antomSettings->maxWidth);
	myAntom.SetMaxParts(antomSettings->maxParts);
	myAntom.SetRelaxationLits(antomSettings->setRelaxationLits);

	// Weighted MaxSAT related stuff
    if (mode == 30)
    {
        myAntom.SetBaseMode(antomSettings->base);
        myAntom.SetPartitionStrategy(antomSettings->partitionStrategy);
        myAntom.SetHeuristic(antomSettings->groupHeuristic);
        myAntom.SetPercentOff(antomSettings->percentOff);
        myAntom.SetPercentOffreinsert(antomSettings->percentOffReinsert);
        myAntom.SetEqualWeight(antomSettings->equalWeight);
        myAntom.SetSolveAtFirst(antomSettings->solveAtFirst);
        myAntom.SetEncodeStrategy(antomSettings->encodeStrategy);
        myAntom.SetOnlyByTares(antomSettings->onlyByTares);
        myAntom.SetCreateGraphFile(antomSettings->createGraphFile);
        myAntom.SetFeatureTest(antomSettings->featureTest);
        myAntom.SetCascadeDivider(antomSettings->cascadeDivider);
        myAntom.SetSolveTareCascadeOnlyByTares(antomSettings->tareCascadeOnlyByTares);
        myAntom.SetMultipleCascadeDivideStrategy(antomSettings->mcDivideStrategy);
        myAntom.SetSepHiWeight(antomSettings->sepHiWeight);
        myAntom.SetWeightPlusOne(antomSettings->weightPlusOne);
        myAntom.SetMaxBucketSize(antomSettings->maxBucketSize);
        myAntom.SetNOfCascades(antomSettings->nOfCasc);
        myAntom.SetInterimResult(antomSettings->interimResult);
    }
  }
#endif

    // Loads a formula from file, returns FALSE if the formula is unsatisfiable.
  // "type" has to be set according to the type of benchmark:
  // type = 0 --> SAT 
  // type = 1 --> partial MaxSAT
  // type = 2 --> partial weighted MaxSAT
  bool Parser::LoadCNF(const std::string& file, uint32_t& maxIndexOrigVars, Antom& antom)
  {
	uint32_t type(0);
	if (mode >= 10 && mode < 20) { type = 1; }
	if (mode >= 20 && mode < 30 ) { type = 2; }
	if (mode >= 30 && mode < 40 ) { type = 3; }

	// Open the file.
	std::ifstream source;
	source.open(file.c_str());

	// Any problems while opening the file?
	if (!source) 
	  {
		// Output. 
		std::cout << "c Unable to open file" << std::endl
				  << "s UNKNOWN" << std::endl; 
      
		// Return UNKNOWN.
		exit(ANTOM_UNKNOWN); 
	  }

	// Variables.
	std::vector<uint32_t> clause;
	uint32_t maxVars(0);
    uint64_t literal(0);
    uint32_t sign(0);
    uint64_t weight(0);
    uint64_t firstWeight(0);
    uint64_t secondWeight(0);
    uint32_t threshold(0);
	bool moreThanTwoWeights(false);
   	bool firstLit(true); 
	char c('0');

	uint64_t topWeight(0);
	uint64_t minWeight((uint64_t) - 1);
	uint64_t maxWeight(0);
	uint64_t sumOfSoftWeights(0);

	if ( type == 2 )
	  {
		topWeight = 2;
	  }

	// Process the file.
	while (source.good())
	  {
		// Get the next character.
		c = (char)source.get();
   
		// No more clauses?
		if (!source.good())
		  { break; }

		// Statistics?
		if (c == 'p')
		  {
			// Get the next char. 
			c = (char)source.get(); 

			// Remove whitespaces.
			while (c == ' ') 
			  { c = (char)source.get(); }

			// In case of a partial MaxSAT benchmark file, the next character has to be "w".
			if (type == 2 || type == 3)
			  { assert(c = 'w'); c = (char)source.get(); }

			// The next three characters have to be "c", "n", and "f".
			assert(c = 'c');
			c = (char)source.get();
			assert(c = 'n');
			c = (char)source.get(); 
			assert(c = 'f');
			c = (char)source.get(); 

			// Remove whitespaces.
			while (c == ' ') 
			  { c = (char)source.get(); }

			// Let's get the number of variables within the current CNF.
            while (c != ' ' && c != '\n' && c != '\r')
			  { maxVars = (maxVars * 10) + (unsigned int) c - '0'; c = source.get(); }

			// Update "maxIndexOrigVariables".
			maxIndexOrigVars = maxVars; 

			// Remove whitespaces.
            while (c == ' ' || c == '\r')
			  { c = source.get(); }

			// Get top weight in case of weighted partial MaxSAT benchmark file.
			if ( type == 3 )
			  {
				
				// Remove the number of clauses within the current CNF.
                while (c != ' ' && c != '\n' && c != '\r')
				  { c = source.get(); }
				
				// Remove whitespaces.
                while (c == ' ' || c == '\r')
				  { c = source.get(); }

				// Let's get the top weight - to differ later on hard & soft clauses.
                //TOBI: Why to set topweight here to 0. why not set it directly?
                //topWeight = 0;
                while (c != ' ' && c != '\n' && c != '\r')
                  {
					topWeight = (topWeight * 10) + (unsigned int) c - '0';
                    c = source.get();
                  }

				if ( topWeight == 0 )
                  {
                    topWeight = (uint64_t)-1;
                  }
			  }

			// Remove whitespaces.
			while (c == ' ') 
			  { c = (char)source.get(); }

			
			// Set the maximum number of variables the solver has to deal
			antom.SetMaxIndex(maxVars);
		  }
		else
		  {
			// Clause? 
			if (c != 'c' && c != '%')
			  {
				// Reset "clause".
				clause.clear();

				// Initialize "firstLit".
				firstLit = true; 
	      
				// Initialize "threshold".
				threshold = 0; 

				// Do we have to solve a MaxSAT benchmark?
				if (type == 1)
				  {
					// Update "threshold".
					threshold = 1; 
				  }
				// Get the next clause.
				while (true)
				  {
					// Initialization.
					literal = 0;
					sign    = 0; 
					// Remove whitespaces.
                    while (c == ' ' || c == '\r')
					  { c = (char)source.get(); }

					// Have we reached the clause stopper?
					if (c == '0')
					  {
                        if (type == 3 && threshold == 1)
                        {
                            if (clause.size() >= threshold && !antom.AddWeightedSoftClause(clause,weight))
                            {
                                source.close();
                                return false;
                            }
                        }
                        else if( threshold == 1 )
                        {
							// Add "clause" to the clause database of "solver".	  
							if (clause.size() >= threshold && !antom.AddSoftClause(clause,weight))
                            {
                                source.close();
                                return false;
                            }
                        }
						// Add "clause" to the clause database of "solver".	  
						else if (clause.size() > threshold && !antom.AddClause(clause))
                        {
                            source.close();
                            return false;
                        }
						break; 
					  }
					// Let's get the next literal.
					while (c != ' ') 
					  {
						if (c == '-')
						  { sign = 1; }
						else
						  { literal = (literal * 10) + (uint32_t) c - '0'; }
						c = (char)source.get();
					  }

					// Consistency check.
					assert(literal != 0); 

					// In case we are loading a partial MaxSAT benchmark, the first literal
					// specifies whether the current clause is a "soft" or "hard" one.
					if ( (type == 2 || type == 3 ) && firstLit)
					  {
						// Reset "firstLit".
						firstLit = false;
						// Consistency check.
						assert(sign == 0);

						// Soft clause?
						if (literal < topWeight )
						  {
                            // Update "threshold".
                            threshold = 1;

                            weight = literal;

                            // Some information to analyze the formula
                            if (type == 3)
                            {
                                // maxWeight
                                if (literal > maxWeight)
                                {
                                    maxWeight = literal;
                                }

                                if (literal < minWeight)
                                {
                                    minWeight = literal;
                                }

                                if (!moreThanTwoWeights)
                                {
                                    if (firstWeight == 0)
                                    {
                                        firstWeight = literal;
                                    }

                                    if (secondWeight == 0 && firstWeight != literal)
                                    {
                                        secondWeight = literal;
                                    }

                                    // check if it is possibly only MaxSAT = type==2
                                    if (secondWeight != firstWeight && literal != firstWeight && literal != secondWeight)
                                    {
                                        moreThanTwoWeights = true;
                                    }
                                }

								sumOfSoftWeights += weight;
                                //std::cout << "   weight: " << weight << std::endl;
                            }
						  }
					  }
					else
					  {
						// Another consistency check, in this case due to (partial) MaxSAT solving. 
						assert((type != 1 && type != 2 && type != 3) || literal <= maxVars); 
						// Add "literal" to "clause".		      
						clause.push_back((literal << 1) + sign); 
						// Update "maxVars" if necessary.
						if (literal > maxVars)
						  {
							maxVars = literal;
							antom.SetMaxIndex(maxVars);
						  }
					  }
				  }
			  }
		  }

		// Go to the next line of the file. 
        while (c != '\n')
		  { c = (char)source.get(); } 	  
	  }

    // Close the file.
	source.close();

    //TOBI - is maxSatFormula...
    //if (!moreThanTwoWeights && (firstWeight == 1 || secondWeight == 1))
	if ( type == 3 )
	  {
        // TOBI: think it over - what is really needed - delete one part after thinking it through.
//        antom.GetMainCascade()->SetParserValues(moreThanTwoWeights, topWeight, minWeight, maxWeight, sumOfSoftWeights);
        // TOBI: if cascade class is finished then delete
		antom.SetMoreThanTwoWeights(moreThanTwoWeights);
		antom.SetTopWeight(topWeight);
		antom.SetMinWeight(minWeight);
		antom.SetMaxWeight(maxWeight);
		antom.SetSumOfSoftWeights(sumOfSoftWeights);
	  }

	// Everything went fine.
	return true;
  }

  void Parser::PrintSettings(void) const
  {
	std::cout << "c benchmark file.........: " << filename << std::endl;
	if (storeResult)
	  { std::cout << "c result file............: " << resultFile << std::endl; }
	// Output.
	switch (mode)
	  {
	  case  0: std::cout << "c operating mode.........: SAT (multi-threaded: portfolio)"                                                                    << std::endl; break; 
	  case  1: std::cout << "c operating mode.........: single-threaded, SATzilla-like"                                                                     << std::endl; break; 
	  case 10: std::cout << "c operating mode.........: MaxSAT, (multi-threaded: internal portfolio)"                                << std::endl; break; 
	  case 11: std::cout << "c operating mode.........: partialMode MaxSAT, (multi-threaded: internal portfolio"                      << std::endl; break; 
	  case 20: std::cout << "c operating mode.........: partial MaxSAT, (multi-threaded: internal portfolio)"                        << std::endl; break; 
	  case 21: std::cout << "c operating mode.........: partialMode partial MaxSAT, (multi-threaded: internal portfolio)"          << std::endl; break;
          case 30: std::cout << "c operating mode.........: CASCADE partial weighted MaxSAT, satisfiability-based (..:..)"	                        					<< std::endl; break;
	  case 31: std::cout << "c operating mode.........: partial weighted MaxSAT, naive Version, satisfiability-based (..:..)"	                        					<< std::endl; break;
	  case 99: std::cout << "c operating mode.........: collecting data for regression analysis (--> SATzilla-like SAT solving)"                            << std::endl; break; 
	  }
	if (mode >= 30 )
	  {
		antomSettings->application = WEIGHTEDMAXSAT;
	  }
	else if (mode >= 20)
	  {
		antomSettings->application = MAXSAT;
	  }
	antomSettings->Print();
  }
  
  void Parser::PrintUsage(void) const
  {
	// Ouput. 
	std::cout << "c usage: ./antom --mode=<0..99> [--result-file=<file>] [options] <wcnf/cnf>"                                        << std::endl
			  << "c"                                                                                                                  << std::endl
			  << "c mode = 0 --> SAT (multi-threaded: portfolio)"                                                                     << std::endl
			  << "c mode = 1 --> SAT (single-threaded, SATzilla-like)"                                                                << std::endl
			  << "c mode = 10 --> MaxSAT (multi-threaded: internal portfolio)"                                                        << std::endl
			  << "c mode = 20 --> partial MaxSAT (multi-threaded: internal portfolio)"                                                << std::endl
			  << "c mode = 30 --> weighted partial MaxSAT (multi-threaded: internal portfolio)"                                       << std::endl
			  << "c mode = 31 --> naiv weighted partial MaxSAT (multi-threaded: internal portfolio)"                                  << std::endl
			  << "c mode = 99 --> collecting data for regression analysis (--> SATzilla-like SAT solving)"                            << std::endl
			  << "c"                                                                                                                  << std::endl
			  << "c general options:"                                                                                                 << std::endl
			  << "c --threads=<1..8>           --> number of threads for multithreaded mode (default 1)"                              << std::endl 
			  << "c --v | --verbose            --> increase verbosity"                                                                << std::endl
			  << "c --verbose=0                --> reset verbosity"                                                                   << std::endl
			  << "c --cpuLimit=[>=0.0]         --> CPU limit in seconds"                                                              << std::endl
			  << "c                                0.0: no CPU limit at all (default)"                                                << std::endl
			  << "c --memLimit=[>=0]           --> Memory limit in MB"                                                                << std::endl
			  << "c                                0: no memory limit at all (default)"                                               << std::endl
			  << "c --verify=<true/false>      --> verifies model in SAT-case with a second antom-instance"                           << std::endl

			  << "c solver option:"                                                                                                   << std::endl
			  << "c --restart=<luby/glucose>   --> sets the restart strategy to either Luby or Glucose (default: glucose)"            << std::endl
			  << "c --unitFactor=<value>       --> sets the unit factor of both restart strategies to 'value' (default: 6)"           << std::endl
			  << "c --decayFactor=<value>      --> sets the decay factor (variable activities) to 'value' (default: 1.05)"            << std::endl
			  << "c --decision=<0..3>          --> sets the decision strategy to mode 'value' (default: 0)"                           << std::endl
			  << "c                                0: use the cached polarity together with antom's 'polarity toggling scheme'"       << std::endl
			  << "c                                1: use the cached polarity only"                                                   << std::endl
			  << "c                                2: the polarity will be set to FALSE regardless of the cached value"               << std::endl
			  << "c                                3: the polarity will be set to TRUE regardless of the cached value"                << std::endl
			  << "c --ccmin=<none/basic/deep>  --> none, basic or deep minimization of conflict clauses (default: deep)"              << std::endl
			  << "c --ternary=<true/false>     --> enables/disables special treatment for ternary clauses (default: false)"           << std::endl
			  << "c --lhbr=<true/false>        --> enables/disables 'Lazy Hyper Binary Resolution' (default: true)"                   << std::endl

			  << "c preprocessor options:"                                                                                            << std::endl
			  << "c --prepro=<true/false>      --> enables/disables preprocessor (default: false)"                                    << std::endl
			  << "c --inpro=<true/false>       --> enables/disables inprocessor during solve (default: false)"                        << std::endl
			  << "c --maxInpro=<true/false>    --> enables/disables inprocessor during maxsolve (default: false)"                     << std::endl
			  << "c --maxloops=<value>         --> sets the maximum number of preprocessing main loops (default: 5)"                  << std::endl
			  << "c --upla=<true/false>        --> enables/disables 'UPLA' during pre/inprocessing (default: true)"                   << std::endl
			  << "c --subsumption=<true/false> --> enables/disables full subsumption check during pre/inprocessing (default: true)"   << std::endl
			  << "c --varElim=<true/false>     --> enables/disables variable elimination during pre/inprocessing (default: true)"     << std::endl
			  << "c --bce=<true/false>         --> enables/disables 'BCE' during pre/inprocessing (default: true)"                    << std::endl
			  << "c --bva=<true/false>         --> enables/disables 'BVA' during pre/inprocessing (default: true)"                    << std::endl
			  << "c --2litdiff=<true/false>    --> enables/disables two literals difference for 'BVA' (default: false)"               << std::endl
			  << "c --vivify=<true/false>      --> enables/disables 'Vivification' during pre/inprocessing (default: true)"           << std::endl
			  << "c --satconst=<0..2>          --> enables/disables constant checks with SAT (default: 0)"                            << std::endl
			  << "c                                0: false"                                                                          << std::endl
			  << "c                                1: true (only deduction)"                                                          << std::endl
			  << "c                                2: full (with solver calls)"                                                       << std::endl 

			  << "c maxSAT options:"                                                                                                  << std::endl
			  << "c --incomplete=<true/false>  --> enables/disables incomplete mode (default: false)"                                 << std::endl 
			  << "c --search=<0..2>            --> search mode (default: 1)"                                                          << std::endl
			  << "c                                0: Unsat-based"                                                                    << std::endl
			  << "c                                1: Sat-based"                                                                      << std::endl
			  << "c                                2: Binary-based"                                                                   << std::endl 
              << "c --network=<0..4>           --> network type (default: 0)"                                                         << std::endl
			  << "c                                0: Bitonic Sorter"                                                                 << std::endl
			  << "c                                1: Odd-Even-Sorter"                                                                << std::endl
			  << "c                                2: __currently unused__"                                                           << std::endl 
			  << "c                                3: Totalizer"                                                                      << std::endl 
              << "c                                4: Warners"                                                                      << std::endl
			  << "c --decstrat=<0..2>          --> special decision strategy for MaxSAT related variables (default: 2)"               << std::endl
			  << "c                                0: For sorter outputs + relaxation + tare variables"                               << std::endl
			  << "c                                1: For sorter outputs"                                                             << std::endl
			  << "c                                2: none"                                                                           << std::endl
			  << "c --encode01=<true/false/lastPos1>   --> encode complete 01-comparators or just half (defualt:true)"   
			  << "c --sortsoft=<value>         --> sort soft clause parts for partialMode (default: 0)"                               << std::endl
			  << "c                                0: no sorting"                                                                     << std::endl
			  << "c                                1: sort soft clauses due to the larger number of conflicting soft clauses"         << std::endl
			  << "c                                2: sort soft clauses due to the smaller number of conflicting soft clauses"        << std::endl 
			  << "c                                3: random sort"                                                                    << std::endl
			  << "c --gridMode=<value>         --> activates bypass grid (default: 0)"                                                << std::endl 
			  << "c                                0: no grid"                                                                        << std::endl
			  << "c                                1: horizontal bypasses"                                                            << std::endl
			  << "c                                2: vertical bypasses"                                                              << std::endl
			  << "c                                3: bypass grid (horizontal+vertical)"                                              << std::endl
			  << "c --bypassWidth=<value>      --> sets width of horizontal grid (default: 4)"                                        << std::endl 
			  << "c --csc=<true/false/all>     --> enables/disables 'conflicting soft clauses' (default: false)"                      << std::endl 
			  << "c --skip=<true/false>        --> enables/disables skipping of comperators (default: false)"                         << std::endl 

			  << "c partial maxSAT options:"                                                                                          << std::endl
			  << "c --partial=<value>          --> partial mode (default: 0)"                                                         << std::endl
			  << "c                                0: none"                                                                           << std::endl
			  << "c                                1: Depth-first"                                                                    << std::endl
			  << "c                                2: Breadth-first"                                                                  << std::endl 
			  << "c --splitWidth=<value>       --> sets splitting width for partialMode (default: 32)"                                << std::endl 
			  << "c --maxWidth=<value>         --> sets maximum splitting width for partialMode (default: 0 = none)"                  << std::endl 
			  << "c --maxParts=<value>         --> sets maximum number of splitting parts for partialMode (default: 0 = none)"        << std::endl 
			  << "c --relax=<true/false>       --> enables/disables setting of relaxation literals (default: false)"                  << std::endl
			  << "c --targetOpt=<value>        --> solve maxSAT instance until target optimum is reached"                             << std::endl
			  << "c --preciseTarget=<value>    --> try to find target as close as possible. Must be used together wird 'targetOpt'"   << std::endl
			  << "c weighted maxSAT options:"                                                                                         << std::endl
			  << "c --base=<value>             --> base for sorter buckets (default: 2)"                                              << std::endl
			  << "c --onlyByTares=<true/false> --> solve whole cascade only by tares, adding therefore as many buckets as necessary." << std::endl
              << "c                                Works only with encodeStrategy=1 (default)."                                       << std::endl
              << "c --analyze=<true/false>     --> converts formulas to ((Partial) Max) SAT if possible (default: true)"              << std::endl
              << "                                 recognizes if formula has a common divisor - and divides all weights by it"        << std::endl
              << "c --solveAtFirst=<true/false>--> starts with a solver call at the beginning. (default: true)"                       << std::endl
              << "c --partitionStrategy=<value>--> sets the strategy how to fill the buckets (default: 0)"                            << std::endl
              << "c                                0: Standard, fill buckets with ungrouped values"                                   << std::endl
              << "c                                1: Grouping weights - depth problematic"                                           << std::endl
              << "c                                2: Grouping weights - balanced depth"                                              << std::endl
              << "c                                3: Group weights then by biggest repeating entry due to a heuristic"               << std::endl
              << "c     --heuristic=<value>    --> 0: Standard combination of other heuristics"                                       << std::endl
              << "c                                1: number of reduced (ternary clauses * 2 + binary clauses)"                       << std::endl
              << "c                                2: sum of bucket sizes the SC occurs"                                              << std::endl
              << "c                                3: closest size in percentage"                                                     << std::endl
              << "c                                4: the one with least occurences in other buckets"                                 << std::endl
              << "c                                5: how many merges are possible after this merge"                                  << std::endl
              << "c                                6: greatest depth of submerges (till depth 3, then adds possible merges)"          << std::endl
              << "c                                ?: number of reduced ternary clauses calculating all possible subMerges"           << std::endl
              << "c     --percentOff=<0..100>  --> The max percentage two weights differ with partitionStrategy=3 (default: 100)"     << std::endl
              << "c     --percentOffReinsert=<true/false --> reinsert if size is reached again (default: true)."                      << std::endl
              << "c     --equalWeight=<value>  --> merges equal weights every value rounds (default: 0 <- not checking)"              << std::endl
              << "c --encodeStrategy=<value>   --> sets strategy of how to encode the buckets (default: 0)"                           << std::endl
              << "c                                0: Encode all."                                                                    << std::endl
              << "c                                1: Encode only if needed (default)."                                               << std::endl
              << "c --weightPlus1=<true/false> --> Solve always weight+1 instead of minimizing tares."                                << std::endl
              << "c Multiple Cascade (mc) Mode, works only with encodeStrategy=1."                                                    << std::endl
              << "c --mcDivideStrategy=<value> --> sets the strategy for dividing the Softclauses into multiple cascades."            << std::endl
              << "c                                0: solves the cascade in normal mode. (default)"                                   << std::endl
              << "c                                1: sorts the SoftClauses, max SCs in one sub cascade <= cascadeDiv."               << std::endl
              << "c                                2: SCs are randomly distributed, max SCs in one sub cascade <= cascadeDiv."        << std::endl
              << "c                                3: the order of the SCs is not changed!."                                          << std::endl
              << "c                                4: sorts the SoftClauses, max bucket size in one cascade <= cascadeDiv."           << std::endl
              << "c                                5: SCs are randomly distributed, max bucket size in one cascade <= cascadeDiv."    << std::endl
              << "c --sepHiWeight=<true/false> --> if highest weight is bigger than sum of all other SCs, process it seperately."     << std::endl
              << "c                                same if difference between two weights is bigger than 10x."                        << std::endl
              << "c --cascDiv=<value>          --> divider of the SCs or max bucket size"                                             << std::endl
              << "c                                it also sets the limit for the max bucket size, if cascade is reunioned."          << std::endl
              << "c                                if not defined extra."                                                             << std::endl
              << "c --maxBucketSize=<value>    --> Sets the max bucket size, if cascade is reunioned."                                << std::endl
              << "c --nOfCasc=<value>          --> Sets the number of cascades in multiple mode. Stronger than cascDiv."              << std::endl
              << "c --interimResult=<value>    --> Interim results are used while processing Cascade."                                << std::endl
              << "c                                0: No Interim Result is used. (standard)"                                          << std::endl
              << "c                                1: Use results to cut at Top."                                                     << std::endl
              << "c                                2: Use results to cut at Bottom (best with sorted weights)."                       << std::endl
              << "c                                3: Use results to cut at both."                                                    << std::endl
              << "s UNKNOWN"                                                                                                          << std::endl;
  }
}
