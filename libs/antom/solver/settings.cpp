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

#include <settings.h>

namespace antom
{
  void Settings::Reset(void)
  {
	// General Settings
	ResetGeneral();
	
	// Core Settings
	ResetCore();

	// Prepro Settings
	ResetPrepro();

	// Antombase Settings
	ResetAntom();
  }

  void Settings::ResetGeneral(void)
  {
	threads = 1;
	verbosity = 0;
	cpuLimit = 0.0;
	memLimit = 0.0;
	application = SAT;	
  }

  void Settings::ResetCore(void)
  {
  	decayFactor = 1.05;
	lubyShift = 6;
	restartStrategy = GLUCOSE;
	restartBlocking = false;
	ccMinimization = DEEP;
	initialPolarity = false;
	simplifyStrategy = ANTOM;
	simplifyActivity = SIMP_LBD;
	decisionStrategy = CACHEANDTOGGLE;
	useTernary = false;
	lhbr = false;
	doPreprocessing = NOPREPRO;
	doInprocessing = false;
  }

  void Settings::ResetPrepro(void)
  {
  	maxLoops = 5;
	doSubsumption = true;
	doUpla = true;
	doVarElimination = true;
	doBce = true;
	doHte = true;
	doHse = true;
	doBva = true;
	bvaTwoLitDiff = false;
	doVivification = true;
	varIncrease = 0;
	satconst = 0;
  }

  void Settings::ResetAntom(void)
  {
  	// MaxSAT Settings
	incrementalMode = 0;
	searchMode = SATBASED;
	partialMode = NONE;
	decStrat = 0;
	splittedWidth = 32;
	maxWidth = 0;
	maxParts = 0;
	horizontalWidth = 4;
	bypassGrid = 0;
	sortSoftClauses = 0;
	setCSC = 0;
	doSkipping = false;
	setRelaxationLits = false;
	incompleteMode = false;
	maxInprocess = false;
	satConst = 0;
	networkType = BITONIC;
	preciseTarget = false;
	targetOpt = -1;

	// Weighted MaxSAT Settings
	analyze = true;
	encode01 = true;
	lastPos1 = false;
	base = 2;
	groupHeuristic = 0;
	percentOff = 100;
	percentOffReinsert = true;
	equalWeight = 0;
	partitionStrategy = NOPARTITION;
	solveAtFirst = true;
	encodeStrategy = ENCODEONLYIFNEEDED;
	createGraphFile = "";
	onlyByTares = false;

	// multiple cascade mode
	mcDivideStrategy = SOLVEINNORMALCASCADEMODE;
	cascadeDivider = 0;
	maxBucketSize = 0;
	nOfCasc = 0;
	tareCascadeOnlyByTares = false;
	sepHiWeight = false;
	weightPlusOne = false;
	  
	// interim results
	interimResult = NOINTERIMRESULT;
	  
	featureTest = false;
  }

  void Settings::Print() const
  {
	std::cout << "c #threads...............: " << threads << std::endl; 
	if( cpuLimit > 0.0 )
	  {
		std::cout << "c CPU limit..............: " << cpuLimit << " s" << std::endl;
	  }
	if( memLimit > 0 )
	  {
		std::cout << "c Memory limit...........: " << memLimit << " MB" << std::endl;
	  }
	std::cout << "c application------------: ";
	switch(application)
	  {
	  case SAT : std::cout << "SAT"; break;
	  case MAXSAT : std::cout << "MaxSAT"; break;
	  case WEIGHTEDMAXSAT : std::cout << "Weighted MaxSAT"; break;	
	  }
	std::cout << std::endl;
	std::cout << "c verbosity..............: " << verbosity << std::endl;
	std::cout << "c ------------------------" << std::endl;
	std::cout << "c core options...........:" << std::endl;

	std::cout << "c restart strategy.......: ";
	if (restartStrategy == LUBY)
	  { std::cout << "luby" << std::endl; }
	else if (restartStrategy == GLUCOSE)
	  { std::cout << "glucose" << std::endl; }
	else
	  { std::cout << "undefined" << std::endl; }

	std::cout << "c restart unit factor....: " << lubyShift << std::endl
			  << "c decay factor...........: " << decayFactor << std::endl;

	std::cout << "c restart blocking.......: " << (restartBlocking?"true":"false") << std::endl;

	std::cout << "c cc minimization........: ";
	if (ccMinimization == NOMIN)
	  { std::cout << "none" << std::endl; }
	else if (ccMinimization == BASIC)
	  { std::cout << "basic" << std::endl; }
	else if (ccMinimization == DEEP)
	  { std::cout << "deep" << std::endl; }
	
	std::cout << "c simplify strategy......: ";
	if (simplifyStrategy == ANTOM)
	  { std::cout << "antom" << std::endl; }
	else if (simplifyStrategy == MINISAT)
	  { std::cout << "minisat" << std::endl; }
	else
	  { std::cout << "undefined" << std::endl; }

	std::cout << "c simplify activity......: ";
	if (simplifyActivity == SIMP_LBD)
	  { std::cout << "lbd" << std::endl; }
	else if (simplifyActivity == SIMP_CONF)
	  { std::cout << "conflicts" << std::endl; }
	else
	  { std::cout << "undefined" << std::endl; }

	std::cout << "c decision strategy......: ";
	switch (decisionStrategy)
	  {
	  case CACHEANDTOGGLE: std::cout << "using cached polarity together with antom's 'polarity toggling scheme'" << std::endl; break;
	  case CACHE: std::cout << "using cached polarity only" << std::endl; break;
	  case ALWAYSFALSE: std::cout << "setting polarity to FALSE regardless of the cached value" << std::endl; break;
	  case ALWAYSTRUE: std::cout << "setting polarity to TRUE regardless of the cached value" << std::endl; break;
	  }

	std::cout << "c initial polarity.......: " << (initialPolarity?"true":"false") << std::endl;

	std::cout << "c use ternary clauses....: " << (useTernary?"true":"false") << std::endl;
	std::cout << "c lhbr...................: " << (lhbr?"true":"false") << std::endl;

	std::cout << "c using preprocessing....: ";
	switch(doPreprocessing)
	  {
	  case NONE: std::cout << "false" << std::endl; break;
	  case PREPROCESS: std::cout << "true" << std::endl; break;
	  case INCREMENTAL: std::cout << "incremental" << std::endl; break;
	  default: assert(false);
	  }
	std::cout << "c using inprocessing.....: " << (doInprocessing?"true":"false") << std::endl;
	std::cout << "c using maxinprocessing..: " << (maxInprocess?"true":"false") << std::endl;
	if( (doPreprocessing != NOPREPRO ) || doInprocessing || maxInprocess )
	  {
		std::cout << "c ------------------------" << std::endl;
		std::cout << "c preprocessing options..:" << std::endl;
		std::cout << "c max prepro loops.......: " << maxLoops << std::endl;
		std::cout << "c UPLA...................: " << (doUpla?"true":"false") << std::endl;
		std::cout << "c subsumption............: " << (doSubsumption?"true":"false") << std::endl;
		std::cout << "c variable elimination...: " << (doVarElimination?"true":"false") << std::endl;
		std::cout << "c bce....................: " << (doBce?"true":"false") << std::endl;
		std::cout << "c hte....................: " << (doHte?"true":"false") << std::endl;
		std::cout << "c hse....................: " << (doHse?"true":"false") << std::endl;
		std::cout << "c bva....................: " << (doBva?"true":"false") << std::endl;
		std::cout << "c two literal diff.......: " << (bvaTwoLitDiff?"true":"false") << std::endl;
		std::cout << "c vivification...........: " << (doVivification?"true":"false") << std::endl;
		std::cout << "c const SAT checks.......: ";
		switch ( satconst )
		  {
		  case 0: std::cout << "false" << std::endl; break;
		  case 1: std::cout << "true" << std::endl; break;
		  case 2: std::cout << "full" << std::endl; break;
		  }
	  }
  
	if ( application == MAXSAT || application == WEIGHTEDMAXSAT )
	  {
		std::cout << "c ------------------------" << std::endl;
		std::cout << "c MaxSAT options.........: " << std::endl;;
		std::cout << "c incomplete mode........: " << (incompleteMode?"true":"false") << std::endl;
		std::cout << "c seach mode.............: ";
		switch ( searchMode )
		  {
		  case UNSATBASED: std::cout << "Unsat-based" << std::endl; break;
		  case SATBASED: std::cout << "Sat-based" << std::endl; break;
		  case BINARYBASED: std::cout << "Binary-based" << std::endl; break;
		  }
		std::cout << "c network type...........: ";
		switch ( networkType )
		  {
		  case BITONIC: std::cout << "Bitonic sorter" << std::endl; break;
		  case ODDEVEN: std::cout << "Odd-Even sorter" << std::endl; break;
		  case TOTALIZER: std::cout << "Totalizer" << std::endl; break;
                  default: break;
		  }
		std::cout << "c decision strategies....: ";
		switch ( decStrat )
		  {
		  case 0: std::cout << "Sorter outputs + relaxation + tare variables" << std::endl; break;
		  case 1: std::cout << "Sorter outputs" << std::endl; break;
		  case 2: std::cout << "None" << std::endl; break;
		  }
        std::cout << "c encode 01..............: " << (encode01?"true":"false") << std::endl;
		std::cout << "c sort soft clauses......: ";
		switch (sortSoftClauses )
		  {
		  case 0: std::cout << "None" << std::endl; break;
		  case 1: std::cout << "Most conflicting first" << std::endl; break;
		  case 2: std::cout << "Least conflicting first" << std::endl; break;
		  case 3: std::cout << "Random" << std::endl; break;
		  }
		std::cout << "c bypasses...............: ";
		switch (bypassGrid)
		  {
		  case 0: std::cout << "None" << std::endl; break;
		  case 1: std::cout << "Horizontal" << std::endl; break;
		  case 2: std::cout << "Vertical" << std::endl; break;
		  case 3: std::cout << "Horizontal and Vertical" << std::endl; break;
		  }
		if ( (bypassGrid == 2 ) || (bypassGrid ==3 ) )
		  {
			std::cout << "c bypass width...........: " << horizontalWidth << std::endl;
		  }
		std::cout << "c find csc...............: ";
		switch ( setCSC )
		  {
		  case 0: std::cout << "false" << std::endl; break;
		  case 1: std::cout << "only for current sorter" << std::endl; break;
		  case 2: std::cout << "true" << std::endl; break;
		  }
		std::cout << "c skip comparators.......: " << (doSkipping?"true":"false") << std::endl;

		std::cout << "c using partial mode.....: ";
		switch ( partialMode )
		  {
		  case NONE: std::cout << "None" << std::endl; break;
		  case DEPTHFIRST: std::cout << "Depth first search" << std::endl; break;
		  case BREADTHFIRST: std::cout << "Breadth first search" << std::endl; break;
		  }

		if ( partialMode != NONE )
		  {
			if( maxWidth != 0 )
			  {
				std::cout << "c max splitting width..: " << maxWidth << std::endl;
			  }
			else
			  {
				std::cout << "c fixed splitting width: " << splittedWidth << std::endl;
			  }

			if( maxParts != 0 )
			  {
				std::cout << "c max splitting parts.: " << maxParts << std::endl;
			  }
			std::cout << "c use fixed relax lits..: " << (setRelaxationLits?"true":"false") << std::endl;
		  }

		if (application == WEIGHTEDMAXSAT)
		  {
			std::cout << "c base...................: " << base << std::endl;
			std::cout << "c analyze................: " << (analyze?"true":"false") << std::endl;
			std::cout << "c partitionStrategy......: " << partitionStrategy << std::endl;
		  }
		if ( partitionStrategy == GROUPBYBIGGESTREPEATINGENTRY )
		  {
			std::cout << "c   heuristic............: " << groupHeuristic << std::endl;
		  }
		if (targetOpt >= 0)
		  {
			std::cout << "c defined target optimum.: " << targetOpt << std::endl;
			std::cout << "c find precise target....: " << preciseTarget << std::endl;
		  }
	  }
  }
}
