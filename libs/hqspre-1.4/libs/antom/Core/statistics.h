/********************************************************************************************
statistiscs.hpp -- Copyright (c) 2016-2017, Tobias Schubert, Sven Reimer

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

#ifndef ANTOM_STATISTICS_H_
#define ANTOM_STATISTICS_H_


namespace antom 
{
  // Container for all solver statistics
  struct Statistics 
  {
	uint32_t decisions; // The number of decisions made so far.
	uint32_t bcps; // The number of BCP operations performed so far. 
	uint64_t implications; // The number of implications found so far.
	uint32_t usedVariables; // The number of acutual used Variables (not assigend on level 0)
	uint32_t totalConflicts; // The number of conflicts encountered so far.
	uint32_t totalLearntUnitClauses; // The number of unit clauses deduced due to conflict analysis.
	uint32_t totalLearntBinaryClauses;  // The number of binary clauses deduced due to conflict analysis.
	uint32_t totalLearntTernaryClauses; // The number of ternary clauses deduced due to conflict analysis.

	// current number of binary, ternary and nary clauses
	// Binary clauses are counted twice
	uint32_t currentBinaryClauses; 
	uint32_t currentTernaryClauses;
	uint32_t currentNaryClauses;
	uint32_t staticClauses; // current number of static clauses in the database (including binary)
	uint32_t learnedBinary;
	uint32_t learnedClauses; // current number of learnt clauses in the database (including binary)
	uint32_t minimizedLiterals; // number of minized literals in conflict analysis
	uint32_t restarts; // The number of restarts performed so far.
	uint32_t blockedRestarts;
	uint32_t simplifications; // The number of database simplifications performed so far.
	uint32_t lhbr; // The number of binary clauses deduced due to "Lazy Hyper Binary Resolution".

	uint32_t inprocessings; // The number of inprocessings performed so far.

	uint32_t synchronizations; // The number of synchronizations performed so far.
	double progress; // The progress after the search process has been stopped due to reaching the limit wrt. synchronizations.

	// Divided by "m_conflicts", "m_globalAvgLBD" represents the average 
	// "Literals Block Distance" of all conflict clauses deduced so far. 
	uint64_t globalLBD; 

	uint64_t localLBD;
	uint64_t localConflicts;

	// Divided by "m_conflicts", "m_avgCCLength" represents the average length 
	// of all conflict clauses deduced so far.
	uint64_t staticLength;
	uint64_t learnedLength;
	uint64_t totalLearnedLength; 

	// Divided by "m_conflicts + m_restarts + 1", "m_avgDL" represents 
	// the solver's average decision level before backtracking.
	uint64_t DL; 

	// Divided by "m_conflicts", "m_avgDLclearedCA" represents the 
	// average number of decision levels cleared during conflict analysis. 
	uint64_t DLclearedCA;
	  
	// Divided by "m_conflicts", "m_avgVarsUnassignedCA" represents the 
	// average number of variables getting unassigned during conflict analysis. 
	uint64_t VarsUnassignedCA; 

	// Some preprocessor statistics
	uint32_t constantVariables;
	uint32_t constantVariablesBySAT;
	uint32_t equivalentVariables;
	uint32_t uplaConstantVariables;
	uint32_t uplaEquivalentVariables;	
	uint32_t resolvedVariables;
	uint32_t resolvedLiterals;
	uint32_t blockedClauses;
	uint32_t hiddenTautologies;
	uint32_t hiddenSubsumptions;
	uint32_t monotoneVariables;
	uint32_t dontcareVariables;
	uint32_t subsumptions;
	uint32_t selfSubsumptions;
	uint32_t bvaVariables;
	uint32_t bvaLiterals;
	uint32_t vivifySubsumptions;
	uint32_t vivifyUnits;
	uint32_t vivifyDiff;
	uint32_t unitPropagations;

	double runtime_upla;
	double runtime_subsumption;
	double runtime_varElim;
	double runtime_bce;
	double runtime_bva;
	double runtime_hte;
	double runtime_vivify;
	double runtime_preprocessing;

	Statistics(void) :
	decisions(0),
	  bcps(0),
	  implications(0),
	  usedVariables(0),
	  totalConflicts(0),
	  totalLearntUnitClauses(0),
	  totalLearntBinaryClauses(0),
	  totalLearntTernaryClauses(0),
	  currentBinaryClauses(0),
	  currentTernaryClauses(0),
	  currentNaryClauses(0),
	  staticClauses(0),
	  learnedBinary(0),
	  learnedClauses(0),
	  minimizedLiterals(0),
	  restarts(0),
	  blockedRestarts(0),
	  simplifications(0),
	  lhbr(0),
	  inprocessings(0),
	  synchronizations(0),
	  progress(0.0),
	  globalLBD(0),
	  localLBD(0),
	  localConflicts(0),
	  staticLength(0),
	  learnedLength(0),
	  totalLearnedLength(0),
	  DL(0),
	  DLclearedCA(0),
	  VarsUnassignedCA(0),
	  constantVariables(0),
	  constantVariablesBySAT(0),
	  equivalentVariables(0),
	  uplaConstantVariables(0),
	  uplaEquivalentVariables(0),
	  resolvedVariables(0),
	  resolvedLiterals(0),
	  blockedClauses(0),
	  hiddenTautologies(0),
	  hiddenSubsumptions(0),
	  monotoneVariables(0),
	  dontcareVariables(0),
	  subsumptions(0),
	  selfSubsumptions(0),
	  bvaVariables(0),
	  bvaLiterals(0),
	  vivifySubsumptions(0),
	  vivifyUnits(0),
	  vivifyDiff(0),
	  unitPropagations(0),
	  runtime_upla(0.0),
	  runtime_subsumption(0.0),
	  runtime_varElim(0.0),
	  runtime_bce(0.0),
	  runtime_bva(0.0),
	  runtime_hte(0.0),
	  runtime_vivify(0.0),
	  runtime_preprocessing(0.0)
	{}

	void ClearClauseStatistics(void)
	{
	  usedVariables = 0;
	  currentBinaryClauses = 0;
	  currentTernaryClauses = 0;
	  currentNaryClauses = 0;
	  staticClauses = 0;
	  staticLength = 0;
	  learnedBinary = 0;
	  learnedClauses = 0;
	  learnedLength = 0;
	}
	void ResetCore(void)
	{
	  decisions = 0;
	  implications = 0;
	  bcps = 0;
	  usedVariables = 0;
	  currentBinaryClauses = 0;
	  currentTernaryClauses = 0;
	  currentNaryClauses = 0;
	  totalConflicts = 0;
	  totalLearntUnitClauses = 0;
	  totalLearntBinaryClauses = 0;
	  totalLearntTernaryClauses = 0;
	  staticClauses = 0;
	  learnedBinary = 0;
	  learnedClauses = 0;
	  minimizedLiterals = 0;
	  restarts = 0;
	  blockedRestarts = 0;
	  simplifications = 0;
	  lhbr = 0;
	  inprocessings = 0;
	  synchronizations = 0;
	  progress = 0.0;
	  globalLBD = 0;
	  localLBD = 0;
	  localConflicts = 0;
	  staticLength = 0;
	  learnedLength = 0;
	  totalLearnedLength = 0;
	  DL = 0;
	  DLclearedCA = 0;
	  VarsUnassignedCA = 0;
	}

	void ResetPrepro(void)
	{
	  constantVariables = 0;
	  constantVariablesBySAT = 0;
	  equivalentVariables = 0;
	  uplaConstantVariables = 0;
	  uplaEquivalentVariables = 0;	
	  resolvedVariables = 0;
	  resolvedLiterals = 0;
	  blockedClauses = 0;
	  hiddenTautologies = 0;
	  hiddenSubsumptions = 0;
	  monotoneVariables = 0;
	  dontcareVariables = 0;
	  subsumptions = 0;
	  selfSubsumptions = 0;
	  bvaVariables = 0;
	  bvaLiterals = 0;
	  vivifySubsumptions = 0;
	  vivifyUnits = 0;
	  vivifyDiff = 0;
	  unitPropagations = 0;

	  runtime_upla = 0.0;
	  runtime_subsumption = 0.0;
	  runtime_varElim = 0.0;
	  runtime_bce = 0.0;
	  runtime_bva = 0.0;
	  runtime_hte = 0.0;
	  runtime_vivify = 0.0;
	  runtime_preprocessing = 0.0;
	}

	uint32_t inline GetNaryLearnedClauses(void) const 
	{
	  return totalConflicts-totalLearntUnitClauses-totalLearntBinaryClauses-totalLearntTernaryClauses;
	}

	double inline LitPerStaticClauses(void) const
	{
	  if (staticClauses == 0)
		{ return 0.0; }
	  return ((double)staticLength/staticClauses);
	}

	double inline LitPerConflicts(void) const
	{
	  if (totalConflicts == 0)
		{ return 0.0; }
	  return ((double)learnedLength/totalConflicts);
	}
	  
	double inline CurrentLitPerConflicts(void) const
	{ 
	  if (learnedClauses == 0)
		{ return 0.0; }
	  return ((double)learnedLength/learnedClauses); 
	}

	double inline AvgLBD(void) const 
	{
	  if (totalConflicts == 0)
		{ return 0.0; }
		
	  return ((double) globalLBD / totalConflicts); 
	}

	double inline AvgCCLength(void) const
	{
	  if (totalConflicts == 0)
		{ return 0.0; }
	  return ((double) totalLearnedLength / totalConflicts); 
	}

	double inline CurrentCCLength(void) const
	{
	  if (learnedClauses == 0)
		{ return 0.0; }
	  return ((double) learnedLength / learnedClauses); 
	}

	double inline AvgDL(void) const 
	{ return ((double) DL / (totalConflicts + restarts + 1)); }

	double inline AvgDLclearedCA(void) const
	{
	  if (totalConflicts == 0)
		{ return 0.0; }
	  return ((double) DLclearedCA / totalConflicts); 
	}

	// Returns the average number of variables getting unassigned during conflict analysis. 
	double inline AvgVarsUnassignedCA(void) const
	{
	  if (totalConflicts == 0)
		{ return 0.0; }
	  return ((double) VarsUnassignedCA / totalConflicts); 
	}

	void PrintClauseStats(void) const
	{
	  std::cout << "c #Used variables        : " << usedVariables << std::endl
				<< "c #Bin clauses           : " << (currentBinaryClauses>>1) << std::endl
				<< "c #Ternary clauses       : " << currentTernaryClauses << std::endl
				<< "c #N-nary clauses        : " << currentNaryClauses << std::endl
				<< "c #Total clauses         : " << ((currentBinaryClauses>>1) + currentTernaryClauses + currentNaryClauses ) << std::endl
				<< "c ----------------------------" << std::endl;
	}
	  
  };
}

#endif
