
/********************************************************************************************
main.cpp -- Copyright (c) 2014, Tobias Schubert, Sven Reimer

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

// Include standard headers.
#include <sys/resource.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <csignal>
#include <string>
#include <vector>
#include <omp.h>

// Include antom related headers.
#include "antom_inc.hpp"

// Terminate antom in case of a segmentation fault.
void SIGSEGV_handler(int signum) 
{
  // Output.
  std::cout << "c segmentation fault (signal " << signum << ")" << std::endl
			<< "s UNKNOWN" << std::endl; 

  // Terminate with UNKNOWN.
  exit(ANTOM_UNKNOWN); 
} 

// An example demonstrating how to use antom.
int main (int argc, char** argv)
{ 
  // Define signal handling functions.
  signal(SIGSEGV,SIGSEGV_handler);

  // Initialization.
  std::cout.precision(2);
  std::cout.setf(std::ios::unitbuf);
  std::cout.setf(std::ios::fixed);
  
  // Output.
  std::cout << "c =================================================================" << std::endl
			<< "c antom, Tobias Schubert, Sven Reimer, University of Freiburg, 2014" << std::endl
			<< "c =================================================================" << std::endl;

  // Initialization.

  // Get a first timestamp.
  double start(omp_get_wtime()); 

  // MaxSAT Mode
  unsigned int mode(1);
  unsigned int result(ANTOM_UNKNOWN); 
  unsigned int optimum(0);
  
  // Initialize an "antom" object.
  antom::Antom myAntom;
 
  // Maximum variable index
  unsigned int maxvar(3);

  // Pre-set maximium variable index
  // This is not needed, but will boost a little bit the "addClause"-Routine
  myAntom.setMaxIndex(maxvar);

  // Use incremental maxSAT mode
  myAntom.setIncrementalMode( true );

  // add some clauses
  std::vector< unsigned int > clause;
  
  // add clause "-1 -2 -3"
  clause.push_back(3);
  clause.push_back(5);
  clause.push_back(7);

  if( !myAntom.addClause(clause) )
	{ return ANTOM_UNSAT; }

  clause.clear();
  // add clause "1 3"
  clause.push_back(2);
  clause.push_back(6);

  if( !myAntom.addClause(clause) )
	{ return ANTOM_UNSAT; }

  clause.clear();
  // add clause "2 3"
  clause.push_back(4);
  clause.push_back(6);

  if( !myAntom.addClause(clause) )
	{ return ANTOM_UNSAT; }

  // Now add the soft clauses
  std::vector< unsigned int > softclause;

  // add softclause "-3"
  softclause.push_back(7);

  if( !myAntom.addSoftClause(softclause) )
	{ return ANTOM_UNSAT; }

  softclause.clear();  

  // add softclause "-1"
  softclause.push_back(3);

  if( !myAntom.addSoftClause(softclause) )
	{ return ANTOM_UNSAT; }

  // Now solve the maxSAT problem
  result = myAntom.maxSolve(optimum, mode);

  // Return optimum (should be equal to 1, i.e., one soft clause can not be satisfied)
  std::cout << "c first run: " << std::endl;
  std::cout << "o " << optimum << std::endl 
			<< "s OPTIMUM FOUND" << std::endl; 

  // Get the satisfying variable assignment.
  const std::vector<unsigned int>& model(myAntom.model());
  
  std::cout << "v ";
  for (unsigned int m = 1; m <= maxvar; ++m)
	{
	  if (model[m] != 0)
		{
		  if ((model[m] & 1) == 1)
			{ std::cout << "-"; }
		  std::cout << (model[m] >> 1) << " "; 
		}
	}
  std::cout << "0" << std::endl; 


  // Now add some more clauses and use maxSAT of antom incrementally
  // Add this point the soft clauses and the pseudo Boolean encoding for the maximization constraints are erased from the database

  maxvar = 5;
  clause.clear();

  // add clause "-3 -4 5"
  clause.push_back(7);
  clause.push_back(9);
  clause.push_back(10);

  if( !myAntom.addClause(clause) )
	{ return ANTOM_UNSAT; }

  clause.clear();
  // add clause "3 -5"
  clause.push_back(6);
  clause.push_back(11);

  if( !myAntom.addClause(clause) )
	{ return ANTOM_UNSAT; }

  clause.clear();
  // add clause "4 -5"
  clause.push_back(8);
  clause.push_back(11);

  if( !myAntom.addClause(clause) )
	{ return ANTOM_UNSAT; }

  // Add new soft clauses

  softclause.clear();  
  // add softclause "-1"
  softclause.push_back(3);

  if( !myAntom.addSoftClause(softclause) )
	{ return ANTOM_UNSAT; }

  softclause.clear();  
  // add softclause "-4"
  softclause.push_back(9);

  if( !myAntom.addSoftClause(softclause) )
	{ return ANTOM_UNSAT; }

  softclause.clear();  
  // add softclause "-5"
  softclause.push_back(11);

  if( !myAntom.addSoftClause(softclause) )
	{ return ANTOM_UNSAT; }

  // Now solve the maxSAT problem again
  result = myAntom.maxSolve(optimum, mode);

  // Return optimum (should be equal to 0, i.e. all soft clauses can be satisfied)
  std::cout << "c second run: " << std::endl;
  std::cout << "o " << optimum << std::endl 
			<< "s OPTIMUM FOUND" << std::endl; 

  std::cout << "v ";
  for (unsigned int m = 1; m <= maxvar; ++m)
	{
	  if (model[m] != 0)
		{
		  if ((model[m] & 1) == 1)
			{ std::cout << "-"; }
		  std::cout << (model[m] >> 1) << " "; 
		}
	}
  std::cout << "0" << std::endl; 


  // Get the CPU time. 
  struct rusage resourcesS;
  getrusage(RUSAGE_SELF, &resourcesS); 
  double timeS((double) resourcesS.ru_utime.tv_sec + 1.e-6 * (double) resourcesS.ru_utime.tv_usec);
  timeS += (double) resourcesS.ru_stime.tv_sec + 1.e-6 * (double) resourcesS.ru_stime.tv_usec;

  // Get the wall clock time.
  double timeW(omp_get_wtime() - start); 
 
  // Output.
  std::cout << "c #ID fastest thread.....: " << myAntom.solvingThread()          << std::endl
			<< "c #variables.............: " << myAntom.variables()              << std::endl
			<< "c #clauses...............: " << myAntom.clauses()                << std::endl
			<< "c #literals..............: " << myAntom.literals()               << std::endl
			<< "c #decisions.............: " << myAntom.decisions()              << std::endl
			<< "c #bcp operations........: " << myAntom.bcps()                   << std::endl
			<< "c #implications..........: " << myAntom.implications()           << std::endl
			<< "c #conflicts.............: " << myAntom.conflicts()              << std::endl
			<< "c #restarts..............: " << myAntom.restarts()               << std::endl
			<< "c #simplifications.......: " << myAntom.simplifications()        << std::endl
			<< "c #synchronizations......: " << myAntom.synchronizations()       << std::endl
			<< "c #lhbr clauses..........: " << myAntom.lhbr()                   << std::endl
			<< "c #learnt unit clauses...: " << myAntom.learntUnitClauses()      << std::endl
			<< "c #learnt binary clauses.: " << myAntom.learntBinaryClauses()    << std::endl
			<< "c #learnt ternary clauses: " << myAntom.learntTernaryClauses()   << std::endl
			<< "c average lbd............: " << myAntom.avgLBD()                 << std::endl
			<< "c average cc length......: " << myAntom.avgCCLength()            << std::endl
			<< "c average dec. level.....: " << myAntom.avgDL()                  << std::endl
			<< "c average lev. cleared...: " << myAntom.avgDLclearedCA()         << std::endl
			<< "c average vars unassigned: " << myAntom.avgVarsUnassignedCA()    << std::endl
			<< "c #inprocessings.........: " << myAntom.inprocessings()          << std::endl
			<< "c cpu time...............: " << timeS << "s"                     << std::endl
			<< "c wall clock time........: " << timeW << "s"                     << std::endl
			<< "c cpu utilization........: " << ((timeS / timeW) * 100.0) << "%" << std::endl
			<< "c ===================================================="          << std::endl;

  if( result == ANTOM_UNKNOWN )
	{
	  std::cout << "s TIMEOUT" << std::endl; 
	  return ANTOM_UNKNOWN; 
	}
  else if( result == ANTOM_UNSAT )
	{
	  std::cout << "s UNSATISFIABLE" << std::endl;
	  return ANTOM_UNSAT; 
	}

  std::cout << "s SATISFIABLE" << std::endl;
  return ANTOM_SAT;
}
