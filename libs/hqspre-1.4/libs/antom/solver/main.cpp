
/********************************************************************************************
main.cpp -- Copyright (c) 2013-2016, Tobias Schubert, Sven Reimer

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
#include <cassert>
#include <csignal>
#include <vector>
#include <fstream>
//#include <omp.h>
#include <unistd.h>

// Include antom related headers.
#include "antom.h"
#include "parser.h"

//#define COMPETITION

// Some more headers. 
void SIGSEGV_handler(int signum);
void printUsage( void );

// An example demonstrating how to use antom.
int main (int argc, char** argv)
{ 
  // Define signal handling functions.
  signal(SIGSEGV,SIGSEGV_handler);
  antom::Parser parser;

  // Initialization.
  std::cout.precision(2);
  std::cout.setf(std::ios::unitbuf);
  std::cout.setf(std::ios::fixed);
  
  // Output.
  std::cout << "c ================================================================================" << std::endl
			<< "c antom; Tobias Schubert, Sven Reimer, Tobias Paxian; University of Freiburg, 2017" << std::endl
			<< "c ================================================================================" << std::endl;

  // Wrong number of command line parameters?
  if (argc < 3)
    {

	  parser.PrintUsage();
      // Return UNKNOWN.
      return ANTOM_UNKNOWN;
    }


  // Initialize an "antom" object.
  antom::Antom myAntom;

  // Parse command lines
  if( parser.ParseCommandline(myAntom, argc,argv) == ANTOM_UNKNOWN )
	{ return ANTOM_UNKNOWN; }

  parser.PrintSettings();

#ifndef NDEBUG
  std::cout << "c running in DEBUG mode" << std::endl;
#endif

  // Get a first timestamp.
#ifdef PARALLEL
  double start(omp_get_wtime());
#endif

  // Set all collected options
  //parser.SetSettings(myAntom);
 
  // Initialization.
  uint32_t maxIndexOrigVars(0); 

  // Load the benchmark file specified by the user.
  // Open the file.

  if (!parser.LoadCNF(parser.filename, maxIndexOrigVars, myAntom))
	{
      // Output.
      std::cout << "s UNSATISFIABLE" << std::endl; 
      
      // Return UNSAT.
      return ANTOM_UNSAT;
    }
  
  // Solve the benchmark specified by the user.
  uint32_t result(ANTOM_UNKNOWN); 

  int64_t optimum(-1);

  // analyze structure if it is really a weighted MaxSAT Formula.
  if (parser.mode == 30 && parser.antomSettings->analyze == true)
	{
      switch ( myAntom.AnalyzeandConvertStructure() )
		{
		case antom::ISSAT:
          std::cout << "Is SAT Formula." << std::endl;
          parser.mode = 0;
          break;
		case antom::ISMAXSAT:
          std::cout << "Is MaxSAT Formula." << std::endl;
          parser.mode = 20;
          break;
		case antom::ISWEIGHTEDMAXSAT: break;
		case antom::CONVERTTOMAXSAT: assert(false);
		case antom::DIVIDEWEIGHTSBYDIVISOR: assert(false);
		}

	}
  else if (parser.mode == 30)
	{
      myAntom.AddSoftClauseVector();
	}

  // First preprocess the formula
  // Do not apply Preprocessing with max-antom, since it's performed within "maxSolve()" incrementally
  if( parser.mode < 10 && parser.antomSettings->doPreprocessing != antom::NOPREPRO )
	{ result = myAntom.Preprocess(); }

  if( result == ANTOM_UNKNOWN )
	{
	  switch(parser.mode)
		{ 
		case  0: result = myAntom.Solve(); break; 
		case  1: myAntom.SetVerbosity(1); result = myAntom.SolveSATzilla(); break; 
		case 10: result = myAntom.MaxSolve(optimum); break;
		case 20: result = myAntom.MaxSolve(optimum); break; 
		case 30: result = myAntom.MaxSolveWeightedPartial(optimum); break;
		case 31: result = myAntom.MaxSolve(optimum); break;
		case 99: myAntom.GetDataRegressionAnalysis(); break; 
		default : std::cout << "unknown mode" << std::endl; return 0; break;
		}
	}

  // Get the CPU time. 
  struct rusage resourcesS;
  getrusage(RUSAGE_SELF, &resourcesS); 
  double timeS((double) resourcesS.ru_utime.tv_sec + 1.e-6 * (double) resourcesS.ru_utime.tv_usec);
  timeS += (double) resourcesS.ru_stime.tv_sec + 1.e-6 * (double) resourcesS.ru_stime.tv_usec;

  // Get the wall clock time.
  int64_t residentset = resourcesS.ru_maxrss*1024;
#ifdef PARALLEL
  double timeW(omp_get_wtime() - start);
#endif

  myAntom.PrintStatistics();
  std::cout << "c cpu time...............: " << timeS << "s"                       << std::endl
#ifdef PARALLEL
			<< "c wall clock time........: " << timeW << "s"                       << std::endl
			<< "c cpu utilization........: " << ((timeS / timeW) * 100.0) << "%"   << std::endl
#endif
			<< "c resident set memory....: " << (residentset/(1024*1024)) << " MB" << std::endl
			<< "c ===================================================="            << std::endl;

  // Satisfiable formula?
  if (result == ANTOM_SAT) 
    {
      // Output.

	  if (!parser.storeResult)
		{ std::cout << "s SATISFIABLE" << std::endl; }
	  
	  if (parser.mode > 1)
		{ 
		  assert(parser.mode >= 10 && parser.mode < 40); 
		  std::cout << "o " << optimum << std::endl 
					<< "s OPTIMUM FOUND" << std::endl; 
		}

      // Get the satisfying variable assignment.
      const std::vector<uint32_t>& model(myAntom.Model());

      // Consistency check.
      assert((myAntom.Variables() == 0) || (maxIndexOrigVars > 0)); 
      assert((myAntom.Variables() == 0) || (maxIndexOrigVars < model.size())); 

      // Should we print or save the model?
      if (!parser.storeResult || parser.mode > 1)
		{
		  /*
		  // Output.
		  std::cout << "v ";
		  for (uint32_t m = 1; m <= maxIndexOrigVars; ++m)
			{
			  if (model[m] != 0)
				{
				  if ((model[m] & 1) == 1)
					{ std::cout << "-"; }
				  std::cout << (model[m] >> 1) << " "; 
				}
			}
		  if (parser.mode <= 1)
			{ std::cout << "0"; } 
		  std::cout << std::endl; 
		  */
		}
      else
		{
		  // Store our solution so that SatELite is able to extend the partial model. 
		  assert(parser.mode <= 1 && parser.storeResult); 
		  std::ofstream satELite(parser.resultFile);
		  satELite << "SAT\n";
		  for (uint32_t m = 1; m <= maxIndexOrigVars; ++m)
			{ 
			  assert(model[m] != 0);
			  if ((model[m] & 1) == 1)
				{ satELite << "-"; }
			  satELite << (model[m] >> 1) << " "; 
			}
		  satELite << "0\n";
		  satELite.close();
		}


	  // Check model with a new and plain antom core
	  if( parser.verify )
		{
		  antom::Antom verify; 
		  verify.SetPreprocessing(antom::NOPREPRO);
		  verify.SetInprocessing(false);

		  //std::vector<uint32_t> verifytriggerVars;
		  uint32_t maxverifyindex(0);
		  if (!parser.LoadCNF(argv[(uint32_t) argc - 1], maxverifyindex, verify))
			{
			  std::cout << "Bad input" << std::endl;
			  assert(false);
			}

		  std::vector< uint32_t > verifyassumptions;
		  for( uint32_t j = 1; j <= maxverifyindex; ++j )
			{
			  if( model[j] != 0 )
				{ 
				  verifyassumptions.push_back(model[j]); 
				}
			}
		  uint32_t verifyresult = verify.Solve(verifyassumptions);

		  if( verifyresult == ANTOM_SAT )
			{
			  std::cout << "Everything okay with model!" << std::endl;
			}
		  else
			{
			  assert( verifyresult == ANTOM_UNSAT );
			  std::cout << "WRONG MODEL" << std::endl;
			  assert(false);
			}
		}

      // Return SAT.
      return ANTOM_SAT; 
    }

  if( result == ANTOM_UNKNOWN )
	{
	  std::cout << "s TIMEOUT" << std::endl; 
	  return ANTOM_UNKNOWN; 
	}
  // Output.
  if (parser.mode != 0 || !parser.storeResult)
    { std::cout << "s UNSATISFIABLE" << std::endl; }

  // Return UNSAT.
  return ANTOM_UNSAT; 
}


// Terminate antom in case of a segmentation fault.
void SIGSEGV_handler(int signum) 
{
  // Output.
  std::cout << "c segmentation fault (signal " << signum << ")" << std::endl
			<< "s UNKNOWN" << std::endl; 

  // Terminate with UNKNOWN.
  exit(ANTOM_UNKNOWN); 
} 
