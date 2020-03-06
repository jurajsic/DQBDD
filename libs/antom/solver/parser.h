/********************************************************************************************
parser.hpp -- Copyright (c) 2013-2017, Tobias Schubert, Sven Reimer, Tobias Paxian

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

#ifndef ANTOM_PARSER_HPP
#define ANTOM_PARSER_HPP

#include "antom.h"
#include <fstream>

namespace antom 
{
  struct Parser 
  {

  public:

	Parser(void);
	~Parser(void) = default;

	uint32_t ParseCommandline(antom::Antom& myAntom, int argc, char** argv);

	void SetSettings(antom::Antom& myAntom) const;
  
	bool LoadCNF(const std::string& file, uint32_t& maxIndexOrigVars, Antom& antom);

	void PrintSettings(void) const;
	void PrintUsage(void) const;

	Settings* antomSettings;
	
	// some "meta" or binary specific settings
	bool storeResult;
	std::ifstream source;
	std::string filename;
	std::string resultFile;
	uint32_t mode;
	bool verify;
	SolverType solver;

	// storing assumptions
	// assumptions for current (incremental) call
	std::vector< uint32_t > assumptions;
	uint32_t instance;
	bool doIncrementalSolve;
	bool doSolve;
	bool doReset;

  private:
	// Copy constructor.
    Parser (const Parser&) = default;

    // Assignment operator.
    Parser& operator = (const Parser&) = default;
  };
}

#endif
