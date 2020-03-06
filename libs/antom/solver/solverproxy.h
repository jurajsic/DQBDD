/********************************************************************************************
solverproxy.h -- Copyright (c) 2016, Sven Reimer

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

#ifndef ANTOM_SOLVERPROXY_H_
#define ANTOM_SOLVERPROXY_H_

#include <vector>

#include <antom.h>
//#include <cryptominisat5/cryptominisat.h>

namespace antom
{
  
  // Datastructure for managing different back-end solvers
  class SolverProxy
  {
  public:

	// Constructor
	explicit SolverProxy(AntomBase* antom) :
	_antom(antom),
	  _variables(0)
		{}
	virtual ~SolverProxy(void) {};

	//Interfaces to both solvers
	virtual bool AddClause(std::vector<uint32_t>& clause) = 0;
	virtual bool AddUnit(uint32_t lit) = 0;
	virtual uint32_t Solve(void) = 0;
	virtual uint32_t Solve(const std::vector<uint32_t>& assumptions) = 0;
	virtual const std::vector<uint32_t>& Model(void) = 0;
	virtual uint32_t NewVariable(void) = 0;
	virtual void SetMaxIndex(uint32_t var) = 0;
	virtual void SetCPULimit(double limit) = 0;
	virtual void SetPreprocessing(bool val) = 0;
	virtual void SetVerbosity(uint32_t val) = 0;

	virtual SolverType GetSolverType(void) const = 0;

  protected:
	AntomBase* _antom;
	uint32_t _variables;

  private :
	// Copy constructor.
    SolverProxy (const SolverProxy&) = default;

    // Assignment operator.
    SolverProxy& operator = (const SolverProxy&) = default;
  };

  class AntomSolverProxy : public SolverProxy
  {
  public:
	explicit AntomSolverProxy(AntomBase* antom);
	~AntomSolverProxy(void);

	virtual bool AddClause(std::vector<uint32_t>& clause);
	virtual bool AddUnit(uint32_t lit);
	virtual uint32_t Solve(void);
	virtual uint32_t Solve(const std::vector<uint32_t>& assumptions);
	virtual const std::vector<uint32_t>& Model(void);
	virtual uint32_t NewVariable(void);
	virtual void SetMaxIndex(uint32_t var);
	virtual void SetCPULimit(double limit);
	virtual void SetPreprocessing(bool val);
	virtual void SetVerbosity(uint32_t val);

	virtual SolverType GetSolverType(void) const;
  };

#if 0
  class CryptoMiniSatSolverProxy : public SolverProxy
  {
  public:
	CryptoMiniSatSolverProxy(AntomBase* antom);
	virtual ~CryptoMiniSatSolverProxy(void);

	virtual bool AddClause(std::vector<uint32_t>& clause);
	virtual bool AddUnit(uint32_t lit);
	virtual uint32_t Solve(void);
	virtual uint32_t Solve(const std::vector<uint32_t>& assumptions);
	virtual const std::vector<uint32_t>& Model(void);
	virtual uint32_t NewVariable(void);
	virtual void SetMaxIndex(uint32_t var);
	virtual void SetCPULimit(double limit);
	virtual void SetPreprocessing(bool val);
	virtual void SetVerbosity(uint32_t val);
	
	virtual SolverType GetSolverType(void) const;

  protected:
	CMSat::SATSolver* _cryptominisat;
	std::vector<uint32_t> _model;
  };
#endif
}

#endif
