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

#include <solverproxy.h>
#include <helper.h>

namespace antom
{
  AntomSolverProxy::AntomSolverProxy(AntomBase* antom) :
	SolverProxy(antom)
  { }
  
  AntomSolverProxy::~AntomSolverProxy(void)
  {}

  bool AntomSolverProxy::AddClause(std::vector<uint32_t>& clause)
  {
	//std::cout << __func__ << "Antom" << std::endl;
	bool result = _antom->AddClause(clause);
	_variables = _antom->Variables();
	return result;
  }

  bool AntomSolverProxy::AddUnit(uint32_t lit)
  {
	//std::cout << __func__ << "Antom" << std::endl;
	std::vector<uint32_t> unitClause {lit};
	bool result = _antom->AddClause(unitClause);
	_variables = _antom->Variables();
	return result;
  }
  
  uint32_t AntomSolverProxy::Solve(void)
  {
	//std::cout << __func__ << "Antom" << std::endl;
	return _antom->Solve();
  }
  uint32_t AntomSolverProxy::Solve(const std::vector<uint32_t>& assumptions)
  {
	//std::cout << __func__ << "Antom" << std::endl;
	return _antom->Solve(assumptions);
  }
  
  const std::vector<unsigned int>& AntomSolverProxy::Model(void)
  {
	//std::cout << __func__ << "Antom" << std::endl;
	return _antom->Model();
  }

  uint32_t AntomSolverProxy::NewVariable(void)
  {
	//std::cout << __func__ << "Antom" << std::endl;
	uint32_t newVar = _antom->NewVariable();
	_variables = newVar;
	return newVar;
  }

  void AntomSolverProxy::SetMaxIndex(uint32_t var)
  {
	//std::cout << __func__ << "Antom" << std::endl;
	_antom->SetMaxIndex(var);
	_variables = _antom->Variables();
  }

  void AntomSolverProxy::SetCPULimit(double limit)
  {
	_antom->SetCPULimit(limit);
  }

  void AntomSolverProxy::SetPreprocessing(bool val)
  {
	_antom->SetPreprocessing(val?PREPROCESS:NOPREPRO);
  }

  void AntomSolverProxy::SetVerbosity(uint32_t val)
  {
	_antom->SetVerbosity(val);
  }

  SolverType AntomSolverProxy::GetSolverType(void) const
  {
	return SolverType::ANTOMSOLVER;
  }

#if 0
  CryptoMiniSatSolverProxy::CryptoMiniSatSolverProxy(AntomBase* antom)
  {
   	_cryptominisat = new CMSat::SATSolver;
	_antom = antom;
  }
  
  CryptoMiniSatSolverProxy::~CryptoMiniSatSolverProxy(void)
  {
	assert(_cryptominisat != NULL);
	delete _cryptominisat;
  }
  
  bool CryptoMiniSatSolverProxy::AddClause(std::vector<uint32_t>& clause)
  {
	//std::cout << __func__ << "Crypto" << std::endl;
	std::vector<CMSat::Lit> cmClause;
	for (uint32_t lit : clause)
	  {
		cmClause.push_back(CMSat::Lit(lit>>1,static_cast<bool>(lit&1)));
		if ((lit>>1) > _variables)
		  {
			_variables = lit>>1;
		  }
	  }
	return _cryptominisat->add_clause(cmClause);
  }

    bool CryptoMiniSatSolverProxy::AddUnit(uint32_t lit)
  {
	//std::cout << __func__ << "Crypto " << helper::Lit(lit) << std::endl;
	std::vector<CMSat::Lit> unitClause {CMSat::Lit(lit>>1,static_cast<bool>(lit&1))};
	if ((lit>>1) > _variables)
	  {
		_variables = lit>>1;
	  }
	bool result = _cryptominisat->add_clause(unitClause);
	return result;
  }
  
  uint32_t CryptoMiniSatSolverProxy::Solve(void)
  {
	std::vector<uint32_t> assumptions;
	return Solve(assumptions);
  }
  
  uint32_t CryptoMiniSatSolverProxy::Solve(const std::vector<uint32_t>& assumptions)
  {
	//std::cout << __func__ << "Crypto" << std::endl;
	std::vector<CMSat::Lit> cmAssumptions;
	for (uint32_t lit : assumptions)
	  {
		cmAssumptions.push_back(CMSat::Lit(lit>>1,static_cast<bool>(lit&1)));
	  }
	CMSat::lbool result = _cryptominisat->solve(&cmAssumptions);
	if (result == CMSat::l_True)
	  {
		return ANTOM_SAT;
	  }
	else if (result == CMSat::l_False)
	  {
		return ANTOM_UNSAT;
	  }
	assert(result == CMSat::l_Undef);
	return ANTOM_UNKNOWN;	
  }

  // TODO: direct model, we may have to change the antom model
  const std::vector<uint32_t>& CryptoMiniSatSolverProxy::Model(void)
  {
	//std::cout << __func__ << "Crypto" << std::endl;
	const std::vector<CMSat::lbool> model = _cryptominisat->get_model();
	_model.resize(model.size(),0);
	for (unsigned int v = 0; v != model.size(); ++v)
	  {
		if(model[v] == CMSat::l_Undef)
		  {
			_model[v] = 0;
		  }
		else
		  {
			_model[v] = (v<<1) + static_cast<uint32_t>(model[v] == CMSat::l_False);
			//std::cout << helper::Lit(_model[v]) << " ";
		  }
	  }
	//std::cout << std::endl;
	return _model;
  }

  // This method does not introduce a physical new variable
  // It just returns the next free index
  uint32_t CryptoMiniSatSolverProxy::NewVariable(void)
  {
	//std::cout << __func__ << "Crypto" << std::endl;
	uint32_t newVar = _antom->NewVariable();
	_cryptominisat->new_var();
	return newVar;
  }

  void CryptoMiniSatSolverProxy::SetMaxIndex(uint32_t var)
  {
	_cryptominisat->new_vars(var-_variables+1);
	_antom->SetMaxIndex(var);
	_variables = _antom->Variables();
  }

  void CryptoMiniSatSolverProxy::SetCPULimit(double limit)
  {
	if (limit == 0.0)
	  {
		_cryptominisat->set_timeout_all_calls(std::numeric_limits<double>::max());
	  }
	else
	  {
		_cryptominisat->set_timeout_all_calls(limit);
	  }
  }

  void CryptoMiniSatSolverProxy::SetPreprocessing(bool val)
  {
	if (!val)
	  {
		_cryptominisat->set_no_simplify();
	  }
  }

  void CryptoMiniSatSolverProxy::SetVerbosity(uint32_t val)
  {
	_cryptominisat->set_verbosity(val);
  }

  SolverType CryptoMiniSatSolverProxy::GetSolverType(void) const
  {
	return SolverType::CRYPTOMINISATSOLVER;
  }
#endif
}
