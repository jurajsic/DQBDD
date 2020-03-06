/********************************************************************************************
varcandidate.h -- Copyright (c) 2016, Sven Reimer

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

#ifndef ANTOM_VARCANDIDATE_H_
#define ANTOM_VARCANDIDATE_H_

#include <limits>

namespace antom 
{
  const static int32_t maxcosts = std::numeric_limits<int>::max();

  // Var candidate class
  // Container for maintaining variable elimination candidates
  class VarCandidate 
  {
  public:
	
	VarCandidate(void) :
	  _recalc(true),
	  _tooLarge(false),
	  _costs(antom::maxcosts),
	  _variable(0)
	{}

	explicit VarCandidate(uint32_t var) :
	  _recalc(true),
	  _tooLarge(false),
	  _costs(antom::maxcosts),
	  _variable(var)
	{}

	~VarCandidate(void)
	{};

	void SetRecalc(bool val) 
	{
	  _recalc = val; 
	}

	bool ToRecalc(void) const
	{ 
	  return _recalc; 
	}

	void SetTooLarge(bool val)
	{ 
	  _tooLarge = val; 
	  _costs = maxcosts;
	}

	bool IsTooLarge(void) const 
	{ return _tooLarge; }

	void SetCosts(int32_t costs) 
	{ _costs = costs; }

	int32_t GetCosts(void) const
	{ return _costs; }

	uint32_t GetVariable(void) const 
	{ return _variable; }

	void Print(void) const 
	{
	  std::cout << "var: " << _variable << " costs: " << _costs << std::endl
				<< "tooLarge: " << _tooLarge << " recalc: " << _recalc << std::endl;
	}

  private:
	bool _recalc;
	bool _tooLarge;
	int32_t _costs;
	uint32_t _variable;
  };

  struct VarCandidateSorter
  {
	inline bool operator()(const VarCandidate* v1, const VarCandidate* v2) const
	{
	  if ( v2->ToRecalc() )
		{ return true; }
	  if ( v1->ToRecalc() )
		{ return false; }
	  if ( v2->IsTooLarge() )
		{ return true; }
	  if ( v1->IsTooLarge() )
		{ return false; }
  
	  if ( v1->GetCosts() < v2->GetCosts() )
		{ 
		  return true; 
		}
	  
	  return false;
	}
  };
}

#endif
