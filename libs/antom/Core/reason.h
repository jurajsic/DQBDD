/********************************************************************************************
reason.h -- Copyright (c) 2016-2017, Sven Reimer

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

#ifndef ANTOM_REASON_H_
#define ANTOM_REASON_H_

#include <clause.h>

namespace antom {

  // Reason handler
  // Stores the reason for the implication of a variable
  // Handles binary, ternary and n-nary clauses 
  // Adapted from Cryptominisat
  class Reason 
  {
  public:

	//TODO: introduce Lit class

	// Constructor for binary and n-nary clauses:
	// Stores the second literal or clause position in data1
  Reason(uint32_t lit, bool isBinary) :
	_data1(lit),
	  _data2(0),
	  _clauseType(isBinary?BINARY:NNARY)
		{}

	// Constructor for ternary clause:
	// Stores the second literal in data1
	// Stores the third literal in data2
  Reason(uint32_t lit1, uint32_t lit2) :
	_data1(lit1),
	  _data2(lit2),
	  _clauseType(TERNARY)
	  
	{}


	// No Reason
  Reason(void) :
	_data1(0),
	  _data2(0),
	  _clauseType(NOTYPE)
	{}

	~Reason() = default;

	bool IsBinary(void) const
	{
	  return (_clauseType == BINARY);
	}

	bool IsTernary(void) const
	{
	  return (_clauseType == TERNARY);
	}

	bool IsClause(void) const
	{
	  return (_clauseType == NNARY);
	}

	void SetClause(CRef& ref) 
	{
	  assert(IsClause());
	  _data1 = ref;
	}

	CRef GetClause(void) const
	{
	  assert(IsClause());
	  return (CRef)_data1;
	}

	uint32_t GetSecondLit(void) const
	{
	  assert( IsBinary() || IsTernary() );
	  return _data1;
	}

	uint32_t GetThirdLit(void) const
	{ 
	  assert( IsTernary() );
	  return _data2;
	}

	// Returns true if no reason is set
	bool NoReason(void) const
	{
	  return (_clauseType == NOTYPE);
	}

	// Returns true if reason is set
	bool HasReason(void) const
	{
	  return (_clauseType != NOTYPE);
	}

	// check whether the clause "cref" forces this reason
	bool ForcedBy(const CRef& cref) const 
	{
	  if( IsBinary() || IsTernary() )
		{
		  // To be safe, return "true"
		  return true;
		}
	  if( NoReason() )
		{ 
		  return false; 
		}
	  return (cref == (CRef)_data1);
	}
	
	void ClearReason(void) 
	{
	  _clauseType = NOTYPE;
	  _data1 = 0;
	  _data2 = 0;
	}

	void Print(const ClauseAllocator& ca) const
	{
	  switch(_clauseType)
		{
		case NOTYPE: std::cout << "no reason"; break;
		case BINARY: std::cout << "binary"; break;
		case TERNARY: std::cout << "ternary"; break;
		case NNARY: std::cout << "n-nary"; break;
		default: std::cout << "unvalid reason type"; assert(false); break;
		}
	  std::cout << ": " << std::endl;
	  if( IsBinary() || IsTernary() )
		{
		  std::cout << helper::Lit(GetSecondLit());
		  if( IsTernary() )
			{
			  std::cout << " " << helper::Lit(GetThirdLit());
			}
		  std::cout << std::endl;
		}
	  else if (IsClause())
		{
		  ca[GetClause()].Print(); 
		  assert(!ca[GetClause()].ToDelete() );
		}
	  else
		{
		  std::cout << "no forced Reason" << std::endl;
		}
	}

  private:

	uint32_t _data1:32;
	uint32_t _data2:30;
	ClauseType _clauseType:2;
  };

  // Complete Reason for "analyze
  class ReasonComplete 
  {
  public:

	// Constructor 
  ReasonComplete(const Reason& original, uint32_t conflictLit, const ClauseAllocator& ca) :
	_clauseType(NOTYPE),
	  _clause(NULL),
	  _clauseSize(0)
	{
	  if(original.IsBinary() || original.IsTernary())
		{
		  _lits[0] = conflictLit;
		  _lits[1] = original.GetSecondLit();
		  if( original.IsTernary() )
			{
			  _lits[2] = original.GetThirdLit();
			  _clauseType = TERNARY;
			  _clauseSize = 3;
			}
		  else
			{
			  _lits[2] = 0;
			  _clauseType = BINARY;
			  _clauseSize = 2;
			}
		}
	  else if( original.IsClause() )
		{
		  _lits[0] = _lits[1] = _lits[2] = 0;
		  _clause = ca.GetClause(original.GetClause());
		  _clauseType = NNARY;
		  _clauseSize = _clause->size();
		  assert((*_clause)[0] == conflictLit);
		}
	  else
		{
		  _lits[0] = _lits[1] = _lits[2] = 0;
		}
	}

	~ReasonComplete() = default;

	uint32_t size(void) const
	{
	  return _clauseSize;
	}

	bool IsBinary(void) const
	{
	  return (_clauseType == BINARY);
	}

	bool IsTernary(void) const
	{
	  return (_clauseType == TERNARY);
	}

	bool IsClause(void) const
	{
	  return (_clauseType == NNARY);
	}

	uint32_t operator[](const uint32_t i) const
	{
	  if ( _clauseType == NNARY )
		{
		  assert(_clause != NULL);
		  assert( i < _clauseSize );
		  return (*_clause)[i];
		}
	  else 
		{
		  assert( _clauseType == BINARY || _clauseType == TERNARY );
		  assert( i <= _clauseSize );
		  return _lits[i];
		}
	}

	void Print(void) const
	{
	  switch(_clauseType)
		{
		case NOTYPE: std::cout << "no reason"; break;
		case BINARY: std::cout << "binary"; break;
		case TERNARY: std::cout << "ternary"; break;
		case NNARY: std::cout << "n-nary"; break;
		default: std::cout << "unvalid reason type"; assert(false); break;
		}
	  std::cout << ": " << std::endl;
	  if( IsBinary() || IsTernary() )
		{
		  std::cout << helper::Lit(_lits[0]) << " " << helper::Lit(_lits[1]);
		  if( IsTernary() )
			{
			  std::cout << " " << helper::Lit(_lits[2]);
			}
		  std::cout << std::endl;
		}
	  else if (IsClause() )
		{
		  _clause->PrintLits(); 
		}
	  else
		{
		  std::cout << "no forced Reason" << std::endl;
		}
	}
	

  private:

	ClauseType _clauseType;
	Clause* _clause;
	uint32_t _clauseSize;
	uint32_t _lits[3];

  };
}

#endif
