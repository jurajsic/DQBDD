/********************************************************************************************
watcher.h -- Copyright (c) 2016-2017, Sven Reimer

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

#ifndef WATCHER_HPP
#define WATCHER_HPP

#include "helper.h"
#include "clause.h"

namespace antom 
{

  // Watcher class
  // Handles binary and ternary clauses 
  // Adapted from Cryptominisat
  class Watcher
  {
  public:

	// Constructor for binary clause:
	// Stores the second literal in data1
	// Stores the information whether the clause was learnt in data2
  Watcher(uint32_t lit, bool learnt) :
	  _data1(lit),
		_data2(learnt),
		_clauseType(BINARY)
	{}

	// Constructor for ternary and n-nary clause:
	// Stores a second/blocking literal in data1
	// Stores the third literal/clauses position in data2
  Watcher(uint32_t lit1, uint32_t lit2, ClauseType type) :
	  _data1(lit2),
		_data2(lit1), 	  // data2 is used for clause position in case of n-nary clause
		_clauseType(type)
	{
	  assert( type == TERNARY || type == NNARY );
	}

  Watcher(void) :
		_data1(0),
		  _data2(0),
		  _clauseType(NOTYPE)
		  {
		  }

	~Watcher() = default;

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

	ClauseType GetType (void) const
	{ 
	  return _clauseType;
	}

	void SetBlockingLiteral(uint32_t blit)
	{
	  assert(IsClause());
	  _data1 = blit;
	}

	uint32_t GetBlockingLiteral(void) const 
	{
	  assert(IsClause() || IsTernary() );
	  return _data1;
	}

	uint32_t GetSecondLit(void) const
	{
	  return _data1;
	}

	void SetSecondLit(uint32_t lit)
	{
	  assert( IsBinary() || IsTernary() );
	  _data1 = lit;
	}

	uint32_t GetThirdLit(void) const
	{ 
	  assert( IsTernary() );
	  return _data2;
	}

	void SetThirdLit(uint32_t lit)
	{
	  assert( IsTernary() );
	  _data2 = lit;
	}

	void SetClause(CRef cr)
	{
	  assert(IsClause());
	  _data2 = cr;
	}

	CRef GetClause(void) const
	{
	  assert(IsClause());
	  return (CRef)_data2;
	}

	bool IsLearnedBinary(void) const
	{
	  assert(IsBinary());
	  return (bool)(_data2 != 0 );
	}

	void SetLearnedBinary(bool learnt)
	{
	  assert(IsBinary());
	  // Updated learning status should always be "false"
	  assert(learnt == false);
	  _data2 = (uint32_t)learnt;
	}

	bool operator==(const Watcher& other) const 
	{
	  return ( ( _data1 == other._data1 ) && ( _data2 == other._data2 ) && ( _clauseType == other._clauseType ) );
	}

	bool operator!=(const Watcher& other) const 
	{
	  return !(*this==other);
	}

	void Print(const ClauseAllocator& ca) const
	{
	  std::cout << "type: ";
	  switch(_clauseType)
		{
		case BINARY: std::cout << "binary"; break;
		case TERNARY: std::cout << "ternary"; break;
		case NNARY: std::cout << "n-nary"; break;
		default: assert(false); break;
		}
	  std::cout << std::endl;
	  
	  if( IsBinary() || IsTernary() )
		{
		  std::cout << "secondlit: " << helper::Lit(_data1);
		  if( IsBinary() )
			{
			  std::cout << " isLearned: " << _data2 << std::endl;;
			}
		  if( IsTernary() )
			{
			  std::cout << " thirdlit: " << helper::Lit(_data2) << std::endl;
			}
		}
	  else
		{
		  std::cout << "blockinglit: " << helper::Lit(_data1) << " cr: " << _data2 << std::endl;
		  std::cout << "the clause: " << std::endl;
		  ca[_data2].Print();
		}
	  
	}

  private:

	uint32_t _data1:32;
	uint32_t _data2:30;
	ClauseType _clauseType:2;
  };

  // Sorter structure, sorts watchlist in binary, ternary, n-nary
  struct WatchedSorter
  {
    inline bool operator()(const Watcher& w1, const Watcher& w2) const 
	{
  	  assert( w2.IsBinary() || w2.IsTernary() || w2.IsClause() );
	  assert( w1.IsBinary() || w1.IsTernary() || w1.IsClause() );
	  if (w2.IsBinary()) 
		{
		  // Sort binaries by index for duplicate checks
		  if( w1.IsBinary() && (w1.GetSecondLit() < w2.GetSecondLit() ) )
			{
			  return true;
			}
		  else
			{ 
			  return false;
			}
		}
	  //w2 is not binary, but w1 is, so w1 must be first
	  if (w1.IsBinary()) 
		{ return true; }

	  //from now on, none is binary.
	  if (w2.IsTernary()) 
		{ return false; }
	  if (w1.IsTernary()) 
		{ return true; }

	  //from now on, none is binary or ternary
	  //don't bother sorting these
	  return false;
	}

  };

  // Contains the watcher and the watched literal
  struct WatcherFull {

	WatcherFull(uint32_t lit, const Watcher& w, const ClauseAllocator& ca):
	wlit(lit),
	  clauseSize(0),
	  type(w.GetType()),
	  clausePointer(NULL)
	{
	  if ( w.GetType() == NNARY )
		{
		  clausePointer = ca.GetClause(w.GetClause());
		  clauseSize = clausePointer->size();
		  lits[0] = lits[1] = lits[2] = 0;
		}
	  else
		{
		  lits[0] = lit;
		  lits[1] = w.GetSecondLit();
		  if ( w.GetType() == TERNARY )
			{
			  lits[2] = w.GetThirdLit();
			  clauseSize = 3;
			}
		  else 
			{
			  lits[2] = 0;
			  clauseSize = 2;
			}
		}
	}

	WatcherFull(void) :
	wlit(0),
	  clauseSize(0),
	  type(NOTYPE),
	  clausePointer(NULL)
	{
	  lits[0] = lits[1] = lits[2] = 0;
	}

	~WatcherFull() = default;

	uint32_t size(void) const
	{
	  return clauseSize;
	}

	uint32_t operator[](const uint32_t i) const
	{
	  if ( type == NNARY )
		{
		  assert(clausePointer != NULL);
		  assert( i < clauseSize );
		  return (*clausePointer)[i];
		}
	  else 
		{
		  assert( type == BINARY || type == TERNARY );
		  assert( i <= clauseSize );
		  return lits[i];
		}
	}

	bool operator==(const WatcherFull& other) const 
	{
	  return ( (wlit==other.wlit ) && ( ( clausePointer == other.clausePointer ) || (lits[0] == other.lits[0] && lits[1] == other.lits[1] && lits[2] == other.lits[2] ) ) ) ;
	}

	bool operator!=(const WatcherFull& other) const 
	{
	  return !(*this==other);
	}

	void Print(void) const
	{
	  std::cout << "wlit: " << helper::Lit(wlit) << std::endl;
	  if( type == NNARY )
		{
		  clausePointer->Print();
		}
	  else
		{
		  std::cout << helper::Lit(lits[0]) << " " << helper::Lit(lits[1]) << " " << helper::Lit(lits[2]) << std::endl;
		}
	}
	
	friend class WatchedSorterFull;

	uint32_t wlit;
	uint32_t clauseSize;
	uint32_t lits[3];
	ClauseType type;
	Clause* clausePointer;
  };

  struct WatchedSorterFull
  {
	inline bool operator()(const WatcherFull& w1, const WatcherFull& w2) const 
	{
	  ClauseType w1type = w1.type;
	  ClauseType w2type = w2.type;

	  // First solve by type
	  if( ( w1type == BINARY ) && ( w2type != BINARY ) )
		{ 
		  return true;
		}
	  if( ( w2type == BINARY ) && ( w1type != BINARY ) )
		{
		  return false;
		}
	  if( ( w1type == TERNARY ) && ( w2type != TERNARY ) )
		{ 
		  return true;
		}
	  if( ( w2type == TERNARY ) && ( w1type != TERNARY ) )
		{
		  return false;
		}

	  assert( w1type == w2type );
	  // Now we have the same type, now sort by clause content

	  if ( w1.size() != w2.size() )
		{
		  return ( w1.size() < w2.size() );
		}
	  
	  for ( uint32_t i = 0; i < w1.size(); ++i )
		{
		  if ( w1[i] != w2[i] )
			{
			  return ( w1[i] < w2[i] );
			}
		}
	  return false;
	}
  };
  
}

#endif
