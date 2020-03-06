/********************************************************************************************
clause.h -- Copyright (c) 2013-2017, Tobias Schubert, Sven Reimer

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


#ifndef ANTOM_CLAUSE_H_
#define ANTOM_CLAUSE_H_

// Include standard headers.
#include <cassert>
#include <algorithm>
#include <iostream>
#include "allocator.h"
//#include "literal.hpp"

namespace antom
{

  enum ClauseType {
	NOTYPE,
	BINARY,
	TERNARY,
	NNARY,
  };

  typedef RegionAllocator<uint32_t>::Ref CRef;

  // The "Clause" class.
  class Clause
  {
  public:

	// Constructor.
	// Attentention! The memory for the clause has to be allocated in advance, otherwise we run into errors
  Clause(uint32_t* lits, uint32_t lbd, float act, uint32_t length) :
	//_header(),
	  _sign(0),
	  _activity(act)
		{
		  assert(length > 0);
		  _header.reloced = 0;
		  _header.lbd = lbd;
		  _header.length = length;

		  // Copy clause
		  for (uint32_t i = 0; i < length; ++i)
			{
			  GetLiterals()[i] = lits[i];
			}
		}

	// Constructor.
	// Attentention! The memory for the clause has to be allocated in advance, otherwise we run into errors
  Clause(const std::vector<uint32_t>& lits, uint32_t& lbd, float& act, uint32_t& length) :
	  //_header(),
		_sign(0),
		_activity(act)
		{
		  // Consistency check.
		  //assert(_literals != NULL);
		  assert(lbd > 0);
		  assert(length > 0);
		  _header.reloced = 0;
		  _header.lbd = lbd;
		  _header.length = length;

		  // Copy clause
		  for (uint32_t i = 0; i < length; ++i)
			{
			  GetLiterals()[i] = lits[i];
			}
		}

	~Clause(void) {}
	
	uint32_t* GetLiterals(void)
	{ return (uint32_t*)((char*)this + sizeof(Clause)); }

	const uint32_t* GetLiterals(void) const
	{ return (uint32_t*)((char*)this + sizeof(Clause)); }

	uint32_t operator[](size_t i) const
	{ return *(GetLiterals()+i); }
	uint32_t& operator[](size_t i)
	{ return *(GetLiterals()+i); }

	bool operator==(const Clause& other) const 
	{ return GetLiterals() == other.GetLiterals();}

	// Clause will be reallocated from time to time
	bool Reloced(void) const 
	{ return (_header.reloced == 1); }
	CRef Relocation(void) const 
	{ assert(_header.reloced == 1); return GetLiterals()[0]; }
	void Relocate(CRef c) 
	{ _header.reloced = 1; GetLiterals()[0] = (uint32_t)c; }

	// Returns the "Literals Blocks Distance".
	uint32_t Lbd(void) const 
	{ return _header.lbd; }

	// Returns whether the clause was learned
	bool IsLearned (void) const 
	{ return (_header.lbd>1); }

	void MarkForDeletion(void)
	{ _header.lbd = 0; }

	// Returns whether clause is marked as to delete
	bool ToDelete(void) const 
	{ return (bool)(_header.lbd == 0); }

	// Returns the clause length.
	uint32_t size(void) const 
	{ return _header.length; }

	// Returns signature of clause.
	/* uint64_t Sign(void) const  */
	uint32_t Sign(void) const 
	{ return _sign; } 

	void IncreaseActivity(float act)
	{ 
	  _activity += act; 
	  assert(_activity > 0.0);
	}

	float Activity(void) const 
	{ return _activity; }

	void ScaleActivity(void) 
	  //{ _activity *= static_cast<float>(1e-20); }
	{ _activity *= 1e-20; }

	// Sets the "Literals Block Distance" to "lbd".
	void SetLBD(uint32_t lbd) 
	{ 
	  assert(lbd < (1 << 16));
	  _header.lbd = lbd;
	}

	// Sets the clause length to "l".
	void SetLength(uint32_t l) 
	{
	  assert(l > 0); 
	  _header.length = l;
	}

	// Sets the signature of clause to "s".
	//void SetSign(uint64_t s) 
	void SetSign(uint32_t s) 
	{ 
	  assert(s > 0); 
	  _sign = s; 
	}

	int32_t HasLitPolarity(uint32_t lit) const
	{
	  const uint32_t* lits = GetLiterals();
	  for( unsigned int pos = 0; pos != _header.length; ++pos )
		{
		  if( lits[ pos ] == lit )
			{
			  return +1;
			}
		  else if( lits[ pos ] == (lit^1) )
			{
			  return -1;
			}
		}
	  return 0;
	}   

	void CalcSign(void)
	{
	  _sign = 0;

	  uint32_t* lits = GetLiterals();
	  for(uint32_t pos = 0; pos != _header.length; ++pos)
		{
		  //_sign |= 1ULL << (lits[pos] % 64);
		  _sign |= 1 << (lits[pos] % 32);
		}
	} 

	void Sort(void)
	{
	  uint32_t* literals = GetLiterals();
	  std::sort(literals, literals + _header.length );
	}

	// prints literals of clause
	void PrintLits(void) const;
	void Print(void) const;

  private:

	friend class ClauseAllocator;

	struct 
	{
	  unsigned reloced : 1;
	  unsigned length  : 31;
	  unsigned lbd     : 32; 
	} _header;

	// The signature of the clause.
	//uint64_t _sign;
	uint32_t _sign;

	float _activity;
  };

  // A helper function to compare two clauses wrt. their LBD value & clause length.
  bool inline CompareLBD(Clause* c1, Clause* c2) 
  { 
	if (c1->Lbd() == c2->Lbd())
	  { return c1->size() > c2->size(); }
	return c1->Lbd() > c2->Lbd(); 
  }

  // A helper function to compare two clauses wrt. their clause activities & clause length.
  bool inline CompareActivity(Clause* c1, Clause* c2) 
  { 
	if (c1->Activity() == c2->Activity())
	  { return c1->size() > c2->size(); }
	return c1->Activity() < c2->Activity(); 
  }

  // Adopted from minisat

  // Undefined reference for dummy references in case a clause was deleted
  const CRef CRef_Undef = RegionAllocator<uint32_t>::Ref_Undef;	
  // ClauseAllocator -- a simple class for allocating memory for clauses:
  class ClauseAllocator : public RegionAllocator<uint32_t>
  {
  public:

	explicit ClauseAllocator(uint32_t initialCap) : 
	  RegionAllocator<uint32_t>(initialCap)
	{}


	void CopyTo(ClauseAllocator& to)
	{
	  RegionAllocator<uint32_t>::CopyTo(to); 
	}

	void MoveTo(ClauseAllocator& to)
	{
	  RegionAllocator<uint32_t>::MoveTo(to); 
	}

	CRef Alloc(Clause& c)
	{
	  assert(sizeof(uint32_t) == 4);

	  CRef cid = RegionAllocator<uint32_t>::Alloc(ClauseWord32Size(c.size()));
	  new (GetAdress(cid)) Clause(c.GetLiterals(), c.Lbd(), c.Activity(), c.size());

	  return cid;
	}

	CRef Alloc(uint32_t* lits, uint32_t lbd, float act, uint32_t length)
	{
	  assert(sizeof(uint32_t) == 4);

	  CRef cid = RegionAllocator<uint32_t>::Alloc(ClauseWord32Size(length));
	  new (GetAdress(cid)) Clause(lits, lbd, act, length);

	  return cid;
	}

	CRef Alloc(const std::vector<uint32_t>& lits, uint32_t lbd, float act, uint32_t length)
	{
	  assert(sizeof(uint32_t) == 4);

	  CRef cid = RegionAllocator<uint32_t>::Alloc(ClauseWord32Size(length));
	  new (GetAdress(cid)) Clause(lits, lbd, act, length);

	  return cid;
	}

	void AllocRegion(uint32_t size)
	{ RegionAllocator<uint32_t>::AllocRegion(size); }

	Clause& operator[](Ref r) 
	{ return (Clause&)RegionAllocator<uint32_t>::operator[](r); }

	const Clause& operator[](Ref r) const 
	{ return (Clause&)RegionAllocator<uint32_t>::operator[](r); }

	inline Clause* GetClause (Ref r)
	{ return reinterpret_cast<Clause*>(RegionAllocator<uint32_t>::GetAdress(r)); }

	inline Clause* GetClause (Ref r) const
	{ return reinterpret_cast<Clause*>(RegionAllocator<uint32_t>::GetAdress(r)); }

	uint32_t GetReference (const Clause* t) 
	{ return RegionAllocator<uint32_t>::GetReference((uint32_t*)t); }

	void Free(CRef cid)
	{
	  Clause& c = operator[](cid);
	  RegionAllocator<uint32_t>::Free(ClauseWord32Size(c.size()));
	}

	void FreeLiterals(uint32_t numberOfLits) 
	{ RegionAllocator<uint32_t>::Free(numberOfLits); }

	uint32_t GetSize(void) const
	{ return RegionAllocator<uint32_t>::GetSize(); }

	uint32_t GetCapacity(void) const
	{ return RegionAllocator<uint32_t>::GetCapacity(); }
  
	void Reloc(CRef& cr, ClauseAllocator& to)
	{
	  Clause& c = operator[](cr);
	  if (c.Reloced()) 
		{
		  cr = c.Relocation(); 
		  return; 
		}
        
	  cr = to.Alloc(c);
	  c.Relocate(cr);
	}

	void Reset(void)
	{ RegionAllocator<uint32_t>::Reset(); }

  private:
	static uint32_t ClauseWord32Size(uint32_t size)
	{
	  return static_cast<uint32_t>((sizeof(Clause) + (sizeof(uint32_t) * size )) / sizeof(uint32_t)); 
	}
	
  };
}

#endif
