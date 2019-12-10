
/********************************************************************************************
bva.h -- Copyright (c) 2016, Sven Reimer

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

#ifndef ANTOM_BVA_H_
#define ANTOM_BVA_H_

#include "antombase.h"
#include "watcher.h"
#include "statistics.h"
#include "clause.h"
#include "varheap.h"

namespace antom {

  class BVA 
  {
  public:
	explicit BVA(Preprocessor* prepro);

	bool BoundedVariableAddition(void);
	
  private:

	// Copy constructor.
    BVA (const BVA&) = default;

    // Assignment operator.
    BVA& operator = (const BVA&) = default;

	struct LitPair 
	{
	  explicit LitPair(uint32_t a) :
	  lit1(a),
		lit2(0)
	  {}

	  bool operator==(const LitPair& other) const 
	  {
		return (lit1 == other.lit1 && lit2 == other.lit2);
	  }

	  bool operator!=(const LitPair& other) const 
	  {
		return !(*this == other );
	  }

	  size_t hash(const size_t n) const 
	  {
		uint64_t h;
		h = lit1;
		
		if( lit2 == 0 )
		  {
			return (h % n);
		  }

		h = (h<<31) + lit2;
		return (h % n);
	  }

	  uint32_t lit1;
	  uint32_t lit2;
	};

	struct PotentialClause 
	{
	  explicit PotentialClause(const LitPair& lits, WatcherFull c) :
		litpair(lits),
		clause(c)
	  {}

	  bool operator<(const PotentialClause& other) const
	  {
		if( litpair == other.litpair )
		  { 
			return false;
		  }

		if( litpair.lit1 != other.litpair.lit1)
		  {
			return (litpair.lit1 < other.litpair.lit1 );
		  }
		return litpair.lit2 < other.litpair.lit2;
	  }

	  LitPair litpair;
	  WatcherFull clause;
	};

	class TouchList 
	{
	public:
	  TouchList(void) :
		_touchedVar(),
		_isTouched()
	  {}

	  void Touch(const std::vector< uint32_t >& clause)
	  {
		for ( uint32_t i = 0; i != clause.size(); ++i )
		  {
			Touch(clause[i]>>1);
		  }
	  }

	  void Touch(uint32_t var)
	  {
		if( _isTouched.size() <= var )
		  {
			_isTouched.resize(var+1,false);
		  }
		if( !_isTouched[var] )
		  {
			_touchedVar.push_back(var);
			_isTouched[var] = true;
		  }
	  }

	  void clear(void)
		{
		  for( uint32_t i = 0; i != _touchedVar.size(); ++i )
			{
			  _isTouched[_touchedVar[i]] = false;
			}
		  _touchedVar.clear();
		}

	  std::vector< uint32_t > GetList(void) const
	  { 
		return _touchedVar;
	  }

	private:
	  std::vector< uint32_t > _touchedVar;
	  std::vector< bool > _isTouched;
	};


	uint32_t GetLeastOccuring(const WatcherFull& c);
	LitPair GetMostOccuringPotential( uint32_t& largest);

	void RemoveDuplicates(void);

	LitPair DiffClauses(const WatcherFull& c1, const WatcherFull& c2);
	bool SmallerFormula(uint32_t occurences) const;
	int32_t FormulaSize(int32_t litsize, int32_t clssize) const;

	void CollectPotentialClauses(uint32_t nextlit);
	void UpdateQueue(void);

	void CollectClausesForRemove(void);
	CRef GetBVAClause(const std::vector< uint32_t>& clause);

	uint32_t GetNewVariable(void);
	
	bool PerformBVA(void);

	Preprocessor* _prepro;
	Settings* _setting;

	VarHeap< helper::DescendingOrder, size_t > _queue;
	std::vector< unsigned char > _seen;
	std::vector< unsigned char > _seenmlit;
	std::vector< LitPair > _mlits;
	std::vector< LitPair > _mlitscurrent;
	std::vector< WatcherFull > _mcls;
	std::vector< PotentialClause > _potential;
	std::vector< std::vector < uint32_t > > _clausesToRemove;
	TouchList _touched;
	
	ClauseAllocator& _ca;
	Statistics& _statistics;
	bool& _okay;

  };

}

#endif
