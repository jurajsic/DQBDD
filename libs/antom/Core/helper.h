/********************************************************************************************
helper.h -- Copyright (c) 2013-2017, Tobias Schubert, Sven Reimer

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

#ifndef ANTOM_HELPER_H_
#define ANTOM_HELPER_H_

// Some useful helper functions

#include "clause.h"
#include <string>

namespace antom 
{
  namespace helper 
  {
	template <typename T>
	struct DescendingOrder
	{
	explicit DescendingOrder(const std::vector<T>& a) :
	  act(a)
	  {}
	  
	  bool operator()(const uint32_t& x, const uint32_t& y) const 
	  {
		return (act[x] > act[y]);
	  }
	  const std::vector<T>& act;
	};

	template <typename T>
	struct AscendingOrder
	{
	  explicit AscendingOrder(const std::vector<T>& a) :
	  act(a)
	  {}
	  
	  bool operator()(const uint32_t& x, const uint32_t& y) const 
	  {
		return (act[x] < act[y]);
	  }
	  const std::vector<T>& act;
	};

	/*
	template <typename T>
	struct LBDOrder
	{
	LBDOrder(const std::vector<T>& c) :
	  act(c)
	  {}
	  
	  bool operator()(const uint32_t& x, const uint32_t& y) const
	  {
		return CompareLBD(act[x], act[y]);
	  }
	  
	  const std::vector<T>& act;
	};

	template <typename T>
	  struct ActivityOrder
	  {
		
		
	  ActivityOrder(const std::vector<T>& c) :
		act(c)
	  {}

	  bool operator()(const uint32_t& x, const uint32_t& y) const
	  {
		return CompareActivity(act[x], act[y]);
	  }
	  
	  const std::vector<T>& act;
	};
	*/
	
	
	int inline Lit(uint32_t lit)
	{ return ( (lit>>1) * (((lit&1)!=0)?-1:1)  ); }

	template<typename T1, typename T2>
	  bool inline SortPairBySecond(const std::pair<T1, T2>& i, const std::pair<T1, T2>& j)
	{
	  if( i.second != j.second )
		{ return ( i.second < j.second ); }
	  else
		{ return ( i.first < j.first ); }
	}

	template<typename T1, typename T2>
	  bool inline SortPairBySecondReverse(const std::pair<T1, T2>& i, const std::pair<T1, T2>& j)
	{
	  if( i.second != j.second )
		{ return ( i.second > j.second ); }
	  else
		{ return ( i.first > j.first ); }
	}

	// Helper function for sorter network generation
	uint64_t inline GreatestPowerOfTwoLessThan(uint64_t n)
	{
	  uint32_t k(1);
	  while (k < n)
		{ k <<= 1; }
	  return (k >> 1);
	}

	// Helper function for sorter network generation
	uint64_t inline SmallestPowerOfTwoLargerThan(uint64_t n)
	{
	  uint32_t k(1);
	  while (k < n)
		{ k <<= 1; }
	  return k;
	}

	// Returns the i.th element of the "Luby sequence". See "Optimal Speedup 
	// of Las Vegas Algorithms" (Luby, Sinclair, and Zuckerman) for more details.
	uint32_t inline Luby(uint32_t i)
	{
	  // By definition, "i" has to be greater than 0.
	  assert(i > 0);
      
	  // Variables.
	  uint32_t p(2); 
	  uint32_t iP1(i + 1); 
      
	  // Calculate the i.th element of the Luby sequence.
	  while (p < iP1) 
		{ p = p << 1; }
      
	  // Return the i.th element.
	  if (p == iP1)
		{ return (p >> 1); }
	  return Luby(i - (p >> 1) + 1);
	}

  }
}

#endif
