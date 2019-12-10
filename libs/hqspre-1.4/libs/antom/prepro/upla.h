
/********************************************************************************************
upla.h -- Copyright (c) 2016, Sven Reimer

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

#ifndef ANTOM_UPLA_H_
#define ANTOM_UPLA_H_

#include "antombase.h"
#include "watcher.h"
#include "statistics.h"
#include "clause.h"
#include "varheap.h"

namespace antom 
{

  class ModelRebuilder;

  class UPLA
  {
  public:
	explicit UPLA(Preprocessor* prepro);

	bool DoUpla(bool incremental = false);

	friend class Preprocessor;
	
  protected:
	// Copy constructor.
    UPLA (const UPLA&) = default;

    // Assignment operator.
    UPLA& operator = (const UPLA&) = default;
	
	void Backtrack(uint32_t declit);
	void AddUnit(uint32_t lit);

	Core* _core;
	Preprocessor* _prepro;
	Settings* _setting;
	ModelRebuilder* _rebuilder;

	std::vector< bool > _nextcandidates;
	std::vector< std::pair< uint32_t, uint32_t> > _newEquivalences;
	std::vector< uint32_t > _posimplications;
	std::vector< uint32_t > _negimplications;

	std::vector< uint32_t > _activity;
	VarHeap< helper::DescendingOrder, uint32_t > _heap;

	uint32_t& _dsEndIndex;
	uint32_t& _dsImplIndex;
	std::vector< unsigned char >& _assignment;

	ClauseAllocator& _ca;
	Statistics& _statistics;
	bool& _okay;
  };

}

#endif
