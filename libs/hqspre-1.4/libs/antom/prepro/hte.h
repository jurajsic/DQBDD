/********************************************************************************************
hte.h -- Copyright (c) 2016, Sven Reimer

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

#ifndef ANTOM_HTE_H_
#define ANTOM_HTE_H_

#include "antombase.h"
#include "watcher.h"
#include "statistics.h"
#include "clause.h"
#include "varheap.h"

namespace antom {

  class HTE 
  {
  public:

	explicit HTE(Preprocessor* prepro);
	void DoHiddenTautologyClauseElimination(void);

  private:

	// Copy constructor.
    HTE (const HTE&) = default;

    // Assignment operator.
    HTE& operator = (const HTE&) = default;

	//bool CheckHiddenTautology(std::vector< uint32_t >& clause, uint64_t& sign);
	bool CheckHiddenTautology(std::vector< uint32_t >& clause, uint32_t& sign);

	Settings* _setting;
	Preprocessor* _prepro;
	Statistics& _statistics;
	bool _okay;

	VarHeap< helper::AscendingOrder, size_t > _candidates;
	std::vector< bool > _seen;
  };

}

#endif
