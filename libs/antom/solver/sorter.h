/********************************************************************************************
sorter.h -- Copyright (c) 2016, Sven Reimer

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

#ifndef ANTOM_SORTER_H_
#define ANTOM_SORTER_H_

#include <vector>

#include "helper.h"
#include "softclause.h"
#include "antom.h"
//#include "totalizerencodetree.h"

namespace antom
{
  class Antom;
  class Control;
  struct TotalizerEncodeTree;

  // Datastructure for managing sorter structures
  class Sorter 
  {

  public:

	// Constructor
	Sorter(uint32_t size, Antom* antom);

	~Sorter(void);

	uint32_t size (void) const 
	{ return _outputs.size(); }
	uint32_t NumberOfSoftclauses(void) const 
	{ return _softClauses.size(); }
	uint32_t NumberOfTares(void) const 
	{ return _tare.size(); }
    uint64_t Weight(void) const
	{ return _sumOfWeights; }

    void AddWeight(uint64_t weight)
	{
	  _sumOfWeights += weight;
	}

	std::vector< uint32_t > GetTares(void) const 
	  {
		return _tare;
	  }

    void AddOutput(uint32_t var)
      {
        _outputs.push_back(var);
      }

    void AddPreprocessingOutput(uint32_t var, uint64_t howOften)
      {
        _preprocessingOutputs.push_back(var);

        // TODO: for later use with higher bases - to encode them more efficently!
        _numberOfOutput.push_back(howOften);
      }

    void AddSortedVector(std::vector<uint32_t> vec)
      {
        _sortedVectors.push_back(vec);
      }

    void SortSortedVectorBySize()
      {
        std::sort(_sortedVectors.begin(), _sortedVectors.end(),
                  [](const std::vector<uint32_t> & a, const std::vector<uint32_t> & b)
                        { return a.size() > b.size(); });
      }

	std::vector< uint32_t > GetOutputs(void) const
	  {
		return _outputs;
	  }

	void AddSoftClause(SoftClause* sc)
	{
	  _softClauses.push_back(sc);
	}

	std::vector< SoftClause* > GetSoftClauses(void) const 
	  {
		return _softClauses;
	  }

	void AddSoftClauseToSorter(SoftClause* sc)
	{
	  for (unsigned int w = 0; w != sc->weight; ++w)
		{
          AddOutput(sc->relaxationLit >> 1);
		}
	  AddWeight(sc->weight);
	  AddSoftClause(sc);
	}

    void AddTare(uint32_t var, uint32_t partitionStrategy)
	{
	  _tare.push_back(var);
      if (partitionStrategy == 2)
      {
        _preprocessingOutputs.push_back(var);
      }
      else
      {
        _outputs.push_back(var);
      }
	}

	void SetTarePosition(uint32_t pos)
	{ 
	  _tarePosition = pos;
	}

	uint32_t GetTarePosition(void) const 
	{ 
	  return _tarePosition;
	}

	void IncreaseTarePosition(void) 
	{
	  ++_tarePosition;
	}

    void SetMinSatisfied(uint32_t val)
	{
	  _minSatisfied = val;
	}

    uint32_t GetMinSatisfied(void) const
	{
	  return _minSatisfied;
	}

	void SetProceeded(bool val)
	{
	  _proceeded = val;
	}

	bool GetProceeded(void) const
	{
	  return _proceeded;
	}

	void SetProceedNext(bool val)
	{
	  _proceedNext = val;
	}

	bool GetProceedNext(void) const
	{
	  return _proceedNext;
	}

	uint32_t GetCurrentTare(void) const;

	uint32_t CurrentGlobalAssumption(void);

	void AddBypassClause(uint32_t lit1, uint32_t lit2);

	void FindCSC(void);
	void SetSorterCSC(void);

	void MergeWithSorter(Sorter& sorter);
    // needed for the bucket version without triggervar check!
    void SimpleMergeWithSorter(Sorter &sorter);

	void SetVerticalBypasses(void);
	void SetHorizontalBypasses(const std::vector< uint32_t >& outputA, const std::vector< uint32_t >& outputB);

    void EncodeSorter(uint64_t weight);
    void EncodeSorterWithGroups();
    void InsertSortedVector(std::vector< uint32_t >& X);
    void DumpSortedVector();

	void BitonicSort(uint32_t lo, uint32_t n, bool dir);
	void BitonicMerge(uint32_t lo, uint32_t n, bool dir);
	void BitonicCompare(uint32_t i, uint32_t j, bool dir);

	void OddEvenSort(uint32_t lo, uint32_t n);
	void OddEvenMerge(uint32_t lo, uint32_t n, uint32_t r);
	void OddEvenCompare(uint32_t i, uint32_t j);

	void MergeTotalizer(std::vector< uint32_t >& X);
	void TotalizerPhi(uint32_t lo, uint32_t hi);
    void TotalizerPhiChosenOutputs(uint32_t lo, uint32_t hi, std::vector<uint32_t>* outputList = NULL);
    // Encodes in ascending order of input A.
    void EncodingUV(uint32_t beginA, uint32_t endA, uint32_t endB);
    // Encodes in ascending order of the output vector.
    void EncodingUVDiagonal(uint32_t beginA, uint32_t endA, uint32_t endB, std::vector<uint32_t>* outputList = NULL);
    inline void AddClausesToIndex(bool direction, uint32_t outputInd, uint32_t sizeA, uint32_t beginA, uint32_t sizeB, uint32_t endA, uint32_t endB, uint32_t outputVar, uint32_t beginB);
    // Totalizer Structure to Encode only requested outputs.
    void CreateTotalizerEncodeTree();
    // Encodes just one position.
    uint32_t TotalizerEncodeOutput(TotalizerEncodeTree* tree, uint32_t position);

	// to generate an unary vector containing only the literals with
    // weight base^(i+1)
	std::vector<uint32_t> BaseRanks(uint32_t baseMode) const;

    /**
     * @brief GetOrEncodeOutput
     *          Depending on encodeStrategy and if it is yet encoded.
     * @param position
     * @return requested _outputs position or encoded _outputTreePosition
     */
    uint32_t GetOrEncodeOutput(uint32_t position, bool encodeOnlyOnes = false);

	void SetSolverProxy(SolverProxy* solver);

	void Print(void) const;

    friend class MultipleCascade;
    friend class Bucket;
    friend class Cascade;
    friend class MultipleCascade;
    uint32_t TotalizerEncodeOnes(TotalizerEncodeTree *tree, uint32_t outputIndex, uint32_t outputVar);
  private:

	// Copy constructor.
    Sorter (const Sorter&) = default;

    // Assignment operator.
    Sorter& operator = (const Sorter&) = default;
	
	Antom* _antom;
	Settings* _setting;

    // current outputs of sorter
	std::vector< uint32_t > _outputs;

    // totalizer needs to save in between outputs if it encodes only necessary outputs.
    // therefore this Tree is needed
    TotalizerEncodeTree* _outputTree;

    // for using grouping strategies for the sorter
    std::vector< uint32_t > _preprocessingOutputs;
    // TOBI: get a pointer solution for _sortedVectors[i]
    // -> better work with pointer - problems with Sorting...
    // unary sorted groups of soft clauses - used in more than one bucket
    std::vector< std::vector< uint32_t > > _sortedVectors;
	//-> std::vector< Sorter* > _sortedVectors;
    // for base > 2 it occurs that one soft clause occurs more than once
    std::vector< uint32_t > _numberOfOutput;
	// List of soft clauses which are currently connected to this sorter
	std::vector< SoftClause* > _softClauses;

	// List of tare variabels for weighted mode
	std::vector< uint32_t > _tare;

	uint32_t _minContra;

	// The local optimum
    uint64_t _depth;

	// Minimal number of unsatisfied soft clauses in this sorter
	uint32_t _minUnsatisfied;
	// Minimal number of satisfied soft clauses
   	uint32_t _minSatisfied;

	uint32_t _fakeLiterals;
	
    uint64_t _sumOfWeights;

	// Calculated local optimum for this sorter?
	bool _proceeded;
	bool _proceedNext;

	uint32_t _tarePosition;

    SorterType _sorterType;
  };

}

#endif // SORTER_H
