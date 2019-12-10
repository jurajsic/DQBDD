/********************************************************************************************
multiplecascade.h -- Copyright (c) 2017, Tobias Paxian

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

#ifndef MULTIPLECASCADE_H
#define MULTIPLECASCADE_H

#include "cascade.h"


#define MULTIPLECASCADE_MINDIV 5

namespace antom
{
  class Antom;
  class Bucket;
  class Control;
  //class Sorter;
  class SolverProxy;
  struct SoftClauseNodes;

class MultipleCascade
{

  public:
    MultipleCascade(Antom* antom, bool onlyByTares, bool solveTareCascadeOnlyByTares, uint32_t cascDiv, uint32_t maxBucketSize, uint32_t nOfCasc, InterimResult interimResult, uint64_t sumOfSoftWeights);

    ~MultipleCascade(void) {}

    /**
     * @brief DivideAndConnectCascades
     *          Divides softClauses into multiple cascades due to divideStrategy and maxSize, after that they are connected to solve them.
     * @param divideStrategy
     *          Sorted SCs
     *          Random SCs
     *          Separate by weight differences.
     *          ...
     * @param softClauses
     *          All SC's to process.
     * @return False    - if it is to be solved in normal cascade mode.
     *         True     - if it is Divided normally.
     */
    bool DivideAndConnectCascades(MultipleCascadeDivideStrategy divideStrategy, std::vector<SoftClause *> *softClauses);

    /**
     * @brief Solve
     *          Solves multiple cascade.
     * @return
     *          Returns result.
     */
    uint32_t Solve();

    /**
     * @brief GetAssumptions
     *          Assumptions of all Additional Tares.
     *          Needed for _mainCascade->SolveLastBucket()
     * @return
     *          Assumption vector.
     */
    std::vector<uint32_t> GetAssumptions();

    /**
     * @brief DumpFinalValues
     *          Dumps some Calculation Information.
     * @param verbosity
     *          Dumps most information only if value higher than verbosity.
     */
    void DumpFinalValues(uint16_t verbosity);

    friend class Cascade;
private:
	
	// Copy constructor.
    MultipleCascade (const MultipleCascade&) = default;

    // Assignment operator.
    MultipleCascade& operator = (const MultipleCascade&) = default;
	
    Antom* _antom;
	Settings* _setting;

    bool _onlyByTares;
    bool _solveTareCascadeOnlyByTares;
    uint32_t _base;
    uint64_t _satWeight;
    uint64_t _highestMultiplicator;
    std::vector<SoftClause*> _tareSoftClauses;
    std::vector<std::vector<SoftClause*>> _partedSoftClauses;
    std::vector<Cascade*> _cascades;
    std::multimap<uint64_t, Cascade*> _cascadesMM;
    Cascade* _mainCascade;
    Cascade* _tareCascade;
    uint64_t _tareWeight;
    uint64_t _tareWeightWithoutAddTs;
    uint64_t _tareWeightOfAddTs;
    uint32_t _maxMainCascadePosition;
    uint32_t _maxTareCascadePosition;
    uint32_t _cascadeDivider;
    uint32_t _maxBucketSize;
    uint32_t _numberOfCascades;
    int64_t _estimatedWeightBoundaries[2];
    int64_t _tareWeightUNSAT;
    int64_t _tareWeightSAT;
    InterimResult _interimResult;
    uint64_t _sumOfSoftWeights;

    bool _solveTareCascade;
    bool _solveLikeNormalCascade;
    std::vector<uint32_t> _tareAssumptions;

    /**
     * @brief CreateSCVectorFromCascade
     *          Creates SC-Vector from cascades to create _tareCascade from it.
     * @param cascades
     *          Cascades to get tares from.
     */
    void CreateSCVectorFromCascades(std::vector<Cascade*> cascades);

    /**
     * @brief CreateSCesFromTares
     *          Creates SC-Vector from given Tares to create _tareCascade from it.
     * @param tares
     *          Given Tares to create vector from.
     * @param base
     *          Overall Base to calc weights.
     */
    void CreateSCVectorFromTares(std::vector<std::vector<uint32_t> > *tares, uint32_t base);

    /**
     * @brief CreateSCForTare
     *          Creates a SC for given values. Creates a new binary clause to change direction of tare value.
     * @param tare
     *          Tare to process.
     * @param weight
     *          Weight to process.
     * @return
     *          Created Softclause.
     */
    SoftClause* CreateSCForTare(uint32_t tare, uint64_t weight);

    /**
     * @brief DivideSoftClauseVector
     *          Saves them into _partedSoftClauses due to divideStrategy.
     * @param divideStrategy
     * @param softClauses
     *          SoftClauseVecotr
     * @return
     *      True - if divided.
     *      False - if not divided - solve in single cascade mode.
     */
    bool DivideSoftClauseVector(MultipleCascadeDivideStrategy divideStrategy, std::vector<SoftClause *> *softClauses);

    /**
     * @brief FillSubCascades
     *          Creates sub Cascades due to _partedSoftClauses and saves them into _cascadesMM MultiMap.
     * @param partedSoftClauses
     *
     */
    void FillSubCascades(std::vector<std::vector<SoftClause *> > partedSoftClauses);

    /**
     * @brief ConnectSubCascades
     *          Only if not SolveByTares is set. Connects the given cascades and assures no Bucket gets bigger than maxSize.
     *          Saves final cascade in _mainCascade.
     */
    void ConnectSubCascades();

    /**
     * @brief AddSolutionAsUnitClauses
     *          For debugging, fast way to add a vector as unit clauses.
     */
    void AddSolutionAsUnitClauses();
    void ConnectLeastWeightedSmallestCascades();
    void ConnectHighestWeightedCascades();
};

}

#endif // MULTIPLECASCADE_H
