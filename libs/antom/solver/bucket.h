/********************************************************************************************
bucket.h -- Copyright (c) 2017, Tobias Paxian

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

#ifndef ANTOM_BUCKET_H
#define ANTOM_BUCKET_H

#include <vector>

#include "sorter.h"
#include "softclause.h"

namespace antom
{
  class Antom;
  class Cascade;
  struct SoftClauseNodes;

  class Bucket
  {
  public:
    Bucket(Antom* antom, Cascade* cascade, uint16_t position = ~0);

    ~Bucket(void);

    /**
     * @brief size of bucket.
     * @param reduceByFactorNthOutput
     *         Important for calculating recurcive depth with n'th Bucket entries taken as subBuckets.
     * @return Sum of all SubBucket sizes + sorter.size.
     */
    uint32_t size(bool reduceByFactorNthOutput = false);

    /**
     * @brief AddSoftClauseNode to the bucket.
     * @param SCNode
     */
    void AddSoftClauseNode(SoftClauseNodes* SCNode);

    /**
     * @brief Adds a Tare to the bucket.
     * @param tare
     */
    void AddTare(uint32_t tare);

    /**
     * @brief EncodeTop
     *          Encodes all subBuckets and merges them into the bucket _sorter.
     * @param howOftenUsed
     *          Indicates how often the bucket at its top is used in other buckets.
     */
    void EncodeTop(uint32_t howOftenUsed = 0);

    /**
     * @brief EncodeTopAddAtLast
     *          Almost same as EncodeTop but merges all subbuckets always into _sorter -> worse balanced but less complicated.
     */
    void EncodeTopAddAtLast(void);

    /**
     * @brief CreateTotalizerEncodeTree
     *          Creates the tree for recursivly for all _subtrees and its sorter.
     */
    void CreateTotalizerEncodeTree(void);

    /**
     * @brief MergeSorterWith
     *          Merges into _sorter the given vector.
     * @param MergeVector
     */
    void MergeSorterWith(std::vector<uint32_t> MergeVector);

    /**
     * @brief SolveTares
     *          Solves all tares of this bucket.
     * @param diffEstimatedToCurrentSatWeight
     * @return Returns currentresult of solver->solve
     */
    uint32_t SolveTares(uint64_t diffEstimatedToCurrentSatWeight, uint32_t currentresult = ANTOM_UNKNOWN);

    /**
     * @brief GetEveryNthOutput of sorter output.
     * @param n
     * @return Every Nth output - if encoded.
     */
    std::vector<uint32_t> GetEveryNthOutput(uint32_t n);


    /**
     * @brief SetSolvingParameters
     *          Necessary for solving with multiple cascade to differenciate between modes.
     * @param assumptions
     * @param upperWeightBoundOfTareCascade
     * @param weightOfAllTareCascadeTares
     * @param groundWeight
     * @param SolvingTares
     */
    void SetSolvingParameters(std::vector<uint32_t> *assumptions,
                              uint64_t upperWeightBoundOfTareCascade,
                              uint64_t weightOfAllTareCascadeTares,
                              int64_t groundWeight,
                              MultipleCascadeSolvingState SolvingTares);

    /**
     * @brief SolveBucketReturnMaxPosition
     *          Solves all necessary positions.
     * @param onlyWithAssumptions
     *          Doesn't set unit clauses for results.
     * @param localCalc
     *          Only due to this cascade.
     * @return
     *          last solvable position.
     */
    uint32_t SolveBucketReturnMaxPosition(bool onlyWithAssumptions, bool localCalc);

    /**
     * @brief Dump current bucket and all subBuckets.
     * @param depth
     *          Current depth of bucket.
     * @return depth + maxdepth of subBuckets.
     */
    uint16_t DumpAndGetMaxDepth(uint16_t depth);

    // indicates how often this bucket is used as subBucket in other buckets.
    uint32_t _howOftenUsed;

    Sorter* _sorter;

    friend class Cascade;
    friend class MultipleCascade;

  private:
	// Copy constructor.
	Bucket (const Bucket&) = default;
	
    // Assignment operator.
	Bucket& operator = (const Bucket&) = default;

  
    Antom* _antom;
	Settings* _setting;
	Control* _control;
    Cascade* _cascade;

    // variables:
    std::vector< Bucket* > _subBuckets;
    std::vector< SoftClause* >* _softClauses;
    std::vector< uint32_t >* _outputs;
    std::vector< uint32_t > _tares;
    std::vector<uint32_t> _bucketAssumptions;

    uint32_t _tarePosition;
    uint16_t _bucketPosition;
    uint32_t _base;
    uint32_t _nthOutputTaken;

    uint32_t _topEntries;
    uint32_t _binaryTopCl;
    uint32_t _binaryTopClEstimated;
    uint32_t _binaryBottomCl;
    uint32_t _binaryBottomClEstimated;
    uint32_t _ternaryTopCl;
    uint32_t _ternaryTopClEstimated;
    uint32_t _ternaryBottomCl;
    uint32_t _ternaryBottomClEstimated;

    // base^position
    uint64_t _multiplicator;
    // needed for multiple cascade weight calculation as base / ground weight
    int64_t _groundWeight;

    uint32_t _localMaxPos;
    uint32_t _positionsCuttedAtBottom;

    // to calculate final boundaries, this accuracy is needed in the multiple main cascade
    uint64_t _weightBoundaryAccuracy;
    uint64_t _upperWeightBoundOfTareCascade;
    uint64_t _sumSoftWeightsOfTareCascade;
    MultipleCascadeSolvingState _solvingState;
    bool _encoded;
    bool _encodeTreeGenerated;
    bool _isLastBucket;

    /**
     * @brief EstimateNumberOfClauses
     *          Binary Clauses, Ternary Clauses, due to totalizer sorting network
     * @param top
     */
    void EstimateNumberOfClauses(bool top);

    /**
     * @brief CalculateNumberOfClauses
     *          Binary Clauses, Ternary Clauses of Top and Bottom Buckets.
     * @param top
     *          For Top or Bottom Buckets.
     * @param setCurrent
     *          Set Current Values as initial Values.
     * @param calcEstimation
     *          Additionally Calculates the same values as an estimation.
     */
    void CalculateNumberOfClauses(bool top, bool setCurrent, bool calcEstimation);

    /**
     * @brief DumpSolveInformation
     *          Used by SolveBucketReturnMaxPosition in every iteration to dump informations.
     * @param head
     *          First Call, special infos.
     * @param localCalc
     *          only this cascade
     * @param currentresult
     *          SAT, UNSAT, UNKNOWN.
     * @param lastPos
     *          Last succesfully solved.
     * @param actualPos
     *          Next to solve.
     */
    void DumpSolveInformation(bool head, bool localCalc, uint32_t currentresult = ANTOM_SAT, uint32_t lastPos = 0, uint32_t actualPos = 0);

    /**
     * @brief EvaluateResult
     *          Used by SolveBucketReturnMaxPosition after solving bucket. Different case handling of last call sat, unsat, unknown.
     * @param currentresult
     *          SAT, UNSAT, UNKNOWN
     * @param actualPos
     *          last treated position.
     * @param lastPos
     *          position before actualPos.
     * @param onlyWithAssumptions
     * @return last solvable position.
     */
    uint32_t EvaluateResult(uint32_t currentresult, uint32_t actualPos, uint32_t lastPos, bool onlyWithAssumptions);

    /**
     * @brief CalculateSatWeight
     *          Update local (cascade) and global (antom) sat weight, if global weight is better than before, update _lastModel.
     * @param localCalc
     * @return due to local Calc return Sat weight of cascade or antom.
     */
    uint64_t CalculateSatWeight(bool localCalc);

    /**
     * @brief SetAsUnitClause
     *          Sets actual Pos as unit clause or assumption, due to currentResult. If lastPos1, then does the necessary coding if UNSAT.
     * @param actualPos
     * @param currentResult
     * @param onlyWithAssumptions
     */
    void SetAsUnitClause(uint32_t actualPos, uint32_t currentResult, bool onlyWithAssumptions);

    /**
     * @brief CalcBoundaries
     *          Upper and Lower Bound, due to the given parameters and cascade type.
     * @param actualPos
     * @param localCalc
     * @param finalCalculation
     */
    void CalcBoundaries(uint32_t actualPos, bool localCalc, bool finalCalculation = false);

    /**
     * @brief CalcNextPositionToSolve
     *          Due to cascade type.
     * @param satWeight
     * @param lastPos
     * @return position
     */
    uint32_t CalcNextPositionToSolve(uint64_t satWeight, uint32_t lastPos);

    /**
     * @brief GetAssumptions
     *          Normal Cascade: - only last Pos negated.
     *          _mainCascade in multiple cascade mode: additional weights + lastPos negated
     *          solve only with assumptions: additional all solved positions.
     * @param actualPos
     * @return assumptionVector
     */
    std::vector<uint32_t> GetAssumptions(uint32_t actualPos);



    /**
     * @brief CutMaxPos Cut at max position
     * @return max position
     */
    uint32_t CutMaxPos(bool solve);
    uint32_t CutMinPos(bool solve);

    void SetMaxPos(uint32_t maxPos);
  };

}

#endif // ANTOM_BUCKET_H
