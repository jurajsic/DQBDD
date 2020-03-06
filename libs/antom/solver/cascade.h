/********************************************************************************************
cascade.h -- Copyright (c) 2017, Tobias Paxian

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

#ifndef CASCADE_H
#define CASCADE_H


#include "sorter.h"
#include <map>

namespace antom
{
  class Antom;
  class Bucket;
  class Control;
  class MultipleCascade;
  class SolverProxy;
  struct SoftClauseNodes;

class Cascade
{

  public:
  Cascade(Antom* antom, MultipleCascade* multipleCascade = nullptr, bool onlyByTares = false);

    ~Cascade(void) {}

    /**
     * @brief Fills the bucket structure due to given multiple strategy
     * @param softClauses
     *          Which softclauses to use for filling.
     */
    void Fill(std::vector<SoftClause*>* softClauses, PartitionStrategy partitionStrategy, EncodeStrategy encodeStrategy);

    /**
     * @brief FillStructure
     *          To call directly from Multiplecascade, this part of Fill is extracted!
     */
    void FillStructure(PartitionStrategy partitionStrategy, EncodeStrategy encodeStrategy);

    /**
     * @brief Encodes the buckets, due to given strategy.
     * @return True if all worked, else false.
     */
    bool Encode();

    /**
     * @brief Solves the cascade.
     *          Solves at first last bucket, then the tares.
     * @return ANTOM_SAT / UNSAT / UNKNOWN
     */
    uint32_t Solve();

    /**
     * @brief DumpSCNodeStructure
     *          Dumps the structure of the SoftClauseNodes
     * @param dumpingSCTree
     *          Which SCTree to dump
     * @param verbosity
     *          Only dumps out if verbosity at least equal to that value.
     */
    void DumpSCNodeStructure(std::vector<SoftClauseNodes*>* dumpingSCTree, uint16_t verbosity = 2);

    /**
     * @brief DumpBucketStructure
     *          If buckets are filled with values, then this function dumps the structure of the BucketTree.
     * @param onlyLastBucket
     *          Only dumps out Last Bucket.
     * @param verbosity
     *          Only dumps out if verbosity at least equal to that value.
     */
    void DumpBucketStructure(bool onlyLastBucket = false, uint16_t verbosity = 2);

    /**
     * @brief CreateSoftClauseTree
     *          Creates one or two SoftClauseTrees filled with SoftClauseNodes.
     * @param softClauses
     *          All Softclauses to be processed.
     * @param split
     *          Shall the values directly splittet into two SoftClauseTrees.
     */
    void CreateSoftClauseTree(std::vector<SoftClause*>* softClauses, bool split = true);

    uint64_t GetHighestMultiplicator();

    uint32_t GetMaxBucketSize();

    friend class Bucket;
    friend class MultipleCascade;


private:
	// Copy constructor.
    Cascade (const Cascade&) = default;

    // Assignment operator.
    Cascade& operator = (const Cascade&) = default;
	
    Antom* _antom;
    Control* _control;
	Settings* _setting;
    MultipleCascade* _multipleCascade;
	uint32_t _base;
    bool _onlyByTares;

    uint64_t _satWeight;
    uint64_t _tareWeight;
    int64_t _estimatedWeightBoundaries[2];
    int64_t _weightToSubstract;
    uint64_t _sumOfSoftWeights;
    bool _softClauseTreeCreated;

    std::vector<uint32_t> _collectedCascadeAssumptions;
    std::vector< Bucket* > _structure;
    uint32_t _numberOfBuckets;
    uint32_t _totalBucketEntriesperWeight;
    uint32_t _totalBucketOccurrences;
    std::vector <uint32_t> _totalBucketEntries;

    uint32_t _maxSorterDepth;
    uint64_t _highestBucketMultiplicator;
    uint64_t _upperWeightBoundAllLowerCascades;
    uint32_t _howManyPositionsCuttedAtBottom;

    std::vector<SoftClause*> _softClauses;
    std::vector< SoftClauseNodes* > _softClauseTree;
    std::vector< SoftClauseNodes* > _processingSoftClauseTree;
    // contains SoftClauseNodes having no match partner at the moment
    // but there is a chance that a match can be found later again.
    std::multimap<uint32_t, SoftClauseNodes*> _processingPercentOffTree;
    uint32_t _howOftenReinsertedFromProcessingPercentOffTree;



    // Functions

    std::vector<std::vector< uint32_t >> GetTareVector();

    /**
     * @brief PartitionSoftClauseTree
     *          Splits the tmpSoftClauseTree into processingSoftClauseTree and softClauseTree
     * @param tmpSoftClausesTree
     *          SCT to split.
     */
    void PartitionSoftClauseTree(std::vector<SoftClauseNodes*>* tmpSoftClausesTree);

    /**
     * @brief PartitionSoftClauses
     *          Due to the given partitionStrategy.
     * @param partitionStrategy
     */
    void PartitionSoftClauses(PartitionStrategy partitionStrategy);

    /**
     * @brief GroupByWeight
     *          Groups all SCN by their weight.
     */
    void GroupByWeight();

    /**
     * @brief GroupByBiggestRepeatingEntry
     *          Groups all SCN always by the given maxNodeIndex due to given heuristic.
     */
    void GroupByBiggestRepeatingEntry();

    /**
     * @brief GetMaxNodeIndex
     * @return Index of SoftClauseTree with most Bucket entries!
     *         If 2 SCN have same # entries, return the index of the one with more occurrences.
     */
    uint16_t GetMaxNodeIndex();

    /**
     * @brief GetBenefitOfMergingNodes
     *          Calculates due to a given heuristic the Benefit of merging two nodes.
     * @param NodesToMerge
     *          A list of currently two nodes to merge
     * @param tmpBucketIndices
     *          Indices of all Buckets - NodesToMerge is in.
     * @param dump
     *          With or without dumping the values.
     * @param heuristic
     *          Due to which heuristic the values shall be chosen.
     * @param percentOff
     *          Max percentage the size of minSCN/maxSCN should be differ.
     *
     * @return benefit  due to heuristic.
     *          0       if no merge is possible.
     *          -size   of SCN if no merge is possible but a merge is maybe later on possible.
     *
     */
    int32_t GetBenefitOfMergingNodes(std::vector< SoftClauseNodes* > NodesToMerge, std::vector<uint16_t> *tmpBucketIndices, bool dump, uint32_t heuristic = 0, int32_t maxPercentOff = 100);

    /**
     * @brief MergeNodes
     *          Merges the two nodes of given indices.
     * @param tmpBucketIndices
     *          All Indices maxBucketIndex has a weight in.
     * @param nodeIndexToMergeWith
     *          If nodeIndexToMergeWith == -1 there is no mo merge possible with maxBucketIndex - move the SCN to the final _softClauseTree.
     *          If nodeIndexToMergeWith == -2 no merge possible due to percentOff, but a later merge can be possible - move the SCN to the _processingPercentOffTree.
     * @param maxBucketIndex
     *          Node wich is in most buckets in.
     */
    void MergeNodes(std::vector<uint16_t> *tmpBucketIndices, int32_t nodeIndexToMergeWith, int32_t maxBucketIndex);

    /**
     * @brief CalculateNumberOfPossibleSubmerges
     *          Used for Heuristic Calculations.
     * @param tmpBucketIndices
     *          All Indices in which the previously merged Bucket has a value in.
     * @param NodesToMerge
     *          NodeList of "not to use again" nodes in calculation.
     * @param depth
     *          If depth > 0, it calculates the depth until maxdepth 2 and returns depth + possible submerges at that depth.
     * @return #possibleSubmerges                   if no depth is given or depth=0.
     *         Depth(max 2) + #possibleSubmerges    if depth is given - can be very time consuming depending on number and weights of given cascade softclauses.
     */
    int32_t CalculateNumberOfPossibleSubmerges(std::vector<uint16_t> *tmpBucketIndices, std::vector< SoftClauseNodes* > NodesToMerge, int32_t depth = 0);

    /**
     * @brief CalculateNodeIndexToMergeWith
     * @param maxBucketIndex
     *          Index of Bucket to process.
     * @param tmpBucketIndices
     *          All Indices in which maxBucketIndex has a value in.
     * @return Index    of _processingSoftClauseTree to merge maxBucketIndexWith.
     *            -1    if no match was found.
     *            -2    if no match was found but maxBucketIndex has to be inserted to the _processingPercentOffTree.
     */
    int32_t CalculateNodeIndexToMergeWith(uint16_t maxBucketIndex, std::vector<uint16_t> *tmpBucketIndices);

    /**
     * @brief FillBuckets
     *          Distribute the SoftClauseNodes to the buckets.
     */
    void FillBuckets();

    /**
     * @brief AddTaresToBuckets
     *          Every TopBucket, but the last, gets _base-1 many variables extra.
     */
    void AddTaresToBuckets(void);

    /**
     * @brief EncodeTopBuckets
     *          Due to the Cascade Structure, encodes all subBuckets and merges them into a single topBucket.
     */
    void EncodeTopBuckets(void);

    /**
     * @brief EncodeBottomBuckets
     *          Due to the Cascade Structure, encodes the union of every _base' entry of BottomBucket[i-1] with TopBucket[i].
     */
    void EncodeBottomBuckets(void);

    /**
     * @brief CreateTotalizerEncodeTree
     *          Creates structure of coding - to encode Buckets later on only if needed.
     */
    void CreateTotalizerEncodeTree(void);

    /**
     * @brief CalculateTotalBucketEntries
     *          Calculates the variables _totalBucketEntriesperWeight, _totalBucketOccurrences and the vector _totalBucketEntries.
     * @param SCTree - which SCTree to process.
     * @param add - add the weight to the existing values .
     *          Important to calculate the total weight at the beginning if vars in two Trees.
     */
    void CalculateTotalBucketEntries(std::vector<SoftClauseNodes*>* SCTree, bool add);

    /**
     * @brief AtLeastTwoBucketsInCommon
     * @param SCN1 <- first SoftClauseNode
     * @param SCN2 <- second SoftClauseNode
     * @return true if at least two buckets in common.
     */
    bool AtLeastTwoBucketsInCommon(SoftClauseNodes* SCN1, SoftClauseNodes* SCN2);

    /**
     * @brief SolveTares
     *          Solves the tares from the second last bucket downwards.
     * @return ANTOM_SAT / UNSAT / UNKNOWN
     */
    uint32_t SolveTares(bool solveTareInLastBucketToo = false);

    /**
     * @brief Calculates Max and Average Entries of the actual buckets and prints them.
     */
    void CalculateBucketEntries(void);

    /**
     * @brief DumpMaxNodeOverlappingsAndHeuristicValues
     *          Almost same as CalculateNodeIndexToMergeWith - but dumping the values.
     * @param maxBucketIndex
     *          Index of SCN to process.
     * @param tmpBucketIndices
     *          All Indices in which maxBucketIndex has a value in.
     */
    void DumpMaxNodeOverlappingsAndHeuristicValues(uint16_t maxBucketIndex, std::vector<uint16_t> *tmpBucketIndices);

    /**
     * @brief DumpModelOfTares
     *          All tares and their last assigned value.
     * @param verbosity
     *          Only dumps out if verbosity at least equal to that value.
     */
    void DumpModelOfTares(uint16_t verbosity = 2);

    /**
     * @brief DumpBucketSolveInformation
     *          Weight Boundaries, sat weights, actual position. Called from Bucket during bucketSolving.
     * @param actualPos
     *          Actual position of solving that bucket.
     * @param _isLastBucket
     *          If it is the last bucket, it outputs additionally the weight boundaries.
     * @param verbosity
     *          Only dumps out if verbosity at least equal to that value.
     */
    void DumpBucketSolveInformation(uint32_t actualPos, bool _isLastBucket, uint16_t verbosity = 3);

    /**
     * @brief CountSatisfiedSoftClauses
     *          Returns number of satisfied softclauses of this bucket and save _satWeight.
     * @param bucket
     *          If bucket is given, it calculates the number of satisfied softclauses of that bucket.
     * @param model
     *          Which model to use for calculation
     * @return Number/Weight of satisfied softclauses.
     */
    uint64_t CountSatisfiedSoftClauses(Bucket* bucket, const std::vector<uint32_t>& model);

    /**
     * @brief CountSatisfiedSoftClauses
     *          Returns number of satisfied softclauses of this bucket.
     * @param bucket
     *          If bucket is given, it calculates the number of satisfied softclauses of that bucket.
     *          If bucket isn't given, it calculates the number of all satisfied softclauses.
     * @param model
     *          Which model to use for calculation
     * @param addWeight
     *          If addWeight is given it calculates the sum of weights of satisfied softclauses.
     * @return Number/Weight of satisfied softclauses.
     */
    uint64_t CountSatisfiedSoftClauses(std::vector<SoftClause *> softclauses, const std::vector<uint32_t> &model, bool addWeight);

    /**
     * @brief DumpNumberOfBucketEntriesOrClauses
     *          Used for dumping several kinds of structure information due to given values.
     * @param top
     * @param bottom
     * @param estimated
     * @param calculated
     * @param binary
     * @param ternary
     */
    void DumpNumberOfBucketEntriesOrClauses(bool top, bool bottom = false, bool estimated = false, bool calculated = false, bool binary = false, bool ternary = false);

    /**
     * @brief DumpNumberOfBucketsAndClauses
     *          Calls 11 times DumpNumberOfBucketEntriesOrClauses with different parameters to dump.
     */
    void DumpNumberOfBucketsAndClauses();

    /**
     * @brief UnionBucketsIntoLast
     *          If divideStrategy==1 put all buckets together to one final bucket.
     */
    void UnionBucketsIntoLast();

    /**
     * @brief CountSumOfSoftWeights
     *          Counts all weights of softClauses.
     * @param softClauses
     */
    void CountSumOfSoftWeights(std::vector<SoftClause *> *softClauses);


    /**
     * @brief AddNewBucketsTillSizeBoundary
     *          If an additional bucket is added, the size is halfed and multiplicator is doubled. Calls AddNewBucketsTillMultiplicatorMatches
     *          such that no code is duplicated.
     * @param maxSize
     *          boundary
     * @param onlySolveWithTares
     *          if true, it continues adding until maxSize of bucket is 2
     * @param addTareToLastBucket
     * @return True if boundary could be matched.
     */
    bool AddNewBucketsTillSizeBoundary(uint32_t maxSize, bool onlySolveWithTares = false, bool addTareToLastBucket = false);

    /**
     * @brief AddNewBucketsTillMultiplicatorMatches
     *          If an additional bucket is added, the size is halfed and multiplicator is doubled.
     * @param maxMultiplicator
     *          boundary
     * @param onlySolveWithTares
     *          if true, it continues adding until maxSize of bucket is 2.
     * @param addTareToLastBucket
     * @return True if multiplicator could be reached.
     */
    bool AddNewBucketsTillMultiplicatorMatches(uint64_t maxMultiplicator, bool onlySolveWithTares = false, bool addTareToLastBucket = false);

    /**
     * @brief AddAdditionalBucket
     *          Only use if unioned into last bucket - add bucket and additional tare
     */
    void AddAdditionalBucket();

    /**
     * @brief AddTare
     *          adds tare to bucket with position.
     * @param position
     */
    void AddTare(unsigned long position);

    /**
     * @brief AddAsManyBucketsAsPossible
     *          Calls AddNewBucketsTillMultiplicatorMatches wit max uint64_t number, to add as many buckets as possible.
     *
     * ATTENTION:   IF SUM OF SOFTCLAUSE WEIGHTS IS CLOSE TO 2^64 THERE CAN BE AN OVERFLOW!!!
     *              TO AVOID THIS USE UINT128_T FOR ALL WEIGHT RELATED STUFF!
     */
    void AddAsManyBucketsAsPossible();

    /**
     * @brief SolveAllTares
     *          solves all tares, doesn't calculate if one has to be fulfilled already. For debug reasons.
     * @return
     */
    uint32_t SolveAllTares();

    std::vector<uint32_t> CalculateAssumptionsFor(int64_t weight, int32_t startingPos);
    uint32_t SolveTareWeightPlusOne();
    int32_t SetUnitClauses(int32_t startingPos);

    /**
     * @brief CutMaxPos Cut at local max position of last bucket
     * @return last solvable position
     */
    uint32_t CutMaxPos(bool solve = true);
    uint32_t CutMinPos(bool solve = true);
};

}

#endif // CASCADE_H
