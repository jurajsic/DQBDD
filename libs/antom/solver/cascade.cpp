/********************************************************************************************
cascade.cpp -- Copyright (c) 2017, Tobias Paxian

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
#include <cmath>
#include <iomanip>
#include <numeric>
#include <algorithm>

// Include antom related headers.
#include "cascade.h"
#include "multiplecascade.h"
#include "bucket.h"
#include "softclausenodes.h"
#include "antom.h"
#include "control.h"
#include "totalizerencodetree.h"
#include "timemeasurement.h"
#include "timevariables.h"

namespace antom
{

// Constructor
  Cascade::Cascade(Antom* antom, MultipleCascade* multipleCascade, bool onlyByTares):
    _antom(antom),
    _control(antom->_control),
	_setting(antom->_antomSetting),
    _multipleCascade(multipleCascade),
    _base(_setting->base),
    _onlyByTares(onlyByTares),
    _satWeight(0),
    _tareWeight(0),
    _weightToSubstract(0),
    _sumOfSoftWeights(0),
    _softClauseTreeCreated(false),
	_collectedCascadeAssumptions(),
	_structure(),
    _numberOfBuckets(0),
	_totalBucketEntriesperWeight(0),
	_totalBucketOccurrences(0),
	_totalBucketEntries(),
    _maxSorterDepth(0),
	_highestBucketMultiplicator(0),
    _upperWeightBoundAllLowerCascades(0),
    _howManyPositionsCuttedAtBottom(0),
    _softClauses(),
	_softClauseTree(),
	_processingSoftClauseTree(),
	_processingPercentOffTree(),
    _howOftenReinsertedFromProcessingPercentOffTree(0)

{
    assert(antom != NULL);
}

void Cascade::Fill(std::vector<SoftClause*>* softClauses, PartitionStrategy partitionStrategy, EncodeStrategy encodeStrategy)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    if (_onlyByTares)
        encodeStrategy = ENCODEONLYIFNEEDED;

    TimeMeasurement timeFillingBuckets(&_antom->_timeVariables->fillingBuckets, true);

    CountSumOfSoftWeights( softClauses);

    FillStructure(partitionStrategy, encodeStrategy);

    if (_setting->verbosity > 3)
        std::cout << std::endl << "Buckets are filled - structure is dumped" << std::endl << std::endl;
}

void Cascade::CountSumOfSoftWeights( std::vector<SoftClause*>* softClauses )
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    _softClauses = *softClauses;
    _sumOfSoftWeights = 0;

    // sum over all SC weights of CASCADE!
    std::for_each(_softClauses.begin(), _softClauses.end(), [&](SoftClause* SC) { _sumOfSoftWeights += SC->weight; });
    _satWeight = CountSatisfiedSoftClauses(NULL, _antom->_lastModel);

    if (_setting->verbosity > 4)
    {
        std::cout << "c SAT Weight of Cascade..: " << _satWeight << std::endl;
        std::cout << "c SAT Weight of Antom....: " << _antom->_satWeight << std::endl;
        std::cout << "c softWeights of Cascade.: " << _sumOfSoftWeights << std::endl;
        std::cout << "c softWeights of Antom...: " << _antom->_sumOfSoftWeights << std::endl;
    }
}

uint32_t Cascade::CutMaxPos(bool solve)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    _structure.back()->_isLastBucket = true;
    return _structure.back()->CutMaxPos(solve);
}

uint32_t Cascade::CutMinPos(bool solve)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    _structure.back()->_isLastBucket = true;
    return _structure.back()->CutMinPos(solve);
}



void Cascade::FillStructure(PartitionStrategy partitionStrategy, EncodeStrategy encodeStrategy)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    if (_onlyByTares)
        encodeStrategy = ENCODEONLYIFNEEDED;

    PartitionSoftClauses(partitionStrategy);
    FillBuckets();

    AddTaresToBuckets();

    if (encodeStrategy == ENCODEONLYIFNEEDED)
    {
        UnionBucketsIntoLast();
        if(_onlyByTares)
            AddAsManyBucketsAsPossible();
        DumpBucketStructure(true, 3);

    } else
    {
        DumpBucketStructure(false, 3);
    }
}

void Cascade::AddAsManyBucketsAsPossible()
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

//    std::cout << "size: " << _structure.back()->size() << std::endl;
//    std::cout << "_structure.back()->_tares[0]: " << _structure.back()->_tares[0] << std::endl;

    if (_structure.back()->_tares.empty())
        AddTare(_structure.size()-1);

    if(_setting->interimResult==CUTATTOP)
    {
        CutMaxPos();
    }
    _structure.back()->_encodeTreeGenerated = false;

    /** ATTENTION NOT ALWAYS GIVEN
     *
     * If e.g. highest number is close to uint64_t...
     */
    if(AddNewBucketsTillMultiplicatorMatches(static_cast<uint64_t>(-1), true))
    {
        std::cout << "EXCEPTION - cannot be solved only by tares because of multiplicator limit of 64 Bit!";
        _onlyByTares = false;
        return;
    }
    DumpBucketStructure(true);

    _estimatedWeightBoundaries[0] = -_highestBucketMultiplicator + 1;
    _estimatedWeightBoundaries[1] = _highestBucketMultiplicator;

//    std::cout << "EWB0: " << _estimatedWeightBoundaries[0] << std::endl;
//    std::cout << "EWB1: " << _estimatedWeightBoundaries[1] << std::endl;


}

bool Cascade::Encode()
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    switch(_setting->encodeStrategy)
    {
    case ENCODEALL:
        EncodeTopBuckets();
        if ( _control->ReachedLimits() )
            return false;

        EncodeBottomBuckets();
        CalculateBucketEntries();
        DumpBucketStructure(false, 4);
        break;
    case ENCODEONLYIFNEEDED:
        if (_onlyByTares)
            return true;
        CreateTotalizerEncodeTree();
        CalculateBucketEntries();
        DumpBucketStructure(true, 4);
        break;
    }
    if ( _control->ReachedLimits() )
        return false;

    if (_setting->verbosity < 4)
        return true;

    std::cout << std::endl << "Buckets are encoded - structure is dumped" << std::endl << std::endl;
    return true;
}

uint32_t Cascade::Solve()
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    uint32_t currentresult(1);

    if (_antom->_satWeight == _antom->_sumOfSoftWeights)
        return ANTOM_SAT;

    //TimeMeasurement TimeSolvingLastBucket(&_antom->_timeVariables->solvingLastBucket);
    if(!_onlyByTares)
    {
        _structure.back()->SolveBucketReturnMaxPosition(false, false);
        if ( _antom->_resultUnknown )
        {
            return ANTOM_UNKNOWN;
        }
        if (_setting->encodeStrategy == ENCODEONLYIFNEEDED && _setting->createGraphFile != "")
            _structure.back()->_sorter->_outputTree->DumpOutputTree(_setting->createGraphFile + "_withOutputs.tgf", true);
    }

    if (_antom->_satWeight == _antom->_sumOfSoftWeights)
        return ANTOM_SAT;

    //TimeSolvingLastBucket.~TimeMeasurement();
    currentresult = SolveTares();
//    currentresult = SolveAllTares();

    return currentresult;
}

void Cascade::CreateSoftClauseTree(std::vector<SoftClause*>* softClauses, bool split)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    _processingSoftClauseTree.clear();
    _softClauseTree.clear();


    //sort without changing order of elements
//    std::stable_sort(softClauses->begin(), softClauses->end(), SoftClause::bigger );
//    _numberOfBuckets = static_cast<uint32_t>(floor(log2(_softClauses[0]->weight)/log2(_base)));

    if (_setting->featureTest)
    {
        //sort with changing order of elements (better worst case runtime)
        std::sort(softClauses->begin(), softClauses->end(), SoftClause::bigger );
        _numberOfBuckets = static_cast<uint32_t>(floor(log2(_softClauses[0]->weight)/log2(_base)));
    } else
    {
        uint64_t maxValue = 0;
        for (auto softclause : *softClauses)
        {
            maxValue = softclause->weight > maxValue ? softclause->weight : maxValue;
        }
        _numberOfBuckets = static_cast<uint32_t>(floor(log2(maxValue)/log2(_base)));
    }

    // SoftClauseNode Structure is created. Extract function!
    for (uint32_t i = 0; i != softClauses->size(); ++i)
    {
        SoftClauseNodes* softClauseNode = new SoftClauseNodes((*softClauses)[i], _base);
        // softclause is only in one bucket
        if (split && softClauseNode->inHowManyBuckets != 1)
        {
            _processingSoftClauseTree.push_back( softClauseNode );
        } else
        {
            _softClauseTree.push_back( softClauseNode );
        }
    }
    if (split)
        _softClauseTreeCreated = true;
}

void Cascade::PartitionSoftClauseTree(std::vector<SoftClauseNodes*>* tmpSoftClausesTree)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    _processingSoftClauseTree.clear();
    _softClauseTree.clear();

    // SoftClauseNode Structure is created. Extract function!
    for (uint32_t i = 0; i != tmpSoftClausesTree->size(); ++i)
    {
        //SoftClauseNodes* softClauseNode = new SoftClauseNodes(_softClauses[i], _base);
        // softclause is only in one bucket

        if ((*tmpSoftClausesTree)[i]->inHowManyBuckets != 1)
        {
            _processingSoftClauseTree.push_back( (*tmpSoftClausesTree)[i] );
        } else
        {
            _softClauseTree.push_back( (*tmpSoftClausesTree)[i] );
        }
    }
    _softClauseTreeCreated = true;
}

std::vector<std::vector< uint32_t >> Cascade::GetTareVector()
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    std::vector<std::vector< uint32_t >> tares;

    uint32_t addNoTareToLastBucket = ( _setting->cascadeDivider > 0 ) ? 0 : 1;
    if (_setting->verbosity > 3)
        std::cout << "addNoTareToLastBucket: " << addNoTareToLastBucket  << std::endl;
    for (uint32_t ind = 0; ind < _structure.size() - addNoTareToLastBucket; ind++)
    {
        tares.push_back(_structure[ind]->_tares);
    }
    if (_setting->verbosity > 3)
        std::cout << "returning tares" << std::endl;
    return tares;
}

void Cascade::PartitionSoftClauses(PartitionStrategy partitionStrategy)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    if (!_softClauseTreeCreated)
        CreateSoftClauseTree(&_softClauses, true);

    switch(partitionStrategy)
    {
    case NOPARTITION:
        _softClauseTree.insert( _softClauseTree.end(), _processingSoftClauseTree.begin(), _processingSoftClauseTree.end() );
        _processingSoftClauseTree.clear();
        break;
    // both cases have the same grouping algorithm, but differ in the way they connect the weights.
    case GROUPBYWEIGHTADDATLAST:
    case GROUPBYWEIGHT:
        GroupByWeight();
        _softClauseTree.insert( _softClauseTree.end(), _processingSoftClauseTree.begin(), _processingSoftClauseTree.end() );
        _processingSoftClauseTree.clear();
        break;
    case GROUPBYBIGGESTREPEATINGENTRY:
        GroupByWeight();

        // actualize values as sum of both trees.
        CalculateTotalBucketEntries( &_processingSoftClauseTree, false);
        CalculateTotalBucketEntries( &_softClauseTree, true);

        DumpSCNodeStructure( &_processingSoftClauseTree, 2 );

        GroupByBiggestRepeatingEntry();

        assert(_processingSoftClauseTree.empty());
        break;
    }

    CalculateTotalBucketEntries( &_softClauseTree, false);

    DumpSCNodeStructure( &_softClauseTree, 2 );
}

void Cascade::GroupByBiggestRepeatingEntry()
{
    if (_processingSoftClauseTree.empty())
    {
        return;
    }
    //std::cout << std::endl << __func__ << std::endl;

    // get indices of sorted Bucket entries
    // generate Indice Vector to sort that vector!
    std::vector<uint16_t> sortedBucketIndices(_processingSoftClauseTree[0]->highestBucket + 1);
    std::size_t n(0);
    std::generate(std::begin(sortedBucketIndices), std::begin(sortedBucketIndices) + _processingSoftClauseTree[0]->highestBucket + 1, [&]{ return n++; });

    // sort Buckets by total entries.
    std::sort(  std::begin(sortedBucketIndices),
                std::end(sortedBucketIndices),
                [&](std::size_t i1, std::size_t i2) { return (_totalBucketEntries[i1] > _totalBucketEntries[i2] ); } );

    uint32_t ind(0);
    while (true)
    {
        if (_processingSoftClauseTree.size() == 1)
        {
            _softClauseTree.push_back(_processingSoftClauseTree.back());
            _processingSoftClauseTree.clear();
            std::cout << "c grouping iterations....: " << ind << std::endl;
            break;
        }

        uint16_t maxNodeIndex = GetMaxNodeIndex();

        assert(maxNodeIndex < _processingSoftClauseTree.size());

        // form a new vector with all indices of buckets containing a max Bucket entry
        std::vector<uint16_t> tmpBucketIndices;
        for (auto v : sortedBucketIndices)
        {
            if (v <= _processingSoftClauseTree[maxNodeIndex]->highestBucket && _processingSoftClauseTree[maxNodeIndex]->occursHowOftenInBucket[v] > 0)
                tmpBucketIndices.push_back(v);
        }

        if ((ind % 500 == 0 && _setting->verbosity > 2) || (ind % 250 == 0 && _setting->verbosity > 3) || _setting->verbosity > 4)
            DumpMaxNodeOverlappingsAndHeuristicValues(maxNodeIndex, &tmpBucketIndices);

        //TOBI: Is there another node with exactly the same overlapping? -> then merge more nodes!
        // Maybe there is even a bigger Set in the same subset to merge first with --> see bwt3cc
        int32_t nodeIndexToMergeWith = CalculateNodeIndexToMergeWith(maxNodeIndex, &tmpBucketIndices);

        MergeNodes(&tmpBucketIndices, nodeIndexToMergeWith, maxNodeIndex);

        if ( _setting->equalWeight > 0 && (ind % _setting->equalWeight == 0))
             GroupByWeight();

        if ((ind % 1000 == 0 && _setting->verbosity > 3)  || _setting->verbosity > 4)
        {
            std::cout << "ProcessingSoftClauseTree" << std::endl;
            DumpSCNodeStructure(&_processingSoftClauseTree, 3);

            std::cout << "FinalSoftClauseTree" << std::endl;
            DumpSCNodeStructure(&_softClauseTree, 3);

            std::cout << "_processingSoftClauseTree.size(): " << _processingSoftClauseTree.size() << std::endl;
        }

        ind++;
    }
    for (std::map<uint32_t, SoftClauseNodes*>::iterator it = _processingPercentOffTree.begin(); it != _processingPercentOffTree.end(); ++it)
    {
        _softClauseTree.push_back(it->second);
    }
    DumpSCNodeStructure(&_processingSoftClauseTree, 3);
    std::cout << "c percentOff reinsertions: " << _howOftenReinsertedFromProcessingPercentOffTree << std::endl;
}

void Cascade::GroupByWeight()
{
    std::vector<SoftClauseNodes*> tmpSCTree;
    // sort Buckets by total entries.
    std::sort(  std::begin(_processingSoftClauseTree),
                std::end(_processingSoftClauseTree),
                [&](SoftClauseNodes* sCN1, SoftClauseNodes* sCN2) { return (sCN1->weight > sCN2->weight ); } );

    for (uint32_t i = 0; i < _processingSoftClauseTree.size(); ++i)
    {
        if (i + 1 == _processingSoftClauseTree.size() || _processingSoftClauseTree[i + 1]->weight != _processingSoftClauseTree[i]->weight)
        {
            tmpSCTree.push_back( _processingSoftClauseTree[i] );
            continue;
        }
        Bucket* tmpBucket = new Bucket(_antom, this);
        while(i + 1 < _processingSoftClauseTree.size() && _processingSoftClauseTree[i + 1]->weight == _processingSoftClauseTree[i]->weight)
        {
            tmpBucket->AddSoftClauseNode(_processingSoftClauseTree[i]);
            i++;
        }
        tmpBucket->AddSoftClauseNode(_processingSoftClauseTree[i]);
        SoftClauseNodes* sCNode = new SoftClauseNodes(tmpBucket, _processingSoftClauseTree[i]->weight, _base);
        tmpSCTree.push_back( sCNode );
    }
    _processingSoftClauseTree = tmpSCTree;
}

uint16_t Cascade::GetMaxNodeIndex()
{
    // generate Indice Vector to get highest index from!
    std::vector<uint16_t> maxSCTIndices(_processingSoftClauseTree.size());
    std::size_t m(0);
    std::generate(std::begin(maxSCTIndices), std::begin(maxSCTIndices) + _processingSoftClauseTree.size(), [&]{ return m++; });

    // Get index of SoftClauseTree with most Bucket entries!
    // If 2 SCN have same # entries, take the one with more occurrences.
    return *std::max_element(std::begin(maxSCTIndices), std::end(maxSCTIndices), [&](std::size_t i1, std::size_t i2) {
        return ((_processingSoftClauseTree[i1]->inHowManyBuckets < _processingSoftClauseTree[i2]->inHowManyBuckets) ||
                 ((_processingSoftClauseTree[i1]->inHowManyBuckets == _processingSoftClauseTree[i2]->inHowManyBuckets) &&
                  (_processingSoftClauseTree[i1]->GetOccurrences() * _processingSoftClauseTree[i1]->size() < _processingSoftClauseTree[i2]->GetOccurrences() * _processingSoftClauseTree[i2]->size()))); } );

    // First Idea - sort whole tree...
    // sort whole SoftClauseTree by entries and if they are equal by occurences.
    // maybe to look later on only to the top elements. Otherwise it is sufficient to find in O(n) the biggest element.
    //std::sort(_processingSoftClauseTree.begin(), _processingSoftClauseTree.end(), SoftClauseNodes::SortByEntriesOccurrences);
}


void Cascade::DumpMaxNodeOverlappingsAndHeuristicValues(uint16_t maxNodeIndex, std::vector<uint16_t> *tmpBucketIndices)
{
    if (_setting->verbosity < 1)
        return;

    std::cout << "_processingSoftClauseTree.size: " << _processingSoftClauseTree.size() << std::endl;
    std::cout << "maxBucketEntryIndex: " << maxNodeIndex << "  maxBucketEntries: " << _processingSoftClauseTree[maxNodeIndex]->inHowManyBuckets << std::endl << std::endl;
    for (auto v : *tmpBucketIndices)
        std::cout << std::setw(4) << v;
    std::cout << std::endl;
    for (uint64_t ind = 0; ind < (*tmpBucketIndices).size(); ind++)
        std::cout << std::setw(4) << "----";
    std::cout << "-------------------------" << std::endl;

    for (auto v : *tmpBucketIndices)
    {
        if (_processingSoftClauseTree[maxNodeIndex]->occursHowOftenInBucket[v] > 0 && v < _processingSoftClauseTree[maxNodeIndex]->highestBucket + 1)
            std::cout << std::setw(4) << _processingSoftClauseTree[maxNodeIndex]->size();
        else
            std::cout << std::setw(4) << "";
    }
    std::cout << "   |" << std::setw(7) << "0"
              << " |" << std::setw(7) << "1"
              << " |" << std::setw(7) << "2"
              << " |" << std::setw(7) << "3"
              << " |" << std::setw(7) << "4"
              << " |" << std::setw(7) << "5" << " |  <-- heuristics" << std::endl;

    for (uint64_t ind = 0; ind < (*tmpBucketIndices).size(); ind++)
        std::cout << std::setw(4) << "----";
    std::cout << "-------------------------" << std::endl;

    int32_t maxCosts(0);
    int32_t estimatedMaxCostIndex(0);

    for (int32_t j=0 ; j < (int)_processingSoftClauseTree.size(); j++)
    {
        if (j == maxNodeIndex)
            continue;

        std::vector< SoftClauseNodes* > NodesToMerge = {_processingSoftClauseTree[maxNodeIndex], _processingSoftClauseTree[j]};
        int32_t usedCosts = GetBenefitOfMergingNodes( NodesToMerge, tmpBucketIndices, true, _setting->groupHeuristic, _setting->percentOff);

        if (usedCosts == 0)
            continue;

        if (maxCosts < usedCosts)
        {
            maxCosts = usedCosts;
            estimatedMaxCostIndex = j;
        }
    }

     for (uint64_t ind = 0; ind < (*tmpBucketIndices).size(); ind++)
        std::cout << std::setw(4) << "----";
    std::cout << "-------------------------" << std::endl;

    for (auto v : *tmpBucketIndices)
        std::cout << std::setw(4) << _totalBucketEntries[v];
    std::cout << std::endl;
    std::cout << std::endl << "MaxCosts: " << maxCosts << "  Index: " << estimatedMaxCostIndex << std::endl;
}

void Cascade::DumpModelOfTares(uint16_t verbosity)
{
    if (_setting->verbosity < verbosity)
        return;

    // if there wasn't a SAT call making the result bigger after constructing the whole cascade, then the tares are 0.
    if (_antom->_lastModel[_structure[0]->_tares[0]] == 0)
        return;
    // output the variables of the Tare T
    std::cout << "Model of Tares: (n...0): ";
    for (int bucketInd = (_structure.size() - 1); bucketInd >= 0; --bucketInd)
    {
        std::cout << "(";
        for (int tareInd = static_cast<int>(_structure[bucketInd]->_tares.size() - 1); tareInd >= 0; --tareInd)
        {
            std::cout << helper::Lit(_antom->_lastModel[_structure[bucketInd]->_tares[tareInd]]);
            if (tareInd != 0)
            {	std::cout << ", "; }
        }
        if (bucketInd != 0)
        { std::cout << "), "; }
        else
        { std::cout << ")"; }
    }
    std::cout << std::endl;
}

void Cascade::DumpBucketSolveInformation(uint32_t actualPos, bool _isLastBucket, uint16_t verbosity)
{
    if (_setting->verbosity < verbosity)
        return;

    if (_isLastBucket)
    {
        _estimatedWeightBoundaries[0] = _highestBucketMultiplicator * actualPos;
        _estimatedWeightBoundaries[1] = _highestBucketMultiplicator * _structure.back()->size();
        std::cout << std::setw(52) << "weight boundaries: ( " << _estimatedWeightBoundaries[0] << " / " << _estimatedWeightBoundaries[1] << " )" << std::endl;
    }
    std::cout << std::setw(52) << "satisfied weights: ( " << _satWeight << " / " << _sumOfSoftWeights << " )" << std::endl;
    std::cout << std::setw(50) << "actualPosition: " << actualPos << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
}

int32_t Cascade::CalculateNodeIndexToMergeWith(uint16_t maxNodeIndex, std::vector<uint16_t> *tmpBucketIndices)
{
    int32_t maxCosts(0);
    // highest possible number to indicate that there is no clause node to merge with!
    int32_t estimatedMaxCostIndex(-1);

    int32_t sumOfSizesOfPercentageFails(0);
    for (uint32_t j = 0 ; j < _processingSoftClauseTree.size(); j++)
    {
        if (j == maxNodeIndex)
            continue;

        std::vector< SoftClauseNodes* > NodesToMerge = {_processingSoftClauseTree[maxNodeIndex], _processingSoftClauseTree[j]};
        int32_t usedCosts = GetBenefitOfMergingNodes(NodesToMerge, tmpBucketIndices, false, _setting->groupHeuristic, _setting->percentOff);
        if (usedCosts == 0)
        {
            continue;
        }
        else if (usedCosts < 0)
        {
            sumOfSizesOfPercentageFails -= usedCosts;
        }

        if (maxCosts < usedCosts)
        {
            maxCosts = usedCosts;
            estimatedMaxCostIndex = j;
        }
    }
    // If no merge is possible
    // check if there is the possibility that the node can be merged later on
    if (estimatedMaxCostIndex == -1 && sumOfSizesOfPercentageFails > 0 && _setting->percentOffReinsert)
    {
        if ((double)_processingSoftClauseTree[maxNodeIndex]->size() * (double)((100 - _setting->percentOff)/100) <= (double)sumOfSizesOfPercentageFails)
        {
            return -2;
        }
    }

    //std::cout << "sumOfSizesOfPercentageFails: " << sumOfSizesOfPercentageFails << std::endl;
    return estimatedMaxCostIndex;
}

int32_t Cascade::GetBenefitOfMergingNodes(std::vector< SoftClauseNodes* > NodesToMerge, std::vector<uint16_t> *tmpBucketIndices, bool dump, uint32_t heuristic, int32_t minPercentOff)
{
    int32_t maxSize = (NodesToMerge[0]->size() > NodesToMerge[1]->size()) ? NodesToMerge[0]->size() : NodesToMerge[1]->size();
    int32_t minSize = (NodesToMerge[0]->size() > NodesToMerge[1]->size()) ? NodesToMerge[1]->size() : NodesToMerge[0]->size();
    int calcPercentOff = 100 - (int)(((double)minSize / (double)maxSize) * 100);
    int32_t binaryClauses = NodesToMerge[1]->size() + NodesToMerge[0]->size();
    int32_t ternaryClauses = NodesToMerge[1]->size() * NodesToMerge[0]->size();
    int32_t clauseCosts = 2 * ternaryClauses + binaryClauses;
    int32_t usedCosts(0);
    int32_t howManyBucketsOverlapping(0);
    int32_t sumOfBucketSizes(1);
    int32_t howOftenBucketsOverlapping(0);
    int32_t bucketSizesFactor(1);
    for (auto v : *tmpBucketIndices)
    {
        if (!(NodesToMerge[1]->occursHowOftenInBucket[v] > 0 && v < NodesToMerge[1]->highestBucket + 1))
            continue;

        bucketSizesFactor+= static_cast<uint32_t>(pow(_totalBucketEntries[v], 1.7) * 0.1);

        howManyBucketsOverlapping += NodesToMerge[1]->size();
        sumOfBucketSizes += _totalBucketEntries[v];
        howOftenBucketsOverlapping++;
    }
    if (howOftenBucketsOverlapping <= 1)
        return 0;

    switch (heuristic)
    {
    // Standard combination of other heuristics
    case 0:
        usedCosts = bucketSizesFactor * static_cast<int32_t>(pow((howOftenBucketsOverlapping-1), 1.7) * pow((clauseCosts), 0.5));
        break;
    // number of reduced (ternary clauses * 2 + binary clauses)
    case 1:
        usedCosts = ( howOftenBucketsOverlapping - 1 ) * clauseCosts;
        break;
    // sum of bucket sizes the SC occurs
    case 2:
        usedCosts = sumOfBucketSizes;
        break;
    // closest size in percentage
    case 3:
        usedCosts = 100 - calcPercentOff;
        break;
    // the one with least occurences in other buckets
    case 4:
        usedCosts = (int)_totalBucketEntries.size() - ((int)NodesToMerge[1]->inHowManyBuckets - howOftenBucketsOverlapping);
        break;
    // how many merges are possible after this merge
    case 5:
        usedCosts = CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge, 0);
        break;
    // greatest depth of submerges (till depth 3, then adds possible merges of depth 4)
    // can be made more effective - but for the ones with big numbers it is way to much work!
    case 6:
        usedCosts = CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge, 1);
        break;
    //
    case 7:
        // usedCosts = CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge[1]);
        break;
    }

    if (minPercentOff != 100)
    {
        // at least one difference - should be always possible, or the calculated percentage.
        if (maxSize - minSize != 1 && calcPercentOff > minPercentOff)
        {
            // if NodesToMerge[0]->size() is bigger, then give back the negative value
            // to calculate if it is theoretically possible to reach that value again
            // if no merge at all is possible.
            if (NodesToMerge[0]->size() > NodesToMerge[1]->size())
            {
                usedCosts = -NodesToMerge[1]->size();
            }
            else
            {
                usedCosts = 0;
            }
        }
    }

    if (!dump)
        return usedCosts;

    int32_t heuristic0 = bucketSizesFactor * static_cast<int32_t>(pow((howOftenBucketsOverlapping-1), 1.7) * pow((clauseCosts), 0.5));
    int32_t heuristic1 = ( howOftenBucketsOverlapping - 1 ) * clauseCosts;
    int32_t heuristic2 = sumOfBucketSizes;
    int32_t heuristic3 = calcPercentOff;
    int32_t heuristic4 = (int)_totalBucketEntries.size() - ((int)NodesToMerge[1]->inHowManyBuckets - howOftenBucketsOverlapping);;
    int32_t heuristic5 =  CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge, 0);
    int32_t heuristic6 =  CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge, 1);


    for (auto v : *tmpBucketIndices)
    {
        if (NodesToMerge[1]->occursHowOftenInBucket[v] > 0 && v < NodesToMerge[1]->highestBucket + 1)
            std::cout << std::setw(4) << NodesToMerge[1]->size();
        else
            std::cout << std::setw(4) << "";
    }
    std::cout << "  ";
    std::cout << " |" << std::setw(7) << heuristic0;
    std::cout << " |" << std::setw(7) << heuristic1;
    std::cout << " |" << std::setw(7) << heuristic2;
    std::cout << " |" << std::setw(7) << heuristic3;
    std::cout << " |" << std::setw(7) << heuristic4;
    std::cout << " |" << std::setw(7) << heuristic5;
    std::cout << " |" << std::setw(7) << heuristic6;
    if (minPercentOff != 100)
    {
        std::cout << " |" << std::setw(7) << calcPercentOff;

        if (usedCosts < 0)
            std::cout << "  <   " << minPercentOff << " <- too big difference";
        else
            std::cout << "  >=  " << minPercentOff;

        std::cout << " |  usedCosts: " << usedCosts << std::endl;
    }
    else
    {
        std::cout << std::endl;
    }


    return usedCosts;
}

int32_t Cascade::CalculateNumberOfPossibleSubmerges(std::vector<uint16_t> *tmpBucketIndices, std::vector< SoftClauseNodes* > NodesToMerge, int32_t depth)
{
    int32_t incomingDepth(depth);
    int32_t possibleSubmerges(0);
    int32_t highestDepth(0);
    std::vector<uint16_t> lastNodeIndices;
    //std::cout << "SNI: ";

    for (auto v : *tmpBucketIndices)
    {
        if (!(NodesToMerge.back()->occursHowOftenInBucket[v] > 0 && v < NodesToMerge.back()->highestBucket + 1))
            continue;
        lastNodeIndices.push_back(v);
         //std::cout << v << ", ";
    }

    for (uint32_t j = 0 ; j < _processingSoftClauseTree.size(); j++)
    {
        // checking if j't element is one of the yet merged ones.
        if (std::find(NodesToMerge.begin(), NodesToMerge.end(), _processingSoftClauseTree[j]) != NodesToMerge.end())
        {
            //std::cout << "sizeofNodesToMerge: " << NodesToMerge.size() << "  Index: " << j << std::endl;
            continue;
        }
        //std::cout << "j: " << j << std::endl;
        //if (_processingSoftClauseTree[j] == NodesToMerge.back())
        //    continue;
        int32_t occurencesInNode(0);
        for (auto v : lastNodeIndices)
        {
            if (_processingSoftClauseTree[j]->highestBucket < v)
                continue;
            //std::cout << "j: " << j << " v: " << v << std::endl;
            if (_processingSoftClauseTree[j]->occursHowOftenInBucket[v] > 0)
                occurencesInNode++;
        }
        if ( depth > 0 && occurencesInNode > 1)
        {
            NodesToMerge.push_back(_processingSoftClauseTree[j]);
            // here we are at depth 3! - go one depth deeper and Calculate then the number of possible submerges!
            if (incomingDepth == 2)
                return (incomingDepth + 1 + CalculateNumberOfPossibleSubmerges(&lastNodeIndices , NodesToMerge, 0));

            int32_t newDepth = CalculateNumberOfPossibleSubmerges(&lastNodeIndices , NodesToMerge, incomingDepth + 1);
            //std::cout << "newDepth: " << newDepth << std::endl;
            if (newDepth > highestDepth)
                highestDepth = newDepth;
        }
        else if ( occurencesInNode > 1)
        {
            //std::cout << "Index: " << j << "  oIN: " << occurencesInNode << std::endl;
            possibleSubmerges++;
        }
    }
    if (highestDepth == 0)
        highestDepth = incomingDepth;
    if (depth > 0)
        return highestDepth;
    else
        return possibleSubmerges;
}

void Cascade::MergeNodes(std::vector<uint16_t> *tmpBucketIndices, int32_t nodeIndexToMergeWith, int32_t maxNodeIndex)
{
    // Exception - there is no node to merge with!!
    if (nodeIndexToMergeWith == -1)
    {
        _softClauseTree.push_back(_processingSoftClauseTree[maxNodeIndex]);
        _processingSoftClauseTree.erase(_processingSoftClauseTree.begin() + maxNodeIndex);
        return;
    }
    // Exception - no node to merge with, but with other merges it can become possible again.
    else if (nodeIndexToMergeWith == -2)
    {
        _processingPercentOffTree.insert(std::make_pair(_processingSoftClauseTree[maxNodeIndex]->size(), _processingSoftClauseTree[maxNodeIndex]));
        _processingSoftClauseTree.erase(_processingSoftClauseTree.begin() + maxNodeIndex);
        return;
    }

    // calc new weight for the two subweights to merge
    uint64_t newWeightforBoth(0);
    for (auto v : *tmpBucketIndices)
    {
        if (_processingSoftClauseTree[nodeIndexToMergeWith]->occursHowOftenInBucket[v] && v < _processingSoftClauseTree[nodeIndexToMergeWith]->highestBucket + 1)
        {
            newWeightforBoth += pow(_base,v);
        }
    }

    Bucket* tmpBucket = new Bucket(_antom, this);
    tmpBucket->AddSoftClauseNode(_processingSoftClauseTree[maxNodeIndex]);
    tmpBucket->AddSoftClauseNode(_processingSoftClauseTree[nodeIndexToMergeWith]);
    SoftClauseNodes* sCNode = new SoftClauseNodes(tmpBucket, newWeightforBoth, _base);
    // at least in two buckets - otherwise no merge
    _processingSoftClauseTree.push_back(sCNode);


    // erase the old SoftClauseNodes.
    uint16_t nodeIndex = maxNodeIndex > nodeIndexToMergeWith ? maxNodeIndex : nodeIndexToMergeWith;

    for (uint16_t i = 0; i < 2; i++)
    {
        if (i == 1)
            nodeIndex = maxNodeIndex > nodeIndexToMergeWith ? nodeIndexToMergeWith : maxNodeIndex;

        if (_processingSoftClauseTree[nodeIndex]->weight - newWeightforBoth == 0)
        {
            _processingSoftClauseTree.erase(_processingSoftClauseTree.begin() + nodeIndex);
        }
        else
        {
            _processingSoftClauseTree[nodeIndex]->setWeight(_processingSoftClauseTree[nodeIndex]->weight - newWeightforBoth);
            // only in one bucket - then move it to the from processingSCT to the SCT
            if (_processingSoftClauseTree[nodeIndex]->inHowManyBuckets == 1)
            {
                _softClauseTree.push_back(_processingSoftClauseTree[nodeIndex]);
                _processingSoftClauseTree.erase(_processingSoftClauseTree.begin() + nodeIndex);
            }
        }
    }

    // if the new weight is bigger than a weight from the _laterToMergeMultimap
    // move the corresponding SCN to the _processing SCT.
    for (std::multimap<uint32_t, SoftClauseNodes*>::iterator it = _processingPercentOffTree.begin(); it != _processingPercentOffTree.end(); ++it)
    {
        if (!AtLeastTwoBucketsInCommon(it->second, _processingSoftClauseTree.back()))
            continue;

        int32_t maxSize = (it->first > _processingSoftClauseTree.back()->size()) ? it->first : _processingSoftClauseTree.back()->size();
        int32_t minSize = (it->first > _processingSoftClauseTree.back()->size()) ? _processingSoftClauseTree.back()->size() : it->first;
        uint32_t calcPercentOff = 100 - (uint32_t)(((double)minSize / (double)maxSize) * 100);

        // Is there a chance of merging again?
        if ( calcPercentOff <= _setting->percentOff )
        {
            if (_setting->verbosity > 3)
            {
                std::cout << "_processingPercentOffTree.size(); " << _processingPercentOffTree.size() << std::endl;
                std::cout << "percentOff: " << calcPercentOff << "   _setting->_percentOff: " << _setting->percentOff << "   minSize: " << minSize << "  maxSize: " << maxSize << std::endl;
            }

            _processingSoftClauseTree.push_back(it->second);
            _processingPercentOffTree.erase(it);
            _howOftenReinsertedFromProcessingPercentOffTree++;
        }
        // because of sorted Multimap - there is no chance that calcPercentOff gets biggerEqual than _antom->groupPercentOff again.
        else if (it->first > _processingSoftClauseTree.back()->size()) {
            break;
        }
    }

    if (_setting->verbosity < 4)
        return;

    std::cout << std::setw(30) << "new Weight: " << newWeightforBoth << std::endl;
}

void Cascade::FillBuckets()
{
    // create the bucket structure.
    for (uint16_t ind = 0; ind <= _numberOfBuckets; ind++)
    {
        Bucket* tmpBucket = new Bucket(_antom, this, ind);
        _structure.push_back(tmpBucket);
    }

    // fill each Bucket according to the _softClauseTree structure
    for (auto sCNode : _softClauseTree)
    {
        for (uint16_t ind = 0; ind <= sCNode->highestBucket; ind++)
        {
//            if (sCNode->occursHowOftenInBucket[ind] == 0)
//                continue;
            //
            for (uint32_t howOften = 0; howOften < sCNode->occursHowOftenInBucket[ind]; howOften++)
            {
                _structure[ind]->AddSoftClauseNode(sCNode);
            }
        }
    }

    _highestBucketMultiplicator = static_cast<uint64_t>(pow( _base, _numberOfBuckets));

    if (_setting->verbosity < 4)
        return;

    std::cout << "Buckets are filled with SoftClauseNodes!" << std::endl;
    DumpBucketStructure(false, 5);
}


void Cascade::AddTaresToBuckets()
{
    assert( _structure.size() > 0 );

//    uint32_t addNoTareToLastBucket = ( _antom->_cascadeDivider > 0 || _onlyByTares ) ? 0 : 1;
    uint32_t addNoTareToLastBucket = 1;

    if (_setting->verbosity > 3)
    {
        std::cout << "addNoTareToLastBucket: " << addNoTareToLastBucket  << std::endl;
        std::cout << "Tare(s): ";
    }
    // add at the last position _base many variables T to each trigger except the last one.
    for (uint32_t i = 0; i < _structure.size() - addNoTareToLastBucket; i++)
    {
        //assert( _antom->_sorterTree[i].size() <= 1 );
        if( !_structure.empty() )
        {
            for (uint32_t j = 0; j < _setting->base - 1; j++)
                AddTare(i);
        }
    }

    if (_setting->verbosity < 2)
        return;

    std::cout << std::endl << "Tares are added to Structure!" << std::endl;
}

void Cascade::EncodeTopBuckets()
{
    TimeMeasurement timeEncoding(&_antom->_timeVariables->encoding, true);
    for (auto bucket : _structure)
    {
        bucket->CalculateNumberOfClauses(true, true, false);

        if (_setting->partitionStrategy == GROUPBYWEIGHTADDATLAST)
        {
            bucket->EncodeTopAddAtLast();
        }
        else
        {
            bucket->EncodeTop();
        }

        bucket->CalculateNumberOfClauses(true, false, true);
    }

    if (_setting->verbosity < 1)
        return;
    if (_setting->verbosity > 3)
        std::cout << std::endl << "Top Buckets are encoded!" << std::endl;

}

void Cascade::EncodeBottomBuckets()
{
    TimeMeasurement timeEncoding(&_antom->_timeVariables->encoding, true);
    for (uint16_t ind = 1; ind < _structure.size(); ind++)
    {
        _structure[ind]->CalculateNumberOfClauses(false, true, false);

        _structure[ind]->MergeSorterWith(_structure[ind-1]->GetEveryNthOutput(_base));

        _structure[ind]->CalculateNumberOfClauses(false, false, true);
    }

    _structure.back()->_isLastBucket = true;

    if (_setting->verbosity < 1)
        return;

    DumpNumberOfBucketsAndClauses();

    if (_setting->verbosity > 3)
        std::cout << "Bottom Buckets are encoded!" << std::endl << std::endl;
}

void Cascade::UnionBucketsIntoLast()
{
    // Push all Buckets to one final Bucket == _structure.back()
    for (uint16_t ind = 1; ind < _structure.size(); ind++)
    {
         _structure[ind - 1]->_nthOutputTaken = _base;
         //std::cout << "_base " << _base << std::endl;
         _structure[ind]->_subBuckets.push_back(_structure[ind -1]);
    }

    if (_setting->verbosity > 3)
        std::cout << "All Buckets in the last bucket!" << std::endl;
}

void Cascade::CreateTotalizerEncodeTree()
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    TimeMeasurement timeEncodeTree(&_antom->_timeVariables->createTree, true);

    _structure.back()->_isLastBucket = true;
    _structure.back()->CreateTotalizerEncodeTree();

    if (_setting->createGraphFile != "")
        _structure.back()->_sorter->_outputTree->DumpOutputTree(_setting->createGraphFile + std::to_string(_structure.back()->size()) + ".tgf", false);


    std::cout << "c #max sorter depth......: " << _structure.back()->_sorter->_outputTree->_depth + 1 << std::endl;
//    std::cout << "SIZE OF TREE: " << _structure.back()->_sorter->_outputTree->_size << std::endl;
    if (_setting->verbosity < 1)
        return;

    if (_setting->verbosity > 3)
        std::cout << std::endl << "Totalizer Tree encoded!" << std::endl;
}

void Cascade::DumpNumberOfBucketsAndClauses()
{
    std::cout << std::endl;
    // Dump BucketID
    DumpNumberOfBucketEntriesOrClauses(false, false, false, false, false, false);
    std::cout << std::endl;

    // Dump TopBucketEntries
    DumpNumberOfBucketEntriesOrClauses(true, false, false, false, false, false);
    // Dump BottomBucketEntries
    DumpNumberOfBucketEntriesOrClauses(false, true, false, false, false, false);
    std::cout << std::endl;

    // Dump Estimated TopBinaryClauses
    DumpNumberOfBucketEntriesOrClauses(true, false, true, false, true, false);
    // Dump Calculated TopBinaryClauses
    DumpNumberOfBucketEntriesOrClauses(true, false, false, true, true, false);
    std::cout << std::endl;

    // Dump Estimated TopTernaryClauses
    DumpNumberOfBucketEntriesOrClauses(true, false, true, false, false, true);
    // Dump Calculated TopTernaryClauses
    DumpNumberOfBucketEntriesOrClauses(true, false, false, true, false, true);
    std::cout << std::endl;

    // Dump Estimated Bottom Binary Clauses
    DumpNumberOfBucketEntriesOrClauses(false, true, true, false, true, false);
    // Dump Calculated Bottom Binary Clauses
    DumpNumberOfBucketEntriesOrClauses(false, true, false, true, true, false);
    std::cout << std::endl;

    // Dump Estimated Bottom Ternary Clauses
    DumpNumberOfBucketEntriesOrClauses(false, true, true, false, false, true);
    // Dump Calculated Bottom Ternary Clauses
    DumpNumberOfBucketEntriesOrClauses(false, true, false, true, false, true);
    std::cout << std::endl;
}

void Cascade::DumpNumberOfBucketEntriesOrClauses(bool top, bool bottom, bool estimated, bool calculated, bool binary, bool ternary)
{
    if (!top && !bottom && !estimated && !calculated)
        std::cout << std::setw(35) << "Bucket ID: ( ";
    else if (top && !estimated && !calculated)
        std::cout << std::setw(35) << "TopBucket Entries: ( ";
    else if (bottom && !estimated && !calculated)
        std::cout << std::setw(35) << "BottomBucket Entries: ( ";
    else if (estimated && top && binary)
        std::cout << std::setw(35) << "Estimated BinaryTopClauses: ( ";
    else if (estimated && top && ternary)
        std::cout << std::setw(35) << "Estimated TernaryTopClauses: ( ";
    else if (calculated && top && binary)
        std::cout << std::setw(35) << "Calculated BinaryTopClauses: ( ";
    else if (calculated && top && ternary)
        std::cout << std::setw(35) << "Calculated TernaryTopClauses: ( ";
    else if (estimated && bottom && binary)
        std::cout << std::setw(35) << "Estimated BinaryBottomClauses: ( ";
    else if (estimated && bottom && ternary)
        std::cout << std::setw(35) << "Estimated TernaryBottomClauses: ( ";
    else if (calculated && bottom && binary)
        std::cout << std::setw(35) << "Calculated BinaryBottomClauses: ( ";
    else if (calculated && bottom && ternary)
        std::cout << std::setw(35) << "Calculated TernaryBottomClauses: ( ";


    uint32_t actualValue(0);
    uint32_t sum(0);

    //for (auto bucket : stuint16_t ind = 0; ind < _structure.size(); ind++)
    for (uint16_t ind = 0; ind < _structure.size(); ind++)
    {
        if (!top && !bottom && !estimated && !calculated)
            actualValue = ind;
        else if (top && !estimated && !calculated)
            actualValue = _structure[ind]->_topEntries;
        else if (bottom && !estimated && !calculated)
            actualValue = static_cast<uint32_t>(_structure[ind]->_outputs->size());
        else if (estimated && top && binary)
            actualValue = _structure[ind]->_binaryTopClEstimated;
        else if (estimated && top && ternary)
            actualValue = _structure[ind]->_ternaryTopClEstimated;
        else if (calculated && top && binary)
            actualValue = _structure[ind]->_binaryTopCl;
        else if (calculated && top && ternary)
            actualValue = _structure[ind]->_ternaryTopCl;
        else if (estimated && bottom && binary)
            actualValue = _structure[ind]->_binaryBottomClEstimated;
        else if (estimated && bottom && ternary)
            actualValue = _structure[ind]->_ternaryBottomClEstimated;
        else if (calculated && bottom && binary)
            actualValue = _structure[ind]->_binaryBottomCl;
        else if (calculated && bottom && ternary)
            actualValue = _structure[ind]->_ternaryBottomCl;

        std::cout << std::setw(5) << actualValue;
        if (ind < _structure.size() - 1)
            std::cout << ", ";
        sum += actualValue;
    }
    if (!top && !bottom && !estimated && !calculated)
        std::cout << " )" << std::endl;
    else
        std::cout << " )   sum =" << std::setw(6) << sum << std::endl;

}

uint32_t Cascade::SolveAllTares()
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    TimeMeasurement timeSolvingTares(&_antom->_timeVariables->solvingTares, true);
    uint32_t currentresult(ANTOM_UNKNOWN);
    std::vector<uint32_t> collectedAssumptions;

    if (_setting->verbosity > 0)
    {
        std::cout << std::endl << std::setw(90) << "Minimize Tares from now on!" << std::endl;
        if (_setting->verbosity > 1)
        {
            std::cout << "----------------temporary solve structure--------------------------------------------------" << std::endl;
        }
    }

    uint16_t sizeMinus = 2;
    if (_onlyByTares)
        sizeMinus = 1;

    //start with second last bucket!
    for (int32_t ind = static_cast<int32_t>(_structure.size() - sizeMinus); ind >= 0; ind--)
    {
        //assert(_estimatedWeightBoundaries[1] >= _satWeight);
        collectedAssumptions.push_back((_structure[ind]->_tares[0] << 1));
        currentresult = _antom->Solve(collectedAssumptions);

//        if (currentresult == ANTOM_SAT && _multipleCascade != nullptr)
//            _multipleCascade->CalculateWeightBoundaries(_structure[ind]->_multiplicator);
//        else if (currentresult == ANTOM_UNSAT && _multipleCascade != nullptr)
//            _multipleCascade->CalculateWeightBoundaries(-_structure[ind]->_multiplicator);

        if (currentresult == ANTOM_SAT)
        {
            if (_setting->verbosity > 3)
                std::cout << "SAT" << std::endl;
//            _antom->_lastModel = _antom->Model();
#ifndef NDEBUG
            bool rst = _antom->AddUnit(collectedAssumptions.back()); assert(rst);
#else
			_antom->AddUnit(collectedAssumptions.back());
#endif
            assert(collectedAssumptions.size() == 1);
            _antom->CalculateOverallOptimum(_satWeight, true);

        }
        else if (currentresult == ANTOM_UNSAT)
        {
            if (_setting->verbosity > 3)
			  {
                std::cout << "UNSAT" << std::endl;
			  }
#ifndef NDEBUG
            bool rst = _antom->AddUnit( collectedAssumptions.back()^1 ); assert(rst);
#else
			_antom->AddUnit( collectedAssumptions.back()^1 );
#endif
        }
        if (_setting->verbosity > 2)
        {
            std::cout << std::setw(50) << "Current Bucket Multiplicator: " << _structure[ind]->_multiplicator << std::endl;
            if (currentresult == ANTOM_SAT)
                std::cout << std::setw(34) << "All Tares of Bucket " << ind << " could be set! " << std::setw(40) << "ANTOM_SAT" << std::endl;
            else if (currentresult == ANTOM_UNSAT)
                std::cout << std::setw(31) << "At least one Tare of Bucket " << ind << " couldn't be set! " << std::setw(40) << "ANTOM_UNSAT" << std::endl;
            else
                std::cout << std::setw(88) << "TIMEOUT: " << currentresult << std::endl;
            currentresult = ANTOM_SAT;
        }
        DumpModelOfTares(5);
        if (_setting->verbosity > 3)
            std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
        collectedAssumptions.pop_back();
    }

    return currentresult;
}

uint32_t Cascade::SolveTares(bool solveTareInLastBucketToo)
{
//    std::cout << __PRETTY_FUNCTION__ << std::endl;
//    std::cout << std::setw(50) << "Weight boundaries: " << _estimatedWeightBoundaries[0] << " / " << _estimatedWeightBoundaries[1] << std::endl;
    if (_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] == 0)
        return ANTOM_SAT;

    TimeMeasurement timeSolvingTares(&_antom->_timeVariables->solvingTares, true);
    uint32_t currentresult(ANTOM_UNKNOWN);

    if (_setting->weightPlusOne)
        return SolveTareWeightPlusOne();

    if (_setting->verbosity > 0)
    {   if (_setting->verbosity > 6)
            std::cout << __PRETTY_FUNCTION__ << std::endl;
        std::cout << std::endl << std::setw(90) << "Minimize Tares from now on!" << std::endl;
        if (_setting->verbosity > 1)
        {
            std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
        }
    }

    uint16_t sizeMinus = 2;
    if (solveTareInLastBucketToo || _onlyByTares)
        sizeMinus = 1;

    //start with second last bucket!
    for (int16_t ind = static_cast<int16_t>(_structure.size() - sizeMinus); ind >= 0; ind--)
    {
        uint64_t diffEstimatedCurWeight;

        diffEstimatedCurWeight = _estimatedWeightBoundaries[1] - _antom->_satWeight;

        currentresult = _structure[ind]->SolveTares(diffEstimatedCurWeight);

        if (_setting->verbosity > 1)
        {
            std::cout << std::setw(50) << "Weight boundaries: " << _estimatedWeightBoundaries[0] << " / " << _estimatedWeightBoundaries[1] << std::endl;
            if (_setting->verbosity > 2)
            {
                std::cout << std::setw(50) << "Current SAT Weight Antom: " << _antom->_satWeight << std::endl;
                std::cout << std::setw(50) << "Current Bucket Multiplicator: " << _structure[ind]->_multiplicator << std::endl;
                if (currentresult == ANTOM_SAT)
                    std::cout << std::setw(34) << "All Tares of Bucket " << ind << " could be set! " << std::setw(40) << "ANTOM_SAT" << std::endl;
                else if (currentresult == ANTOM_UNSAT)
                    std::cout << std::setw(31) << "At least one Tare of Bucket " << ind << " couldn't be set! " << std::setw(40) << "ANTOM_UNSAT" << std::endl;
                else
                    std::cout << std::setw(88) << "TIMEOUT: " << currentresult << std::endl;
            }
            DumpModelOfTares(5);
            if (ind != 0)
                std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
            else
                std::cout << "===========================================================================================" << std::endl;
        }


        // CASE TIMEOUT
        if ( currentresult == ANTOM_UNKNOWN || _control->ReachedLimits() ) {
            return ANTOM_UNKNOWN;
        }
//        std::cout << _estimatedWeightBoundaries[0] << std::endl;
//        std::cout << _antom->_satWeight << std::endl;
        assert(_estimatedWeightBoundaries[0] <= static_cast<int64_t>(_antom->_satWeight));
        //assert(static_cast<int64_t>(_antom->_satWeight) <= _estimatedWeightBoundaries[1]);

    }
    if (_setting->verbosity < 3)
        return ANTOM_SAT;
    DumpModelOfTares(3);
    if (_setting->verbosity > 0)
        std::cout << "All Tares are solved!" << std::endl << std::endl;
    return ANTOM_SAT;
}

void Cascade::CalculateBucketEntries()
{
    uint32_t sumOfSizes(0);
    uint32_t maxBucketEntries(0);

    for (auto bucket : _structure)
    {
        sumOfSizes += bucket->size();
        if (bucket->size() > maxBucketEntries)
            maxBucketEntries = bucket->size();
    }
    std::cout << "c #buckets...............: " << _structure.size() << std::endl;
    std::cout << "c max Bucket entries.....: " << GetMaxBucketSize() << std::endl;
    std::cout << "c average Bucket entries.: " << sumOfSizes / _structure.size() << std::endl;
}


uint64_t Cascade::CountSatisfiedSoftClauses( std::vector< SoftClause* > softclauses, const std::vector<uint32_t>& model, bool addWeight)
{
    uint64_t result(0);
    // Proceed all soft clauses
    for( uint32_t i = 0; i != softclauses.size(); ++i )
    {
        uint32_t relaxlit = softclauses[i]->relaxationLit;

        // Just proceed satisfied triggers
        if( model[relaxlit>>1] == relaxlit )
        {
            std::vector< uint32_t> clause( softclauses[i]->clause );
            uint32_t pos = 0;
            for(; pos != clause.size(); ++pos )
            {
                // clause satisfied without trigger?
                if( model[clause[pos]>>1] == clause[pos] )
                {
                    if (addWeight)
                        result += softclauses[i]->weight;
                    else
                        result += 1;
                    break;
                }
            }
        }
        else if ( model[relaxlit>>1] != 0 )
        {
            assert( model[relaxlit>>1] == (relaxlit^1));
            if (addWeight)
                result += softclauses[i]->weight;
            else
                result += 1;
        }
    }
    return result;
}

uint64_t Cascade::CountSatisfiedSoftClauses(Bucket* bucket, const std::vector<uint32_t>& model)
{
    bool addWeight;
    std::vector< SoftClause* > softclauses;

    if( bucket == NULL || bucket->_isLastBucket )
    {
        addWeight = true;
        softclauses = _softClauses;
    }
    else
    {
        addWeight = false;
        softclauses = *bucket->_softClauses;
//        std::cout << "How many Softclauses: " << softclauses.size() << std::endl;
    }

    uint64_t result = CountSatisfiedSoftClauses(softclauses, model, addWeight);

//    std::cout << "Fulfilled Softclauses: " << result << std::endl;
    if (addWeight)
        _satWeight = result > _satWeight ? result : _satWeight;

    return result;
}

void Cascade::CalculateTotalBucketEntries(std::vector<SoftClauseNodes*>* SCTree, bool add)
{
    _totalBucketEntriesperWeight = 0;
    _totalBucketOccurrences = 0;
    if (!add)
    {
        for (uint32_t i = 0; i <= _numberOfBuckets; ++i)
        {
            if (i >= _totalBucketEntries.size())
                _totalBucketEntries.push_back(0);
            else
                _totalBucketEntries[i] = 0;
        }
    }
    for (uint32_t i = 0; i < SCTree->size(); ++i)
    {
        _totalBucketEntriesperWeight += (*SCTree)[i]->inHowManyBuckets;
        _totalBucketOccurrences += (*SCTree)[i]->GetOccurrences() * (*SCTree)[i]->size();
        for (uint32_t j = 0; j <= (*SCTree)[i]->highestBucket; ++j)
        {
            if ((*SCTree)[i]->hasSubBuckets && (*SCTree)[i]->occursHowOftenInBucket[j])
                _totalBucketEntries[j] += (*SCTree)[i]->occursHowOftenInBucket[j] * (*SCTree)[i]->size();
            else
                _totalBucketEntries[j] += (*SCTree)[i]->occursHowOftenInBucket[j];
        }
    }

}

bool Cascade::AtLeastTwoBucketsInCommon(SoftClauseNodes *SCN1, SoftClauseNodes *SCN2)
{
    uint32_t minSize = (SCN1->occursHowOftenInBucket.size() < SCN2->occursHowOftenInBucket.size()) ? SCN1->occursHowOftenInBucket.size() : SCN2->occursHowOftenInBucket.size();
    uint32_t overlappings(0);
    for ( uint32_t ind = 0; ind < minSize; ind++)
    {
        if (SCN1->occursHowOftenInBucket[ind] > 0 && SCN2->occursHowOftenInBucket[ind] > 0)
        {
            overlappings++;
            if (overlappings > 1)
                return true;
        }
    }
    return false;
}

void Cascade::DumpSCNodeStructure(std::vector<SoftClauseNodes*>* dumpingSCTree, uint16_t verbosity)
{
    if (_setting->verbosity < verbosity)
        return;
    std::cout << std::endl;
    std::cout << std::setw(4) << "" << std::setw(8) << "#" << std::setw(15) << "" << std::setw(3) << "|" << std::endl;
    std::cout << std::setw(4) << "" << std::setw(8) << "O" << std::setw(15) << "" << std::setw(3) << "|" << std::endl;
    std::cout << std::setw(4) << "" << std::setw(8) << "c" << std::setw(15) << "" << std::setw(3) << "|" << std::endl;
    std::cout << std::setw(4) << "" << std::setw(8) << "c" << std::setw(15) << "" << std::setw(3) << "|" << std::endl;
    std::cout << std::setw(4) << "#" << std::setw(8) << "u" << std::setw(15) << "" << std::setw(3) << "|" << std::endl;
    std::cout << std::setw(4) << "B" << std::setw(8) << "r" << std::setw(15) << "" << std::setw(3) << "|" << std::endl;
    std::cout << std::setw(4) << "u" << std::setw(8) << "r" << std::setw(15) << "" << std::setw(3) << "|" << std::endl;
    std::cout << std::setw(4) << "c" << std::setw(8) << "e" << std::setw(15) << "" << std::setw(3) << "|" << std::endl;
    std::cout << std::setw(4) << "k" << std::setw(8) << "n" << std::setw(15) << "" << std::setw(3) << "|" << std::endl;
    std::cout << std::setw(4) << "e" << std::setw(8) << "c" << std::setw(15) << "" << std::setw(3) << "|" << "   " << _base << " up to the power of:" << std::endl;
    std::cout << std::setw(4) << "t" << std::setw(8) << "e" << std::setw(15) << "" << std::setw(3) << "|" << std::endl;
    std::cout << std::setw(4) << "s" << std::setw(8) << "s" << std::setw(15) << "Weight" << std::setw(3) << "|";
    for (uint32_t i = 0; i <= _numberOfBuckets; ++i)
    {
        std::cout << std::setw(4) << i;
    }
    std::cout << std::endl;
    std::cout << "--#occurrences---------------|";

    for (uint32_t i = 0; i <= _numberOfBuckets; ++i)
    {
        std::cout << "----";
    }
    std::cout << "----" << std::endl;

    for (uint32_t i = 0; i != (*dumpingSCTree).size(); ++i)
    {
        (*dumpingSCTree)[i]->dumpStructure(true);
    }
    std::cout << "-----------------------------|";
    for (uint32_t i = 0; i <= _numberOfBuckets; ++i)
    {
        std::cout << "----";
    }
    std::cout << "----" << std::endl;

    std::cout << std::setw(4) << _totalBucketEntriesperWeight << std::setw(8) << _totalBucketOccurrences << std::setw(15) << _sumOfSoftWeights << std::setw(3) << "|";
    for (uint32_t i = 0; i < _totalBucketEntries.size(); ++i)
    {
        std::cout << std::setw(4) << _totalBucketEntries[i];
    }
    std::cout << std::endl;
}

void Cascade::DumpBucketStructure(bool onlyLastBucket, uint16_t verbosity)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    if (_setting->verbosity < verbosity)
        return;
    uint16_t depth(0);
    uint16_t maxDepth(0);
    std::cout << "(SoftClauseEntries, SorterEntries, Multiplicator)" << std::endl;
    if (onlyLastBucket)
    {
        std::cout << "Last Bucket: " << std::endl;
        depth = _structure.back()->DumpAndGetMaxDepth(0);
        std::cout << "LastBucket has a depth of: " << depth << std::endl;
        std::cout << std::endl;
        return;
    }
    for (uint32_t i = 0; i < _structure.size(); i++)
    {
        std::cout << "Bucket " << i << ":" << std::endl;
        depth = _structure[i]->DumpAndGetMaxDepth(0);
        if (depth > maxDepth)
            maxDepth = depth;
        std::cout << "Bucket " << i << " has a depth of: " << depth << std::endl;
        std::cout << std::endl;
    }
    std::cout << "Max depth of all buckets is: " << maxDepth << std::endl << std::endl;
}

uint64_t Cascade::GetHighestMultiplicator()
{
    assert(_highestBucketMultiplicator != 0);
    return _highestBucketMultiplicator;
}

uint32_t Cascade::GetMaxBucketSize()
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    uint32_t maxSize = 0;
    for (auto bucket : _structure)
        maxSize = bucket->size(true) * bucket->_nthOutputTaken > maxSize ? bucket->size(true) * bucket->_nthOutputTaken : maxSize;
//    for (uint32_t ind=0; ind < _structure.size(); ind++)
//    {
//        std::cout << "Bucket " << ind << " has a size of " << _structure[ind]->size(true) * _structure[ind]->_nthOutputTaken << " and every " << _structure[ind]->_nthOutputTaken << " output is taken." << std::endl;
//        maxSize = ((_structure[ind]->size(true) * _structure[ind]->_nthOutputTaken) > maxSize) ? (_structure[ind]->size(true) * _structure[ind]->_nthOutputTaken) : maxSize;
//    }
    return maxSize;
}

bool Cascade::AddNewBucketsTillSizeBoundary(uint32_t maxSize, bool onlySolveWithTares, bool addTareToLastBucket)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    if (_setting->verbosity > 3)
        std::cout << "Add another bucket because of size boundary: " << maxSize << std::endl;


    while (_structure.back()->size(true) > maxSize)
    {
        if (!AddNewBucketsTillMultiplicatorMatches(_highestBucketMultiplicator * 2, onlySolveWithTares, addTareToLastBucket))
            return false;
    }
    return true;
}

bool Cascade::AddNewBucketsTillMultiplicatorMatches(uint64_t maxMultiplicator, bool onlySolveWithTares, bool addTareToLastBucket)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    if (_setting->verbosity > 3)
        std::cout << "_highestBucketMultiplicator: " << _highestBucketMultiplicator << std::endl;
    assert(_structure.back()->_multiplicator <= maxMultiplicator);

    while (_structure.back()->_multiplicator < maxMultiplicator)
    {
        if (_setting->verbosity > 3)
        {
            std::cout << "_structure.back()->size(true): " << _structure.back()->size(true) << std::endl;
            std::cout << "onlySolveWithTares: " << onlySolveWithTares << std::endl;
        }


        // This bucket is too small to connect with another bucket.
        // It contains only tare and one result from the bucket before!
        // Calculate the results only by the tares!
        //TOBI: VERIFY!!! - Should be possible to calculate all results only by tares
        if (_structure.back()->size(true) == 1 && onlySolveWithTares)
        {
//            if (_setting->verbosity > 3)
//                std::cout << "c MORE THAN TWO ENTRIES..: TRUE" << std::endl;
            if (_setting->verbosity > 3)
                std::cout << "Add one tare to last bucket with size one!" << std::endl;
            AddTare(_structure.size() - 1);

            // here is the last position the bucket can be dumped
            // otherwise the subbuckets are encoded and not dumpable anymore :-)
            DumpBucketStructure(true, 4);

            CreateTotalizerEncodeTree();
            CalculateBucketEntries();

            // at least the tare should be 0, then by solving max tares to 1 is asked for!
            uint32_t unitClauseVar = (_structure.back()->_sorter->GetOrEncodeOutput(0) << 1) ^ 1;

            if (_setting->verbosity > 3)
                std::cout << "Unit Clause for first entry in last Bucket Added: " << unitClauseVar << std::endl;
            _antom->AddUnit(unitClauseVar);
            return false;
        } else
        {
            //if (_structure.back()->size(true) >= 1)
            AddAdditionalBucket();
        }

        // at least one element + one tare!!
        assert(_structure.back()->size(true) > 0);


    }

    if(addTareToLastBucket && _structure.back()->_tares.empty())
    {
        if (_setting->verbosity > 3)
            std::cout << std::endl << "There is no tare at the last bucket! Add one!" << std::endl;
        AddTare(_structure.size() - 1);
    }

    return true;
}

void Cascade::AddTare(unsigned long position)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    uint32_t tare(_antom->NewVariable());
    _antom->SetDontTouch(tare);
    _structure[position]->AddTare(tare);
    _tareWeight += _structure[position]->_multiplicator;
    if (_setting->verbosity < 3)
        return;

    std::cout  << "   Tare added: " << tare << std::endl;
}

void Cascade::AddAdditionalBucket()
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    if (_setting->verbosity > 3)
    {
        std::cout << "_highestBucketMultiplicator: " << _highestBucketMultiplicator << std::endl;
        std::cout << "Add another Bucket with multiplicator: " << _highestBucketMultiplicator;
    }
    _structure.back()->_isLastBucket = false;
    if(_structure.back()->_tares.empty())
    {
        if (_setting->verbosity > 3)
            std::cout << std::endl << "There is no tare at the last bucket! Add one!" << std::endl;
        AddTare(_structure.size() - 1);
    }

    _structure.push_back(new Bucket(_antom, this));
    _structure.back()->_isLastBucket = true;
    _highestBucketMultiplicator = static_cast<uint64_t>(pow(_base, _structure.size() - 1));
    _structure.back()->_multiplicator = _highestBucketMultiplicator;
    _structure[_structure.size()-2]->_nthOutputTaken = _base;
    _structure.back()->_subBuckets.push_back(_structure[_structure.size()-2]);
    if (_structure.back()->_subBuckets.back()->_localMaxPos != static_cast<uint32_t>(-1))
        _structure.back()->_localMaxPos = static_cast<uint32_t>(floor(_structure.back()->_subBuckets.back()->_localMaxPos / 2));
}


std::vector<uint32_t> Cascade::CalculateAssumptionsFor(int64_t weight, int32_t startingPos)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    assert(weight <= _estimatedWeightBoundaries[1]);
    assert(weight > _estimatedWeightBoundaries[0]);
    assert(weight > static_cast<int64_t>(_antom->_satWeight));

    std::vector<uint32_t> collectedAssumptions;

    int64_t upperDiff = _estimatedWeightBoundaries[1] - weight;
    int64_t lowerDiff = weight - _estimatedWeightBoundaries[0];


    //start with second last bucket!
    for (int32_t ind = startingPos; ind >= 0; --ind)
    {
        int64_t actualMult = static_cast<int64_t>(_structure[ind]->_multiplicator);

        if (_setting->verbosity > 5)
        {
            std::cout << std::setw(50) << "Actual Diffs: " << lowerDiff << " / " << upperDiff << std::endl;
            std::cout << std::setw(50) << "actualMult: " << actualMult << std::endl;
        }

        assert(upperDiff + lowerDiff == 2 * actualMult - 1 );

        if (lowerDiff >= actualMult)
        {
            assert(upperDiff < actualMult);
            collectedAssumptions.push_back(_structure[ind]->_tares[0] << 1);
            lowerDiff -= actualMult;
        }
        else if (upperDiff >= actualMult)
        {
            assert(lowerDiff < actualMult);
            collectedAssumptions.push_back(_structure[ind]->_tares[0] << 1 ^ 1);
            upperDiff -= actualMult;
        }
        else
        {
            assert(false);
        }
    }


    if (_setting->verbosity < 3)
        return collectedAssumptions;

    std::cout << std::setw(51) << "Assumptions: (";
    for(auto assumption : collectedAssumptions)
    {
         std::cout << assumption << ", ";
    }
    std::cout << ")" << std::endl;
    return collectedAssumptions;
}


int32_t Cascade::SetUnitClauses(int32_t startingPos)
{
    if (_setting->verbosity > 6)
	  {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
	  }

    //start with second last bucket!
    for (int32_t ind = startingPos; ind >= 0; ind--)
    {
        int64_t actualMult = static_cast<int64_t>(_structure[ind]->_multiplicator);

        if (_setting->verbosity > 5)
        {
            std::cout << std::setw(50) << "Weight boundaries: " << _estimatedWeightBoundaries[0] << " / " << _estimatedWeightBoundaries[1] << std::endl;
            std::cout << std::setw(50) << "actualMult: " << actualMult << std::endl;
        }

        assert(_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] == 2 * actualMult - 1);
        // The actual result implies that the currentTare and maybe more tares can be set directly to TRUE
        if (_estimatedWeightBoundaries[1] - static_cast<int64_t>(_antom->_satWeight) < actualMult)
        {
//            std::cout << _structure[ind]->_tares[0]<< std::endl;
#ifndef NDEBUG
            bool rst = _antom->AddUnit(_structure[ind]->_tares[0] << 1); assert(rst);
#else
			_antom->AddUnit(_structure[ind]->_tares[0] << 1);
#endif
            _estimatedWeightBoundaries[0] += actualMult;
            continue;
        }
        // Corner Case!
        // If new result is already larger as maximal possible weight, we do not need to try,
        // -> the corresponding tare can be directly set to FALSE
        else if (_estimatedWeightBoundaries[1] - actualMult >= static_cast<int64_t>(_antom->_sumOfSoftWeights))
        {
#ifndef NDEBUG
            bool rst = _antom->AddUnit((_structure[ind]->_tares[0] << 1) ^ 1); assert(rst);
#else
			_antom->AddUnit((_structure[ind]->_tares[0] << 1) ^ 1);
#endif
            
            _estimatedWeightBoundaries[1] -= actualMult;
            continue;
        }
        // set assumption vector.
        else
        {
            if (_setting->verbosity > 2)
            {
                std::cout << std::setw(50) << "Weight boundaries: " << _estimatedWeightBoundaries[0] << " / " << _estimatedWeightBoundaries[1] << std::endl;
                std::cout << std::setw(50) << "actualMult: " << actualMult << std::endl;
            }
            return ind;
        }
    }

    if (_estimatedWeightBoundaries[1] - _antom->_satWeight > 0)
        return 0;
    else
        // result is found.
        return -1;
}


uint32_t Cascade::SolveTareWeightPlusOne()
{
    TimeMeasurement timeSolvingTares(&_antom->_timeVariables->solvingTares, true);
    uint32_t currentresult(ANTOM_SAT);
    std::vector<uint32_t> collectedAssumptions;

    if (_setting->verbosity > 0)
    {
        if (_setting->verbosity > 6)
            std::cout << __PRETTY_FUNCTION__ << std::endl;
        std::cout << std::endl << std::setw(90) << "Try to set tares for actual weight + 1!" << std::endl;
        if (_setting->verbosity > 1)
            std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
    }

    int32_t startingPos = _structure.size() - 2;

    if (_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] > static_cast<int64_t>(_highestBucketMultiplicator))
        startingPos++;
//    else if (_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] == 1 && _highestBucketMultiplicator <= 2)
//        //border case, only 2 buckets... est weight diffs = 1
//        startingPos++;
    else
        assert(_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] >= static_cast<int64_t>(_highestBucketMultiplicator / 2));

//    std::cout << "startingPos: " << startingPos << std::endl;

    while (currentresult == ANTOM_SAT)
    {
        startingPos = SetUnitClauses(startingPos);

        if (startingPos == -1 || static_cast<int64_t>(_antom->_satWeight) == _estimatedWeightBoundaries[1])
            break;

        collectedAssumptions = CalculateAssumptionsFor(static_cast<int64_t>(_antom->_satWeight) + 1, startingPos);

        if( _setting->maxInprocess && _antom->Preprocess(INPROCESS) != ANTOM_UNKNOWN)
            break;

        currentresult = _antom->Solve(collectedAssumptions);
        if (_setting->verbosity > 1)
        std::cout << "-------------------------------------------------------------------------------------------" << std::endl;

//        std::cout << "CR: " << currentresult << std::endl;
        if (currentresult == ANTOM_SAT)
            _antom->CalculateOverallOptimum(0,true);
    }

    if (_setting->verbosity > 1)
        std::cout << "===========================================================================================" << std::endl;

    if (currentresult == ANTOM_UNKNOWN || _control->ReachedLimits() )
    {
        _antom->_resultUnknown = true;
        std::cout << std::setw(88) << "TIMEOUT: " << currentresult << std::endl;
        return ANTOM_UNKNOWN;
    }

    if (_setting->verbosity > 0)
        std::cout << "All Tares are solved!" << std::endl << std::endl;
    return ANTOM_SAT;
}

}
