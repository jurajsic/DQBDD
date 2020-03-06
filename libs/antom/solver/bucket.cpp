/********************************************************************************************
bucket.cpp -- Copyright (c) 2017, Tobias Paxian

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
#include <iomanip>
#include <cmath>
#include <map>

#include "bucket.h"
#include "cascade.h"
#include "multiplecascade.h"
#include "sorter.h"
#include "control.h"
#include "softclausenodes.h"
#include "totalizerencodetree.h"
#include "timemeasurement.h"
#include "timevariables.h"


namespace antom
{

// Constructor
Bucket::Bucket(Antom* antom, Cascade* cascade, uint16_t position) :
    _howOftenUsed(1),
	_sorter(NULL),
    _antom(antom),
	_setting(antom->_antomSetting),
	_control(antom->_control),
    _cascade(cascade),
	_subBuckets(),
	_softClauses(),
	_outputs(),
	_tares(),
	_bucketAssumptions(),
    _tarePosition(0),
    _bucketPosition(position),
	_base(0),
    _nthOutputTaken(1),
    _topEntries(0),
    _binaryTopCl(0),
    _binaryTopClEstimated(0),
    _binaryBottomCl(0),
    _binaryBottomClEstimated(0),
    _ternaryTopCl(0),
    _ternaryTopClEstimated(0),
    _ternaryBottomCl(0),
    _ternaryBottomClEstimated(0),
	_multiplicator(0),
    _groundWeight(0),
    _localMaxPos(static_cast<uint32_t>(-1)),
    _positionsCuttedAtBottom(0),
    _weightBoundaryAccuracy(0),
    _upperWeightBoundOfTareCascade(0),
    _sumSoftWeightsOfTareCascade(0),
    _solvingState(SINGLECASCADE),
    _encoded(false),
    _encodeTreeGenerated(false),
    _isLastBucket(false)
{
    assert(antom != NULL);
    // As long as we take the same base for the whole cascade.
    _base = _cascade->_base;
    if (position != static_cast<uint16_t>(~0))
        _multiplicator = static_cast<uint64_t>(pow(_base, position));
    else
        _multiplicator = 0;

    _sorter = new Sorter(0, _antom);
    _softClauses = &_sorter->_softClauses;
    _outputs = &_sorter->_outputs;
}

// Destructor
Bucket::~Bucket(void)
{}

uint32_t Bucket::size(bool reduceByFactorNthOutput)
{
    if (_encodeTreeGenerated)
        return _sorter->_outputTree->_size / _sorter->_outputTree->_everyNthOutput;

    uint32_t recursiveSorterSizes = 0;

    if (_sorter->_outputTree != nullptr)
        recursiveSorterSizes += _sorter->_outputTree->_size / _sorter->_outputTree->_everyNthOutput;
    else
        recursiveSorterSizes = static_cast<uint32_t>(_sorter->size());

    for (uint32_t i = 0; i < _subBuckets.size(); i++)
    {
        recursiveSorterSizes += _subBuckets[i]->size(reduceByFactorNthOutput);
    }
    if (reduceByFactorNthOutput && _nthOutputTaken > 1)
    {
        recursiveSorterSizes /= _nthOutputTaken;
    }
    return recursiveSorterSizes;

}

void Bucket::AddSoftClauseNode(SoftClauseNodes *SCNode)
{
    if (SCNode->hasSoftClause)
    {
        assert(!SCNode->hasSubBuckets);
        _sorter->AddSoftClause(SCNode->softClause);
        _sorter->AddOutput(SCNode->softClause->relaxationLit >> 1);
    }
    else if (SCNode->hasSubBuckets)
    {
        assert(!SCNode->hasSoftClause);
        for (auto v : SCNode->subBuckets)
        {
            _subBuckets.push_back(v);
        }
    }
    else
    {
        assert(false);
    }
}

void Bucket::AddTare(uint32_t tare)
{
    _sorter->AddOutput(tare);
    _antom->SetDontTouch(tare);
    _tares.push_back(tare);
}

uint16_t Bucket::DumpAndGetMaxDepth(uint16_t depth)
{
    uint16_t currentDepth(depth);
    uint16_t maxDepth(depth);

    for (uint16_t i = 0; i < depth; i++)
        std::cout << std::setw(10) << "";

    // strange multiplicator result happening over and over again. Just reset it to 0
//    if (_multiplicator > 1000000000000)
//        _multiplicator = 0;

    std::cout << "(" << size(true) << ", " << _sorter->size() << ", " << _multiplicator << ", " << _tares.empty() << ")" << std::endl;

    for (auto subBucket : _subBuckets)
    {
        currentDepth = subBucket->DumpAndGetMaxDepth(depth + 1);
        if (currentDepth > maxDepth)
            maxDepth = currentDepth;
    }
    return maxDepth;
}

void Bucket::EncodeTop(uint32_t howOftenReused)
{
    //std::cout << __FUNCTION__ << std::endl;

    /**
     * @brief sorterAllSizes
     * A list of all sizes with bool to indicate if sorter can be used to merge with.
     * false -> cannot be merged with (lookup in bucket list)
     *  true -> can be merged with (lookup in merge list)
     */
    std::multimap<uint32_t, bool> sorterAllSizes;

    /**
     * @brief sorterBucketSizes
     * This sorter cannot be merged with directly because they are still used in other buckets.
     */
    std::multimap<uint32_t, Sorter*> sorterBucketSizes;

    //
    /**
     * @brief sorterMergeSizes
     * With these sorters can be directly merged with.
     */
    std::multimap<uint32_t, Sorter*> sorterMergeSizes;

    assert((int)_howOftenUsed - (int)howOftenReused >= 0);
    _howOftenUsed -= howOftenReused;

    if (_encoded)
    {
        return;
    }
    if (_sorter->size() != 0)
    {
        //encode main Sorter
        _sorter->EncodeSorter(0);
        if (_subBuckets.size() == 1)
        {
            _subBuckets[0]->EncodeTop(_howOftenUsed);
            _sorter->SimpleMergeWithSorter(*_subBuckets[0]->_sorter);
            _subBuckets.clear();
            _encoded = true;
            return;
        }
        else if (_subBuckets.size() == 0)
        {
            _encoded = true;
            return;
        }

        sorterMergeSizes.insert( std::make_pair( _sorter->size(), _sorter ) );
        sorterAllSizes.insert( std::make_pair( _sorter->size(), true ) );

    } else {
        assert(_subBuckets.size() != 0);
        if (_subBuckets.size() == 1)
        {
            _subBuckets[0]->EncodeTop(_howOftenUsed);
            _sorter->SimpleMergeWithSorter(*_subBuckets[0]->_sorter);
            _subBuckets.pop_back();
            _encoded = true;
            return;
        }
    }

    //encode Sub Sorter
    for (uint64_t ind = _subBuckets.size() - 1; ind < _subBuckets.size(); ind--)
    {
        _subBuckets[ind]->EncodeTop(_howOftenUsed);

        if(_subBuckets[ind]->_howOftenUsed > 0 || _base > 2) {
            sorterBucketSizes.insert( std::make_pair( _subBuckets[ind]->size(), _subBuckets[ind]->_sorter ) );
            sorterAllSizes.insert( std::make_pair( _subBuckets[ind]->size(), false ) );
        } else if (_subBuckets[ind]->_howOftenUsed == 0) {
            sorterMergeSizes.insert( std::make_pair( _subBuckets[ind]->size(), _subBuckets[ind]->_sorter ) );
            sorterAllSizes.insert( std::make_pair( _subBuckets[ind]->size(), true ) );
        }
        _subBuckets.pop_back();
    }
    assert(_subBuckets.empty());
//    for (std::map<uint32_t,Sorter*>::iterator it=sorterMergeListbySizes.begin(); it!=sorterMergeListbySizes.end(); ++it)
//        std::cout << it->first << " => " << &it->second << '\n';

    while (!sorterBucketSizes.empty())
    {
        //std::cout << "sorterAllSizes.size(): " << sorterAllSizes.size() << std::endl;
        //std::cout << "sorterBucketSizes.size(): " << sorterBucketSizes.size() << std::endl;
        //std::cout << "sorterMergeSizes.size(): " << sorterMergeSizes.size() << std::endl;

        assert(sorterAllSizes.size() == sorterBucketSizes.size() + sorterMergeSizes.size());
        Sorter* mergeIntoSorter;
        Sorter* mergeSorter;
        // into the smallest sorter can be merged into!
        if ( sorterAllSizes.begin()->second )
        {
            assert(sorterAllSizes.begin()->first == sorterMergeSizes.begin()->first);
            mergeIntoSorter = sorterMergeSizes.begin()->second;
            sorterMergeSizes.erase(sorterMergeSizes.begin());
            sorterAllSizes.erase(sorterAllSizes.begin());
        }
        // into the second smallest sorter can be merged into!
        else if ( (++sorterAllSizes.begin())->second )
        {
            assert((++sorterAllSizes.begin())->first == sorterMergeSizes.begin()->first);
            mergeIntoSorter = sorterMergeSizes.begin()->second;
            sorterMergeSizes.erase( sorterMergeSizes.begin() );
            sorterAllSizes.erase( ++sorterAllSizes.begin() );
        }
        //  the sorter has to be copied in a new sorter.
        else
        {
            assert(sorterAllSizes.begin()->first == sorterBucketSizes.begin()->first);
            mergeIntoSorter = new Sorter(0, _antom);
            mergeIntoSorter->SimpleMergeWithSorter(*sorterBucketSizes.begin()->second);
            sorterBucketSizes.erase( sorterBucketSizes.begin() );
            sorterAllSizes.erase( sorterAllSizes.begin() );
        }


        if ( sorterAllSizes.begin()->second )
        {
            assert(sorterAllSizes.begin()->first == sorterMergeSizes.begin()->first);
            mergeSorter = sorterMergeSizes.begin()->second;
            sorterMergeSizes.erase(sorterMergeSizes.begin());
        }
        else
        {
            assert(sorterAllSizes.begin()->first == sorterBucketSizes.begin()->first);
            mergeSorter = sorterBucketSizes.begin()->second;
            sorterBucketSizes.erase(sorterBucketSizes.begin());
        }
        sorterAllSizes.erase(sorterAllSizes.begin());

        mergeIntoSorter->SimpleMergeWithSorter(*mergeSorter);
        sorterMergeSizes.insert( std::make_pair( mergeIntoSorter->size(), mergeIntoSorter ) );
        sorterAllSizes.insert( std::make_pair( mergeIntoSorter->size(), true ) );
    }

    while (sorterMergeSizes.size() > 1)
    {
        // merge the first Sorter into the second.
        (++sorterMergeSizes.begin())->second->SimpleMergeWithSorter(*sorterMergeSizes.begin()->second);

        // erase the smallest one.
        sorterMergeSizes.erase(sorterMergeSizes.begin());

        // include the second sorter again with the correct size.
        sorterMergeSizes.insert( std::make_pair( sorterMergeSizes.begin()->second->size(), sorterMergeSizes.begin()->second ) );
        sorterAllSizes.insert( std::make_pair( sorterMergeSizes.begin()->second->size(), true ) );

        // delete the first element which has the wrong size because the first Sorter was merged into it.
        sorterMergeSizes.erase(sorterMergeSizes.begin());
    }

    if (_sorter->size() != sorterMergeSizes.begin()->first)
    {
        _outputs->clear();
        _softClauses->clear();
        _sorter->SimpleMergeWithSorter(*sorterMergeSizes.begin()->second);
    }
    _encoded = true;
}

void Bucket::EncodeTopAddAtLast(void)
{
    if (_encoded)
        return;
    if (_sorter->size() != 0)
        _sorter->EncodeSorter(0);

    for (uint64_t ind = _subBuckets.size() - 1; ind < _subBuckets.size(); ind--)
    {
        _subBuckets[ind]->EncodeTopAddAtLast();
        //TOASK: does Bitonic Sort create additional vars to get two subBuckets with same size?
        // the same number of clauses!
        //_sorter->MergeWithSorter(*_subBuckets[ind]->_sorter);
        //std::vector<uint32_t> tmpoutputs = _subBuckets[ind]->_sorter->GetOutputs();
        //_sorter->MergeTotalizer(tmpoutputs);
        // simple merger doesn't do the bypass stuff!
        _sorter->SimpleMergeWithSorter(*_subBuckets[ind]->_sorter);

        _subBuckets.pop_back();
    }
    _encoded = true;
}

void Bucket::CreateTotalizerEncodeTree()
{
//    std::cout << __PRETTY_FUNCTION__ << std::endl;

    /**
     * @brief sorterBucketSizes
     * This sorter cannot be merged with directly because they are still used in other buckets.
     */
    std::multimap<uint32_t, TotalizerEncodeTree*> sorterSizes;

    if(_encodeTreeGenerated)
        return;

    if (_sorter->size() != 0)
    {
        _sorter->CreateTotalizerEncodeTree();
//        _sorter->_outputs.clear();
        if (_subBuckets.size() > 0)
            sorterSizes.insert( std::make_pair( _sorter->size(), _sorter->_outputTree ) );
        else if (_sorter->_outputTree->_size == 1)
        {
            _sorter->_outputTree->_everyNthOutput = _nthOutputTaken;
            _sorter->_outputTree->_hasBeenBucketBefore = true;
            return;
        }
    } else {
        //std::cout << "Sorter is empty!" << std::endl;
        assert( _subBuckets.size() > 0 );
    }


    for (uint64_t ind = _subBuckets.size() - 1; ind < _subBuckets.size(); ind--)
    {
        _subBuckets[ind]->CreateTotalizerEncodeTree();

        //sorterSizes.insert( std::make_pair( _subBuckets[ind]->size(true), _subBuckets[ind]->_sorter->_outputTree ) );
        sorterSizes.insert( std::make_pair( _subBuckets[ind]->_sorter->_outputTree->_size, _subBuckets[ind]->_sorter->_outputTree ) );

        _subBuckets.pop_back();
    }

//    for(auto const &element : sorterSizes)
//        std::cout << element.first << " => " << element.second->_size << std::endl;
//    std::cout << std::endl;

    while (sorterSizes.size() > 1)
    {
        TotalizerEncodeTree* TotTree = new TotalizerEncodeTree(0);
        TotTree->_child1 = sorterSizes.begin()->second;
        sorterSizes.erase(sorterSizes.begin());
        TotTree->_child2 = sorterSizes.begin()->second;
        sorterSizes.erase(sorterSizes.begin());
        TotTree->ActualizeValues();
        sorterSizes.insert(std::make_pair( TotTree->_size, TotTree ));
    }

    if (sorterSizes.size() == 1)
        _sorter->_outputTree = sorterSizes.begin()->second;
    _sorter->_outputTree->_everyNthOutput = _nthOutputTaken;
    _sorter->_outputTree->_howOftenUsed = _howOftenUsed;
    _sorter->_outputTree->ActualizeValues();
    _sorter->_outputTree->_hasBeenBucketBefore = true;

    if (size() > _setting->maxWidth && _setting->maxWidth != 0 && !_isLastBucket)
    {
        _encodeTreeGenerated = true;
        CutMaxPos(true);
    }
    _encodeTreeGenerated = true;
}

void Bucket::MergeSorterWith(std::vector<uint32_t> MergeVector)
{
    // TOASK: Here only Totalizer - or better use MergeWithSorter?!
    _sorter->MergeTotalizer(MergeVector);
}
//void Bucket::SetSolvingParameters(std::vector<uint32_t>* assumptions, uint64_t weightAccuracy, uint64_t weightToSubstract, uint64_t _upperWeightBoundOfTareCascade, uint64_t _weightOfAllTareCascadeTares, bool SolvingTares)
void Bucket::SetSolvingParameters(std::vector<uint32_t>* assumptions,
                                  uint64_t upperWeightBoundOfTareCascade,
                                  uint64_t weightOfAllTareCascadeTares,
                                  int64_t groundWeight,
                                  MultipleCascadeSolvingState solveState)
{
    if (_setting->verbosity > 6)
        std::cout << std::endl << __PRETTY_FUNCTION__ << std::endl;

    _bucketAssumptions = *assumptions;
//    _weightBoundaryAccuracy = weightAccuracy;
    _groundWeight = groundWeight;
    _upperWeightBoundOfTareCascade = upperWeightBoundOfTareCascade;
    _sumSoftWeightsOfTareCascade = weightOfAllTareCascadeTares;
    _solvingState = solveState;

    if (_setting->verbosity < 2)
        return;

    std::cout << std::setw(50) << "_weightBoundaryAccuracy: " << _weightBoundaryAccuracy << std::endl;
    std::cout << std::setw(50) << "_groundWeight: " << _groundWeight << std::endl;
    std::cout << std::setw(50) << "_upperWeightBoundOfTareCascade: " << _upperWeightBoundOfTareCascade << std::endl;
    std::cout << std::setw(50) << "_sumSoftWeightsOfTareCascade: " << _sumSoftWeightsOfTareCascade << std::endl;
}

uint32_t Bucket::CutMaxPos(bool solve)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    if (solve)
    {
        std::cout << "MaxPos before: " << _localMaxPos << std::endl;
        if (!_encodeTreeGenerated)
            CreateTotalizerEncodeTree();
        _localMaxPos = SolveBucketReturnMaxPosition(true, true);
        std::cout << "MaxPos after: " << _localMaxPos << std::endl;
    }
    if (_setting->createGraphFile != "")
        _sorter->_outputTree->DumpOutputTree(_setting->createGraphFile + "_beforeSettingMaxPos"+ std::to_string(_localMaxPos)+".tgf", true);
    SetMaxPos(_localMaxPos + 1);

    // because of local calc, that other buckets can be added afterwards!
    _encodeTreeGenerated = false;
    return _localMaxPos;
}

void Bucket::SetMaxPos(uint32_t maxPos)
{
//    while(_sorter->size() > maxPos)
//    {
//        _sorter->_outputs.pop_back();
//    }
//    std::cout << "maxPos: " << maxPos << std::endl;
//    std::cout << "SORTERSIZE: " << _sorter->size() << std::endl;


    std::cout << "      BUCKETSIZE: " << size() << std::endl;
    std::cout << "          MaxPos: " << maxPos << std::endl;

    if (maxPos + 1 < size())
        SetAsUnitClause(maxPos + 1, ANTOM_UNSAT, false);

    if (_sorter->_outputTree != nullptr)
        _sorter->_outputTree->CutVectorAbove(maxPos + 1);

//    std::cout << "BUCKETSIZE TRUE: " << size(true) << std::endl;
//    std::cout << "BUCKETSIZE FALSE: " << size(false) << std::endl;
}

uint32_t Bucket::CutMinPos(bool solve)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    if (solve)
    {
        std::cout << "MaxPos before: " << _localMaxPos << std::endl;
        CreateTotalizerEncodeTree();
        _localMaxPos = SolveBucketReturnMaxPosition(true, true);
        std::cout << "MaxPos after: " << _localMaxPos << std::endl;
    }
    int32_t positionsToCutFromMaxPos = static_cast<int32_t>(ceil(_cascade->_upperWeightBoundAllLowerCascades / _multiplicator)) + 1;

    std::cout << "_upperWeightBoundAllLowerCascades: " << _cascade->_upperWeightBoundAllLowerCascades << std::endl;
    std::cout << "ANTOM SoftWeights: " << _antom->_sumOfSoftWeights << std::endl;
    std::cout << "positionsToCutFromMaxPos: " << positionsToCutFromMaxPos << std::endl;
    int32_t minPosition = static_cast<int32_t>(_localMaxPos) - positionsToCutFromMaxPos;
    std::cout << "                  minPos: " << minPosition << std::endl;

    _encodeTreeGenerated = false;

    if (minPosition <= 0)
//    if (true)
        return _localMaxPos;

    uint32_t minPos = static_cast<uint32_t>(minPosition);
    assert(minPos <= size() );

    SetAsUnitClause(minPos, ANTOM_SAT, false);

    // the unit clause minPos should be cutted! it is fullfilled in either case!
    if (_setting->createGraphFile != "")
        _sorter->_outputTree->DumpOutputTree(_setting->createGraphFile + "_beforeSettingMinPos" + std::to_string(minPos) + ".tgf", true);
    if (_sorter->_outputTree != nullptr)
    {
        _sorter->_outputTree->CutVectorBelow(minPos);
        _positionsCuttedAtBottom += minPos;
    }

    _sorter->_outputTree->DumpOutputTree(_setting->createGraphFile + "_afterSettingMinPos" + std::to_string(minPos) + ".tgf", true);
    _localMaxPos = _localMaxPos - minPos > 0 ? _localMaxPos - minPos : static_cast<uint32_t>(-1);

    std::cout << "Positions Cutted at Bottom: " << _positionsCuttedAtBottom << std::endl;
    std::cout << "BUCKETSIZE TRUE: " << size(true) << std::endl;
    std::cout << "BUCKETSIZE FALSE: " << size(false) << std::endl;


    return _localMaxPos;
}

uint32_t Bucket::SolveBucketReturnMaxPosition(bool onlyWithAssumptions, bool localCalc)
{
//    uint32_t actualPos = _localMaxPos == static_cast<uint32_t>(-1) ? 0 : _localMaxPos;
//    uint32_t lastPos = _localMaxPos == static_cast<uint32_t>(-1) ? 0 : _localMaxPos;
    uint32_t actualPos = 0;
    uint32_t lastPos = 0;
    uint32_t currentresult(ANTOM_SAT);
    uint64_t currSatWeight(0);
    TimeMeasurement TimeSolvingLastBucket(&_antom->_timeVariables->solvingLastBucket, true);
    uint32_t solverCallsStartSolveBucket = _antom->_satSolverCalls;
    std::vector<uint32_t> collectedAssumptions;

    if (_cascade->_multipleCascade != nullptr)
        _bucketAssumptions = _cascade->_multipleCascade->GetAssumptions();

    if (_setting->verbosity > 6)
        std::cout << std::endl << __PRETTY_FUNCTION__ << std::endl;

    DumpSolveInformation(true, localCalc);

    // solve at least once without assumptions to know if it is solvable at all.
    if (solverCallsStartSolveBucket == 0)
	  currentresult = _antom->Solve();

    while ( currentresult == ANTOM_SAT )
    {
        if (_isLastBucket)
            currSatWeight = CalculateSatWeight(localCalc);


        // all entries of last bucket are satisfied -> start with minimizing tare
        if (actualPos == size()-1)
            break;
        if(_isLastBucket)
            CalcBoundaries(actualPos, localCalc);

        lastPos = actualPos;

        if (_setting->interimResult != NOINTERIMRESULT)
            actualPos++;
        else
            actualPos = CalcNextPositionToSolve(currSatWeight, lastPos);

        if (lastPos == 0 && actualPos == 1)
        {
            lastPos = actualPos == 0 ? 0 : actualPos - 1;
            SetAsUnitClause(lastPos, currentresult, onlyWithAssumptions);
        }

        collectedAssumptions = GetAssumptions(actualPos);


        if( _setting->maxInprocess && _antom->Preprocess(INPROCESS) != ANTOM_UNKNOWN)
            break;

        if (_setting->verbosity > 2)
            std::cout << std::setw(50) << "TRY TO SOLVE POSITION: " << actualPos << std::endl;
        currentresult = _antom->Solve(collectedAssumptions);
        if (_setting->verbosity > 3)
        {
            if (currentresult == ANTOM_SAT)
                std::cout << std::endl << "ANTOM_SAT" << std::endl;
            else if (currentresult == ANTOM_UNSAT)
                std::cout << std::endl << "ANTOM_UNSAT" << std::endl;
            else if (currentresult == ANTOM_UNKNOWN)
                std::cout << std::endl << "ANTOM_UNKNOWN" << std::endl;
            else
                std::cout << std::endl << currentresult << std::endl;
        }

        SetAsUnitClause(actualPos, currentresult, onlyWithAssumptions);

        DumpSolveInformation(false, localCalc, currentresult, lastPos, actualPos);
//        std::cout << currentresult << std::endl;
    }

    actualPos = EvaluateResult(currentresult, actualPos, lastPos, onlyWithAssumptions);
    if      (currentresult == ANTOM_UNKNOWN ||            // case UNSAT, pos 0 couldn't be fulfilled
            (actualPos == 0 && _cascade->_estimatedWeightBoundaries[0] == 0 && currentresult == ANTOM_UNSAT))
    {
        if (localCalc)
            _localMaxPos = actualPos;
       return actualPos;
    }
    if(!localCalc)
        CalcBoundaries(actualPos, localCalc, true);

    std::cout << "c #solver calls bucket...: " << _antom->_satSolverCalls - solverCallsStartSolveBucket << std::endl;
    if (localCalc)
        _localMaxPos = actualPos;
    return actualPos;
}

std::vector<uint32_t> Bucket::GetAssumptions(uint32_t actualPos)
{
    if (_setting->verbosity > 6)
        std::cout << std::endl << __PRETTY_FUNCTION__ << std::endl;

    std::vector<uint32_t> currentAssumptions = _bucketAssumptions;
    currentAssumptions.push_back(( _sorter->GetOrEncodeOutput(actualPos) << 1) ^ 1);

    if (_setting->verbosity < 4)
        return currentAssumptions;

    std::cout << std::setw(51) << "Assumptions: (";
    for(auto assumption : currentAssumptions)
    {
         std::cout << assumption << ", ";
    }
    std::cout << ")" << std::endl;

    return currentAssumptions;
}

uint32_t Bucket::CalcNextPositionToSolve(uint64_t satWeight, uint32_t lastPos)
{
    if (_setting->verbosity > 3)
        std::cout << std::endl << __PRETTY_FUNCTION__ << std::endl;

    if (!_isLastBucket)
        return lastPos + 1;

    // is -1 necessary?
//    uint32_t nextPos = ((static_cast<int64_t>(satWeight) - static_cast<int64_t>(_groundWeight) - 1) / static_cast<int64_t>(_multiplicator)) + 1 > 0 ? ((satWeight + _groundWeight - 1) / _multiplicator) + 1: 0;
    uint32_t nextPos = (satWeight + _cascade->_tareWeight - _groundWeight ) / _multiplicator > 0 ?
                ( satWeight + _cascade->_tareWeight - _groundWeight ) / _multiplicator :
                0;

    nextPos = static_cast<int32_t>(nextPos - _positionsCuttedAtBottom) >= 0 ? nextPos - _positionsCuttedAtBottom: 0;

    nextPos = ( nextPos <= lastPos ) ? lastPos + 1 : nextPos;

    if (nextPos >= size())
        nextPos = size() - 1;

    return nextPos;
}

void Bucket::CalcBoundaries(uint32_t actualPos, bool localCalc, bool finalCalculation)
{
    if (_setting->verbosity > 6)
        std::cout << std::endl << __PRETTY_FUNCTION__ << std::endl;

//    std::cout << "actualPos: " << actualPos << std::endl;
//    std::cout << "multiplic: " << _multiplicator << std::endl;
//    std::cout << "groundWei: " << _groundWeight << std::endl;
    // set boundaries only with knowledge of this cascade
    if (_cascade->_multipleCascade == nullptr || _solvingState == TARECASCADE)
    {
        _cascade->_estimatedWeightBoundaries[0] = (actualPos + _positionsCuttedAtBottom) * _multiplicator + _groundWeight + 1;

        if (finalCalculation)
            _cascade->_estimatedWeightBoundaries[1] = (actualPos + 1) * _multiplicator + _groundWeight;
        else
            _cascade->_estimatedWeightBoundaries[1] = size() * _multiplicator + _groundWeight;

//        if (actualPos == 0)
//        {
//            _cascade->_estimatedWeightBoundaries[0]--;
//            _cascade->_estimatedWeightBoundaries[1]--;
//        }

    }
    else if (_solvingState == MAINCASCADE || localCalc)
    // main Cascade last bucket.
    {
//        _cascade->_estimatedWeightBoundaries[0] =
//                ((static_cast<int64_t>(actualPos) + 1) * static_cast<int64_t>(_multiplicator) - static_cast<int64_t>(_sumSoftWeightsOfTareCascade) > 0) ?
//                    (actualPos + 1) * _multiplicator - _sumSoftWeightsOfTareCascade:
//                    0;
//        std::cout << "actualPos: " << actualPos << std::endl;

        _cascade->_estimatedWeightBoundaries[0] = (actualPos + 1 +_positionsCuttedAtBottom) * _multiplicator - _sumSoftWeightsOfTareCascade;

        if (finalCalculation)
            _cascade->_estimatedWeightBoundaries[1] = _cascade->_estimatedWeightBoundaries[0] + _upperWeightBoundOfTareCascade;
        else
            _cascade->_estimatedWeightBoundaries[1] = (size() + _positionsCuttedAtBottom) * _multiplicator;
    }

    if (_setting->verbosity < 2)
        return;

//    std::cout << "c size of last bucket!!: " << size() << std::endl;
//    std::cout << "c mult!!: " << _multiplicator << std::endl;
    std::cout << std::setw(52) << "weight boundaries: ( " << _cascade->_estimatedWeightBoundaries[0] << " / " << _cascade->_estimatedWeightBoundaries[1] << " )" << std::endl;
//    assert(_cascade->_estimatedWeightBoundaries[0] <= static_cast<int64_t>(satWeight) || satWeight == 0);
//    assert(static_cast<int64_t>(satWeight) <= _cascade->_estimatedWeightBoundaries[1]);
}

uint64_t Bucket::CalculateSatWeight(bool localCalc)
{
    if (_setting->verbosity > 6)
        std::cout << std::endl << __PRETTY_FUNCTION__ << std::endl;

    // update overall sat weight
    // TODO: change to multipleCascade calculation if calculated for one multiple cascade...



//    std::cout << "LSATW: " << localSatWeight << std::endl;
//    std::cout << "GSATW: " << globalSatWeight << std::endl;

    if (localCalc)
    {
        uint64_t localSatWeight = _cascade->CountSatisfiedSoftClauses(this, _antom->Model());
        return localSatWeight;
    } else {
        uint64_t globalSatWeight = _antom->CalculateOverallOptimum(0, true);
        return globalSatWeight;
    }

}

void Bucket::SetAsUnitClause(uint32_t actualPos, uint32_t currentresult, bool onlyWithAssumptions)
{
    if (_setting->verbosity > 6)
        std::cout << std::endl << __PRETTY_FUNCTION__ << std::endl;

    bool negateLiteral;
    if (currentresult == ANTOM_SAT)
    {
//        std::cout << "ANTOM SAT" << std::endl;
        negateLiteral = true;
    }
    else if (currentresult == ANTOM_UNSAT)
    {
//        std::cout << "ANTOM UNSAT" << std::endl;
        negateLiteral = false;

        // encode this position with 01 mode
        // to efficiently set it true
        // can be very costly
        if (_setting->lastPos1)
            _sorter->GetOrEncodeOutput(actualPos, true);
    }
    else
    {
        _antom->_resultUnknown = true;
        return;
    }

    if (onlyWithAssumptions)
    {
        _bucketAssumptions.push_back(( _sorter->GetOrEncodeOutput(actualPos) << 1) ^ negateLiteral);
    }
    else
    {
#ifndef NDEBUG
        bool rst = _antom->AddUnit( ( _sorter->GetOrEncodeOutput(actualPos) << 1) ^ negateLiteral );
        assert(rst);
#else
        _antom->AddUnit( ( _sorter->GetOrEncodeOutput(actualPos) << 1) ^ negateLiteral );
#endif
    }


    if (_setting->verbosity < 3)
        return;

    uint32_t literal = ( _sorter->GetOrEncodeOutput(actualPos) << 1) ^ negateLiteral;

    if (onlyWithAssumptions)
        std::cout << std::setw(50) << "Set following Literal as assumption: " << _sorter->GetOrEncodeOutput(actualPos) << " | " << literal << std::endl;
    else
        std::cout << std::setw(50) << "Set following Literal as unit clause: " << _sorter->GetOrEncodeOutput(actualPos) << " | " << literal << std::endl;
}


void Bucket::DumpSolveInformation(bool head, bool localCalc, uint32_t currentresult, uint32_t lastPos, uint32_t actualPos)
{
    if (_setting->verbosity == 0)
        return;

    if (_setting->verbosity > 6)
        std::cout << std::endl << __PRETTY_FUNCTION__ << std::endl;

    if (head)
    {
        std::cout << std::setw(50) << "Bucket Multiplicator: " << std::setw(13) <<  _multiplicator << std::endl;
        //std::cout << std::endl << std::setw(30) << "         Bucket Size: " << std::setw(13) << size(false) << std::setw(47) << "Solve Last Bucket!" << std::endl;
        std::cout << std::setw(50) << "         Bucket Size: " << std::setw(13) << size(true) << std::endl;

        if (localCalc)
            std::cout << std::endl << std::setw(80) << "Only calculate bucket/cascade boundaries!" << std::endl;
        else
            std::cout << std::endl << std::setw(80) << "Calculate global weight boundaries!" << std::endl;

        if (_isLastBucket)
            std::cout << std::endl << std::setw(80) << "Solve Last Bucket: " << std::endl;
        else
            std::cout << std::endl << std::setw(80) << "Solve Bucket: " << std::endl;

        if (_setting->verbosity > 1)
            std::cout << "-------------------------------------------------------------------------------------------" << std::endl;

        return;
    }

    if (_setting->verbosity > 1)
    {
        std::cout << std::setw(52) << "SATweight Casc/Antom/#CascWeight/#AntomWeight: ( " << _cascade->_satWeight << " / " << _antom->_satWeight
                  << " / " << _cascade->_sumOfSoftWeights << " / " << _antom->_sumOfSoftWeights << " )" << std::endl;
        std::cout << std::setw(50) << "Last Position / Actual Position / Size: " << lastPos << " / " << actualPos << " / " << size() - 1 << std::endl;

        if (currentresult != ANTOM_SAT)
            return;

        std::cout << std::setw(90) << "ANTOM SAT" << std::endl;
        std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
    }
}

uint32_t Bucket::EvaluateResult(uint32_t currentresult, uint32_t actualPos, uint32_t lastPos, bool onlyWithAssumptions)
{
  if (currentresult == ANTOM_UNKNOWN || _control->ReachedLimits())
    {
        if (_setting->verbosity > 1)
            std::cout << std::setw(90) << "ANTOM UNKNOWN" << std::endl;

        actualPos = (actualPos > 0) ? actualPos - 1 : 0;
        //actualPos--;

        _antom->_resultUnknown = true;
        actualPos = lastPos;
    }
    else if (currentresult == ANTOM_UNSAT)
    {
        if (actualPos == 0)
        {
            _cascade->_estimatedWeightBoundaries[0]=0;
            _cascade->_estimatedWeightBoundaries[1]=0;
            return 0;
        }

        actualPos--;

        // Last satisfiable position!
        if (lastPos != actualPos)
            SetAsUnitClause(actualPos, ANTOM_SAT, onlyWithAssumptions);

        if (_setting->verbosity > 1)
            std::cout << std::setw(90) << "ANTOM UNSAT" << std::endl;

    }
    else if (currentresult == ANTOM_SAT)
    {
        if (_setting->verbosity > 1)
            std::cout << "Last Position could be solved!" << std::endl;
    }

    if (_setting->verbosity > 0)
        std::cout << "===========================================================================================" << std::endl;

    std::cout << "c last solvable position.: " << actualPos << "/" << size() - 1 << std::endl;

    return actualPos;
}

uint32_t Bucket::SolveTares(uint64_t diffEstimatedToCurrentSatWeight, uint32_t currentresult)
{
    std::vector<uint32_t> collectedAssumptions;

    //assert(_multiplicator * (_base - _tarePosition) >= _cascade->_estimatedWeightBoundaries[1] - _cascade->_estimatedWeightBoundaries[0]);

    // The solver shouldn't access this function if there is no tare in that bucket.
    assert(_tares.size() != 0);

    // all tares are set -> return
    if( _tarePosition == _base - 1 )
    {
        // process next Bucket;
        return currentresult;
    }


    // The actual result implies that the currentTare and maybe more tares can be set directly to TRUE
    if (diffEstimatedToCurrentSatWeight < _multiplicator)
    {
#ifndef NDEBUG
        assert(_tarePosition < _base - 1);
        bool rst = _antom->AddUnit(_tares[_tarePosition] << 1); assert(rst);
#else
        _antom->AddUnit(_tares[_tarePosition] << 1);
#endif

        if (_setting->verbosity > 1)
        {
            std::cout << std::setw(50) << "DiffToSatWeight < BucketMultiplicator: " << diffEstimatedToCurrentSatWeight << " < " << _multiplicator << std::endl;
            std::cout << std::setw(50) << "Set Tare: " << helper::Lit(_tares[_tarePosition] << 1) << std::endl;
            //std::cout << std::setw(50) << "Set: " << _tares[_tarePosition] * 2 << std::endl;
        }

        _tarePosition++;
        _cascade->_estimatedWeightBoundaries[0] += _multiplicator;

        diffEstimatedToCurrentSatWeight = _cascade->_estimatedWeightBoundaries[1] - _cascade->_satWeight;

        // continue in same bucket with next tare.
        return SolveTares(diffEstimatedToCurrentSatWeight, ANTOM_SAT);
    }
    // Corner Case!
    // If new result is already larger as maximal possible weight, we do not need to try,
    // -> the corresponding tare can be directly set to FALSE
    else if( (_cascade->_estimatedWeightBoundaries[1] - _multiplicator) > _antom->_sumOfSoftWeights ) //_cascade->_sumOfSoftWeights )
    {
#ifndef NDEBUG
        bool rst = _antom->AddUnit((_tares[_tarePosition] << 1)^1); assert(rst);
#else
        _antom->AddUnit((_tares[_tarePosition] << 1)^1);
#endif

        if (_setting->verbosity > 1)
        {
            std::cout << std::setw(50) << "Estimation larger than largest possible result." << std::endl;
            std::cout << std::setw(50) << "Set Tare: " << helper::Lit((_tares[_tarePosition] << 1)^1) << std::endl;
            //std::cout << std::setw(50) << "Set: " << _tares[_tarePosition]*2+1  << std::endl;
        }

        _tarePosition++;
        _cascade->_estimatedWeightBoundaries[1] -= _multiplicator;
        diffEstimatedToCurrentSatWeight = _cascade->_estimatedWeightBoundaries[1] - _antom->_satWeight;

        return SolveTares(diffEstimatedToCurrentSatWeight, ANTOM_UNSAT);
    }
    // At this point we have have found a bucket for which we try to solve the tare with TRUE
    else
    {
        collectedAssumptions.push_back(_tares[_tarePosition] << 1);
        // increase tare position
        _tarePosition++;
        if (_setting->verbosity > 1)
        {
            std::cout << std::setw(50) << "DiffToSatWeight >= BucketMultiplicator: " << diffEstimatedToCurrentSatWeight << " >= " << _multiplicator << std::endl;
            std::cout << std::setw(50) << "Try to set Tare: " << helper::Lit(collectedAssumptions.back()) << std::endl;
            //std::cout << std::setw(50) << "Collected assumptions back: " << collectedAssumptions.back() << std::endl;
        }
    }


    // Inprocessing
    if( _setting->maxInprocess )
    {
        // Do inprocessing
        currentresult = _antom->Preprocess(INPROCESS);

        // Have we found the optimum?
        if (currentresult != ANTOM_UNKNOWN)
        { return currentresult; }
    }

    // Solve next tare
    currentresult = _antom->Solve(collectedAssumptions);

    if (currentresult == ANTOM_SAT)
    {
        // save satisfied model!
        // TOBI: - only if actual Model is better than last one!?
        //_antom->_lastModel = _antom->Model();

#ifndef NDEBUG
        bool rst = _antom->AddUnit(collectedAssumptions.back()); assert(rst);
        assert(collectedAssumptions.size() == 1);
#else
        _antom->AddUnit(collectedAssumptions.back());
#endif
        _cascade->_estimatedWeightBoundaries[0] += _multiplicator;
        collectedAssumptions.pop_back();

        _cascade->CountSatisfiedSoftClauses(NULL, _antom->Model());
        _antom->CalculateOverallOptimum(0, true);
//        _antom->CountSatisfiedSoftClauses(NULL, _antom->Model());
//        _antom->CalculateOverallOptimum(_cascade->_satWeight, true);
        diffEstimatedToCurrentSatWeight = _cascade->_estimatedWeightBoundaries[1] - _antom->_satWeight;

        if (_setting->verbosity > 2)
        {
            if (_base > 2)
            {
                std::cout << std::setw(50) << "Actual tare position: " << _tarePosition - 1 << std::endl;
            }
            std::cout << std::setw(50) << "Calculated SAT Weight: " << _cascade->_satWeight << std::endl;
        }

        // continue in same bucket with next tare.
        uint32_t nextTare = SolveTares(diffEstimatedToCurrentSatWeight, ANTOM_SAT);

        return nextTare;

    }
    else if (currentresult == ANTOM_UNSAT)
    {
        if (_setting->verbosity > 1)
        {
            std::cout << std::setw(50) << "Set Tare: " << helper::Lit((collectedAssumptions.back()^1)) << std::endl;
        }
#ifndef NDEBUG
        bool rst = _antom->AddUnit( collectedAssumptions.back()^1 ); assert(rst);
#else
        _antom->AddUnit( collectedAssumptions.back()^1 );
#endif
        uint32_t remainingTares = _base - _tarePosition;
        //std::cout << "base: " << _base << "  _tarePosition: " << _tarePosition << "  remainingTares: " << remainingTares << std::endl;
        _cascade->_estimatedWeightBoundaries[1] -= ( remainingTares * _multiplicator );
    }
    return currentresult;
}


std::vector<uint32_t> Bucket::GetEveryNthOutput(uint32_t n)
{
    assert(_encoded);
    std::vector< uint32_t > nthOutputVector;
    for (uint32_t i = 0; i < _outputs->size(); i++)
    {
        if (i % n == n - 1)
        {
            nthOutputVector.push_back((*_outputs)[i]);
        }
        // Since we merge the outputs of this sorter next, we can deactive the don't touch status
        _antom->SetDontTouch((*_outputs)[i], false);
    }
    // The tares should be still don't touch
    for( uint32_t i = 0; i != _tares.size(); ++i )
    {
        _antom->SetDontTouch(_tares[i], true);
    }
    return nthOutputVector;
}

void Bucket::EstimateNumberOfClauses(bool top)
{
    if (top)
    {
        _binaryTopClEstimated = static_cast<uint32_t>( size() * log2 ( size() ) );
        _ternaryTopClEstimated = static_cast<uint32_t>( pow( size(), 2) / 2 - ( size() / 2 ) );
    }
    else
    {
        _binaryBottomClEstimated = static_cast<uint32_t>( size() );
        _ternaryBottomClEstimated = static_cast<uint32_t>( ( size() - _topEntries ) * _topEntries );
    }
}

void Bucket::CalculateNumberOfClauses(bool top, bool setCurrent, bool calcEstimation)
{
    if (top && setCurrent)
    {
        _binaryTopCl = _antom->CurrentBinaryClauses();
        _ternaryTopCl = _antom->CurrentTernaryClauses();
    }
    else if (!top && setCurrent)
    {
        _binaryBottomCl = _antom->CurrentBinaryClauses();
        _ternaryBottomCl = _antom->CurrentTernaryClauses();
    }
    else if (top && !setCurrent)
    {
        _topEntries = size();
        _binaryTopCl = _antom->CurrentBinaryClauses() - _binaryTopCl;
        _ternaryTopCl = _antom->CurrentTernaryClauses() - _ternaryTopCl;
    }
    else if (!top && !setCurrent)
    {
        _binaryBottomCl = _antom->CurrentBinaryClauses() - _binaryBottomCl;
        _ternaryBottomCl = _antom->CurrentTernaryClauses() - _ternaryBottomCl;
    }

    if (calcEstimation)
    {
        EstimateNumberOfClauses(top);
    }
}


}

