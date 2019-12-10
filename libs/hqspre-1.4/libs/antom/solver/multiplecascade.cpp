/********************************************************************************************
multiplecascade.cpp -- Copyright (c) 2017, Tobias Paxian

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

#include "multiplecascade.h"
#include "antom.h"
#include "bucket.h"
#include "cascade.h"
#include "softclausenodes.h"
#include "sorter.h"
#include "softclause.h"
#include "totalizerencodetree.h"
//#include "timemeasurement.h"
//#include "timevariables.h"


namespace antom
{

// Constructor
MultipleCascade::MultipleCascade(Antom* antom, bool onlyByTares, bool solveTareCascadeOnlyByTares, uint32_t cascDiv, uint32_t maxBucketSize, uint32_t nOfCasc, InterimResult interimResult, uint64_t sumOfSoftWeights):
    _antom(antom),
	_setting(antom->_antomSetting),
    _onlyByTares(onlyByTares),
    _solveTareCascadeOnlyByTares(solveTareCascadeOnlyByTares),
    _base(_setting->base),
	_satWeight(0),
	_highestMultiplicator(0),
	_tareSoftClauses(),
	_partedSoftClauses(),
	_cascades(),
	_cascadesMM(),
	_mainCascade(NULL),
	_tareCascade(NULL),
    _tareWeight(0),
	_tareWeightWithoutAddTs(0),
    _tareWeightOfAddTs(0),
    _maxMainCascadePosition(0),
    _maxTareCascadePosition(0),
    _cascadeDivider(cascDiv),
    _maxBucketSize(maxBucketSize),
    _numberOfCascades(nOfCasc),
	_tareWeightUNSAT(0),
	_tareWeightSAT(0),
    _interimResult(interimResult),
    _sumOfSoftWeights(sumOfSoftWeights),
    _solveTareCascade(false),
    _solveLikeNormalCascade(false),
	_tareAssumptions()
{
    assert(antom != NULL);
}

bool MultipleCascade::DivideAndConnectCascades(MultipleCascadeDivideStrategy divideStrategy, std::vector<SoftClause*>* softClauses)
{
    if (_setting->verbosity > 3)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    // if there is only one cascade, then solve it with the normal cascade mode!
    if (!DivideSoftClauseVector(divideStrategy, softClauses))
        return false;

    FillSubCascades(_partedSoftClauses);

    if (!_onlyByTares)
        ConnectSubCascades();

    CreateSCVectorFromCascades(_cascades);

    // Add additional tares and SCs due to tareWeight!
    if (!_onlyByTares)
    {
//        AddCompensationTaresToLastBucket();

        _mainCascade->DumpBucketStructure(true, 3);
        _mainCascade->Encode();
//        std::cout << "SIZE OF  M A I N T R E E : " << _mainCascade->_structure.back()->_sorter->_outputTree->_size << std::endl;

        if (_setting->createGraphFile != "")
            _mainCascade->_structure.back()->_sorter->_outputTree->DumpOutputTree(_setting->createGraphFile + std::to_string(_mainCascade->_structure.back()->size()) + "MainCascadeOutputs.tgf", true);
    }

    _tareCascade = new Cascade(_antom, this, _solveTareCascadeOnlyByTares);

    //solve tares like cascade
    _tareCascade->Fill(&_tareSoftClauses, NOPARTITION, ENCODEONLYIFNEEDED);

    _tareCascade->Encode();
    _tareCascade->DumpBucketStructure(true, 3);

    return true;
}

uint32_t MultipleCascade::Solve()
{

    if (_setting->verbosity > 3)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    _tareWeightUNSAT = 0;
    _tareWeightSAT = 0;
    uint32_t currentresult(ANTOM_UNKNOWN);

    if (_setting->verbosity > 3)
    {
        std::cout << std::endl << "All Cascades Bucket Structure: " << std::endl;
        for (auto cascade : _cascadesMM)
            cascade.second->DumpBucketStructure(true, 3);
    }

    _tareWeightWithoutAddTs = _tareWeight - _tareWeightOfAddTs;

    if (!_onlyByTares)
    {
        uint64_t upperWeightBoundOfTC = _tareCascade->_structure.back()->size() * _tareCascade->_highestBucketMultiplicator;

        if (_solveTareCascadeOnlyByTares)
            upperWeightBoundOfTC = _tareCascade->_highestBucketMultiplicator;

        _mainCascade->_structure.back()->SetSolvingParameters(&_tareAssumptions,
                                                              upperWeightBoundOfTC,
                                                              _tareWeight,
                                                              0,
                                                              MAINCASCADE);
        if (_setting->verbosity > 3)
            std::cout << "Solve Main Cascade:" << std::endl;

        _maxMainCascadePosition = _mainCascade->_structure.back()->SolveBucketReturnMaxPosition(false, false);
        if (_setting->createGraphFile != "")
            _mainCascade->_structure.back()->_sorter->_outputTree->DumpOutputTree(_setting->createGraphFile + std::to_string(_maxMainCascadePosition) + "MainCascadeAfterSolvingOutputs.tgf", true);

        if (_antom->_resultUnknown)
            return ANTOM_UNKNOWN;

        _estimatedWeightBoundaries[0] = _mainCascade->_estimatedWeightBoundaries[0];
        _estimatedWeightBoundaries[1] = _mainCascade->_estimatedWeightBoundaries[1];
    } else {
        _estimatedWeightBoundaries[0] = 0 - _tareWeightWithoutAddTs;
        _estimatedWeightBoundaries[1] = 0 - _tareWeightWithoutAddTs;
    }

    _solveTareCascade = true;

     if (!_solveTareCascadeOnlyByTares)
    {
         _tareCascade->_structure.back()->SetSolvingParameters(&_tareAssumptions,
                                                               _tareCascade->_structure.back()->size() * _highestMultiplicator,
                                                               _tareWeight,
                                                               _estimatedWeightBoundaries[0],
                                                               TARECASCADE);

         if (_setting->verbosity > 3)
             std::cout << "Solve Tare Cascade:" << std::endl;

         _maxTareCascadePosition = _tareCascade->_structure.back()->SolveBucketReturnMaxPosition(false, false);
         if (_antom->_resultUnknown)
             return ANTOM_UNKNOWN;
    } else {
         _tareCascade->_estimatedWeightBoundaries[0] += _estimatedWeightBoundaries[0];
         if (_onlyByTares)
            _tareCascade->_estimatedWeightBoundaries[1] += _estimatedWeightBoundaries[1];
         else
            _tareCascade->_estimatedWeightBoundaries[1] += _estimatedWeightBoundaries[0];
    }

    if (_interimResult != NOINTERIMRESULT)
    {
//        std::cout << "InterimResult: " << _interimResult << std::endl;
//        std::cout << "_antomFT: " << _antom->_featureTest << std::endl;
        currentresult = _tareCascade->SolveAllTares();
    }
    else
    {
        currentresult = _tareCascade->SolveTares();
    }

    if (currentresult == ANTOM_UNKNOWN)
        return ANTOM_UNKNOWN;

    _antom->CalculateOverallOptimum(_satWeight, true);

    DumpFinalValues(3);

    return currentresult;
}

void MultipleCascade::DumpFinalValues(uint16_t verbosity)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    std::cout << "c Max weight.............: " << _antom->_satWeight << std::endl;
    std::cout << "c Est. Bd. after TareCasc: " << _tareCascade->_estimatedWeightBoundaries[0] << " / " << _tareCascade->_estimatedWeightBoundaries[1] << std::endl;
    std::cout << "c Diff. of TC boundaries.: " << _tareCascade->_estimatedWeightBoundaries[1] - _tareCascade->_estimatedWeightBoundaries[0] << std::endl;

    if  (_setting->verbosity < verbosity)
        return;

    std::cout << "c Tare Weights...........: " << _tareWeight << std::endl;
    std::cout << "c Additional Tare Weight.: " << _tareWeightOfAddTs << std::endl;
    std::cout << "c Tare Weights wo add....: " << _tareWeight - _tareWeightOfAddTs << std::endl;

    uint64_t actualTareWeight = _tareCascade->CountSatisfiedSoftClauses(NULL, _antom->_lastModel);
    uint64_t maxTareWeight(0);
    maxTareWeight = actualTareWeight > maxTareWeight ? actualTareWeight : maxTareWeight;
    std::cout << "c TareWeight being one...: " << maxTareWeight << std::endl;
    std::cout << "c TareWeight being zero..: " << _tareWeight - maxTareWeight << std::endl;
}

bool MultipleCascade::DivideSoftClauseVector(MultipleCascadeDivideStrategy divideStrategy, std::vector<SoftClause*>* softClauses)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

//    std::cout << "c #Cascades..............: " << _numberOfCascades << std::endl;
//    std::cout << "c softClauses->size()....: " << softClauses->size() << std::endl;
//    std::cout << "c _cascadeDivider........: " << _cascadeDivider << std::endl;
    // DO HERE THE SEPARATION STUFF
    if (_numberOfCascades <= 1 && (_cascadeDivider >= softClauses->size() || _cascadeDivider < 2) )
    {
        _setting->cascadeDivider = 0;
        return false;
    }
    if (_numberOfCascades == 0)
        _numberOfCascades = static_cast<uint32_t>(ceil(static_cast<double>(softClauses->size()) / static_cast<double>(_cascadeDivider)));
    uint32_t sizeOfOnePart = static_cast<uint32_t>(ceil(static_cast<double>(softClauses->size()) / static_cast<double>(_numberOfCascades)));
    // TOBI: maybe better to just
    if (_maxBucketSize <= 2)
        _maxBucketSize = static_cast<uint32_t>(-1);


    if (_cascadeDivider <= 2)
    {
        _cascadeDivider = static_cast<uint32_t>(ceil(static_cast<double>(softClauses->size()) / static_cast<double>(_numberOfCascades)));
        _setting->cascadeDivider = _cascadeDivider;
    }

    _mainCascade = new Cascade(_antom, this, false);

    _highestMultiplicator = 0;

    if (_setting->verbosity > 3)
    {
        std::cout << "c cascadeDivider: " << _cascadeDivider << std::endl;
        std::cout << "c softClauses->size(): " << softClauses->size() << "    maxBucketSize for additional buckets: " << _maxBucketSize << std::endl;
        std::cout << "c howManyCascades: " << _numberOfCascades << std::endl;
        std::cout << "c sizeOfOnePart: " << sizeOfOnePart << std::endl << std::endl;

        std::cout << "c All Softclauses, size: " << softClauses->size() << std::endl;
        for (auto clause = softClauses->begin(); clause != softClauses->end(); ++clause)
            std::cout << (*clause)->weight << ' ';
        std::cout << std::endl << std::endl;
    }

    std::cout << "c #Cascades..............: " << _numberOfCascades << std::endl;

    std::vector<SoftClause*> sClauses;
    sClauses = *softClauses;

    if (divideStrategy == SORTEDNUMBEROFSOFTCLAUSES)
    {
        if (_setting->verbosity > 3)
            std::cout << std::endl << "c Sort!" << std::endl;
        std::stable_sort(sClauses.begin(), sClauses.end(), SoftClause::bigger );

        uint64_t sumOfWeights = 0;
        for (unsigned long ind = sClauses.size(); ind > 0; ind--)
        {
            if (sumOfWeights < sClauses[ind - 1]->weight && sumOfWeights > 0)
            {
                std::cout << "c Number Of Processed SCs: " << sClauses.size() - ind + 1<< std::endl;
                std::cout << "c Sum of weights tillthat: " << sumOfWeights << std::endl;
                std::cout << "c Next higher SoftWeight.: " << sClauses[ind - 1]->weight << std::endl;
            }
            if (ind > 1 && sClauses[ind - 1]->weight * 10 < sClauses[ind - 2]->weight)
            {
                std::cout << "c great diff in SC factor: " << sClauses[ind - 2]->weight / sClauses[ind - 1]->weight << std::endl;
                std::cout << "c actual / next weight...: " << sClauses[ind - 1]->weight << "/" << sClauses[ind - 2]->weight << std::endl;
            }

//                std::cout << ind << ", " << sumOfWeights << ", " << sClauses[ind - 1]->weight << std::endl;
            sumOfWeights += sClauses[ind - 1]->weight;
        }
    }
    else if (divideStrategy == RANDOMNUMBEROFSOFTCLAUSES)
    {
//          not really random - generates always the same sequence!
        std::cout << std::endl << "c Shuffle SoftClauses!" << std::endl;
        std::random_shuffle(sClauses.begin(), sClauses.end());
    }
    else if (divideStrategy == SOFTCLAUSESINORDER)
    {
        std::cout << std::endl << "c Don't change order of SC's!" << std::endl;
    }

    if (_setting->verbosity > 3)
    {
        std::cout << "c All SCs size: " << softClauses->size() << std::endl;
        for (auto clause = sClauses.begin(); clause != sClauses.end(); ++clause)
            std::cout << (*clause)->weight << ' ';
        std::cout << std::endl;
    }


    // divide SoftClauses
    for (uint32_t ind = 0; ind < _numberOfCascades; ind++)
    {
        uint32_t startPoint = ind * sizeOfOnePart;
        //uint32_t endPoint = ind < _numberOfCascades-1 ? (ind+1) * sizeOfOnePart : sClauses.size();
        uint32_t endPoint = (ind+1) * sizeOfOnePart < sClauses.size() ? (ind+1) * sizeOfOnePart : sClauses.size();

        if (endPoint == startPoint)
            break;

        std::vector<SoftClause*> softClausesPart(sClauses.begin()+startPoint, sClauses.begin()+endPoint);
        _partedSoftClauses.push_back(softClausesPart);

        if (_setting->verbosity > 3)
        {
            std::cout << std::endl << "startPoint: " << startPoint << "    endPoint: " << endPoint << std::endl;
            std::cout << "PartedSoftClauses, Part " << ind << ", size: " << _partedSoftClauses.back().size() << std::endl;
            for (auto clause = _partedSoftClauses.back().begin(); clause != _partedSoftClauses.back().end(); ++clause)
                std::cout << (*clause)->weight << ' ';
            std::cout << std::endl;
        }
        if (endPoint == sClauses.size())
            break;
    }
    return true;
}

void MultipleCascade::FillSubCascades(std::vector<std::vector<SoftClause*>> partedSoftClauses)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    uint64_t multiplicator = 0;
    for (auto softClauseVector: partedSoftClauses)
    {
        Cascade* tmpCascade = new Cascade(_antom, this, _onlyByTares);
        tmpCascade->Fill(&softClauseVector, _setting->partitionStrategy, ENCODEONLYIFNEEDED);

        if (_onlyByTares)
            _tareWeightOfAddTs += tmpCascade->_highestBucketMultiplicator;

        _cascades.push_back(tmpCascade);
        _highestMultiplicator = _highestMultiplicator < tmpCascade->_highestBucketMultiplicator ? tmpCascade->_highestBucketMultiplicator : _highestMultiplicator;


        multiplicator = multiplicator == 0 ? tmpCascade->_highestBucketMultiplicator : multiplicator;
        _cascadesMM.insert(std::make_pair(tmpCascade->_highestBucketMultiplicator, tmpCascade));
    }

    if (_setting->verbosity > 3)
    {
        std::cout << "Highest Multiplicator: " << _highestMultiplicator << std::endl << std::endl;
    }
}

void MultipleCascade::ConnectLeastWeightedSmallestCascades()
{
    while(_cascadesMM.size() != 1)
    {
        if (_setting->verbosity > 3)
            std::cout << std::endl << "_cascadesMM.size(): " << _cascadesMM.size() << std::endl;
        assert(_cascadesMM.size() > 0);
        Cascade* firstCascade = _cascadesMM.begin()->second;
        Cascade* secondCascade = (++_cascadesMM.begin())->second;
        // the first two elements have the same bucket highest multiplicator
        if (_cascadesMM.begin()->first == (++_cascadesMM.begin())->first)
        {
            if (_setting->verbosity > 3)
                std::cout << "Same highest bucket multiplicator." << std::endl;
            // sum of sizes is too big, make them smaller!
            if (_cascadesMM.begin()->second->_structure.back()->size(true) + (++_cascadesMM.begin())->second->_structure.back()->size(true) > _maxBucketSize)
            {
                if (_setting->verbosity > 3)
                    std::cout << "Sizes are too big! Size first cascade: " << _cascadesMM.begin()->second->_structure.back()->size(true)
                              << "  Size second Cascade: " << (++_cascadesMM.begin())->second->_structure.back()->size(true) << std::endl;
                if (firstCascade->_structure.back()->size(true) > secondCascade->_structure.back()->size(true))
                {
                    firstCascade->AddNewBucketsTillSizeBoundary(_maxBucketSize / 2, false, true);
                    _cascadesMM.erase(_cascadesMM.begin());
                    _cascadesMM.insert(std::make_pair(firstCascade->_highestBucketMultiplicator,firstCascade));
                } else {
                    secondCascade->AddNewBucketsTillSizeBoundary(_maxBucketSize / 2, false, true);
                    _cascadesMM.erase((++_cascadesMM.begin()));
                    _cascadesMM.insert(std::make_pair(secondCascade->_highestBucketMultiplicator,secondCascade));
                }
                assert(firstCascade->_highestBucketMultiplicator != secondCascade->_highestBucketMultiplicator);
                continue;
            }
            // sum of sizes is small enough, add the first cascades last bucket to the second cascades last bucket _subbuckets
            // erase then the _cascadeMM.begin() entry!
            if (_setting->verbosity > 3)
                std::cout << "Sum of sizes are small enough push firstCascade into secondCascade." << std::endl;
            if (firstCascade->_structure.back()->_tares.empty())
                firstCascade->AddTare(firstCascade->_structure.size()-1);


            if (_setting->interimResult == CUTATTOP)
            {

                if (secondCascade->_structure.back()->_tares.empty())
                    secondCascade->AddTare(secondCascade->_structure.size()-1);

                _tareAssumptions.push_back(firstCascade->_structure.back()->_tares[0] << 1);
                _tareWeightOfAddTs += firstCascade->_highestBucketMultiplicator;

                /*uint32_t maxPos = */firstCascade->CutMaxPos();
                secondCascade->CutMaxPos();
                //secondCascade->_structure.back()->_localMaxPos = static_cast<uint32_t>(-1);
            } else
            {
                _tareAssumptions.push_back(firstCascade->_structure.back()->_tares[0] << 1);
                _tareWeightOfAddTs += firstCascade->_highestBucketMultiplicator;
            }

            secondCascade->_structure.back()->_subBuckets.push_back(firstCascade->_structure.back());
            _cascadesMM.erase(_cascadesMM.begin());

            continue;
        } else {
//            std::cout << "Bucket multiplicator differs by more than one bucket." << std::endl;
            if (_setting->verbosity > 3)
                std::cout << "Bucket multiplicator differs." << std::endl;
            if (!firstCascade->AddNewBucketsTillMultiplicatorMatches(secondCascade->_highestBucketMultiplicator, false, true))
            {
                //this should never happen!
                assert(true);
            }
            _cascadesMM.erase(_cascadesMM.begin());
            _cascadesMM.insert(std::make_pair(firstCascade->_highestBucketMultiplicator,firstCascade));
            continue;
        }
        assert(true);
    }
}

void MultipleCascade::ConnectHighestWeightedCascades()
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    uint64_t upperBoundAllLowerCascades = 0;

    _cascades[0]->AddTare(_cascades[0]->_structure.size()-1);

    //cut all buckets at top at the beginning - and get a better upper weight bound for the smaller cascades!
    for(int32_t i =_cascades.size() - 1; i > 0; i--)
    {
        std::cout << "Next Cascade to process: " << i << std::endl;
        uint32_t maxPos(0);
        uint64_t actualHighestMult(0);
        uint64_t actualMaxWeight(0);

        if(_cascades[i]->_structure.back()->_tares.empty())
            _cascades[i]->AddTare(_cascades[i]->_structure.size() - 1);

        if (_setting->interimResult == CUTATTOP || _setting->interimResult == CUTBOTH)
        {
            maxPos = _cascades[i]->CutMaxPos(true);
            actualHighestMult = _cascades[i]->_highestBucketMultiplicator;
            actualMaxWeight = (maxPos + 1) * actualHighestMult;
            std::cout << "                  ActualMaxWeight: " << actualMaxWeight << std::endl;
        }

        if (_cascades[i]->_sumOfSoftWeights + _cascades[i]->_highestBucketMultiplicator - 1 < actualMaxWeight || actualMaxWeight == 0)
        {
            // sum of weights + sum of tares
            actualMaxWeight = _cascades[i]->_sumOfSoftWeights + _cascades[i]->_highestBucketMultiplicator - 1;
        }


        std::cout << "                   ActualMaxWeight: " << actualMaxWeight << std::endl;
        std::cout << "                Weight Of all SCs: " << _cascades[i]->_sumOfSoftWeights << std::endl;
        std::cout << "_upperWeightBoundAllLowerCascades: " << upperBoundAllLowerCascades << std::endl;
        std::cout << "    Weight of all Higher Cascades: " << _sumOfSoftWeights - upperBoundAllLowerCascades << std::endl;
        std::cout << "           Weight of all Cascades: " << _sumOfSoftWeights << std::endl;
        _cascades[i]->_upperWeightBoundAllLowerCascades = upperBoundAllLowerCascades;

        upperBoundAllLowerCascades += actualMaxWeight;
        assert(upperBoundAllLowerCascades > 0);
    }

    _cascades[0]->_upperWeightBoundAllLowerCascades = upperBoundAllLowerCascades;

#ifndef NDEBUG
    bool test = _cascades[0]->AddNewBucketsTillMultiplicatorMatches(_highestMultiplicator, false, false); assert(test);
#else
	_cascades[0]->AddNewBucketsTillMultiplicatorMatches(_highestMultiplicator, false, false);
#endif

    for(uint32_t ind=1; ind < _cascades.size(); ind++)
    {

        std::cout << "new iteration: " << ind << std::endl;

        if (_setting->interimResult == CUTATTOP || _setting->interimResult == CUTBOTH)
		  {
            _cascades[0]->CutMaxPos(true);
		  }
        std::cout << "AFTER CUTMAXPOS" << std::endl;
        if (_setting->interimResult == CUTATBOTTOM)
		  {
            _cascades[0]->CutMinPos(true);
		  }
        if (_setting->interimResult == CUTBOTH)
		  {
            _cascades[0]->CutMinPos(false);
		  }

        std::cout << "AFTER CUTMINPOS" << std::endl;
//        if (_setting->interimResult != NOINTERIMRESULT)
//            _cascades[0]->_structure.back()->_localMaxPos = static_cast<uint32_t>(-1);

        if (_setting->verbosity > 3)
		  {
            std::cout << "Bucket multiplicator differs." << std::endl;
		  }

        std::cout << "_highestMultiplicator: " << _highestMultiplicator << std::endl;

#ifndef NDEBUG
        bool test = _cascades[ind]->AddNewBucketsTillMultiplicatorMatches(_highestMultiplicator, false, false); assert(test);
#else
		_cascades[ind]->AddNewBucketsTillMultiplicatorMatches(_highestMultiplicator, false, false);
#endif

        if(_cascades[ind]->_structure.back()->_tares.empty())
		  {
            _cascades[ind]->AddTare(_cascades[ind]->_structure.size() - 1);
		  }

         _tareAssumptions.push_back(_cascades[ind]->_structure.back()->_tares[0] << 1);

        if(_cascades[0]->_structure.back()->_tares.empty())
		  {
            _cascades[0]->AddTare(_cascades[0]->_structure.size() - 1);
		  }

        std::cout << "ConnectCascades: 0 + " << ind << std::endl;
        std::cout << "Highest Multipl: " << _cascades[0]->_highestBucketMultiplicator << ", " << _cascades[ind]->_highestBucketMultiplicator << std::endl;
        std::cout << "          Sizes: " << _cascades[0]->_structure.back()->size(true) << ", " << _cascades[ind]->_structure.back()->size(true) << std::endl;
        std::cout << "          Sizes: " << _cascades[0]->_structure.back()->size(false) << ", " << _cascades[ind]->_structure.back()->size(false) << std::endl;
        //secondCascade->_structure.back()->_subBuckets.push_back(firstCascade->_structure.back());
        _cascades[0]->_upperWeightBoundAllLowerCascades = _cascades[ind]->_upperWeightBoundAllLowerCascades;
        _cascades[0]->_structure.back()->_subBuckets.push_back(_cascades[ind]->_structure.back());
        _cascades[0]->_sumOfSoftWeights += _cascades[ind]->_sumOfSoftWeights;
    }
//    _cascades[0]->_structure.back()->_localMaxPos = 0;
}

void MultipleCascade::ConnectSubCascades()
{
    if (_setting->verbosity > 3)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    // put the cascades together again!!
    // most useful with cutting at bottom
    if (_setting->interimResult == CUTATBOTTOM || _setting->interimResult == CUTBOTH)
//    if (true)
    {
        ConnectHighestWeightedCascades();

        _highestMultiplicator = _cascades[0]->_highestBucketMultiplicator;
        std::cout << "c Highest multiplicator..: " << _highestMultiplicator << std::endl;

        _mainCascade->_structure.push_back(_cascades[0]->_structure.back());
        _mainCascade->_structure.back()->_cascade = _mainCascade;
        if (_mainCascade->_structure.back()->_tares.empty())
            _mainCascade->AddTare(_mainCascade->_structure.size()-1);

    }
    else
    {
        ConnectLeastWeightedSmallestCascades();


        _highestMultiplicator = _cascadesMM.begin()->first;
        std::cout << "c Highest multiplicator..: " << _highestMultiplicator << std::endl;

        _mainCascade->_structure.push_back(_cascadesMM.begin()->second->_structure.back());
        _mainCascade->_structure.back()->_cascade = _mainCascade;

//        _tareAssumptions.push_back(_mainCascade->_structure.back()->_tares[0] << 1);
//        _mainCascade->_highestBucketMultiplicator = _cascadesMM.begin()->first;
//        _tareWeightOfAddTs += _mainCascade->_highestBucketMultiplicator;
    }
    _tareAssumptions.push_back(_mainCascade->_structure.back()->_tares[0] << 1);
    _mainCascade->_highestBucketMultiplicator = _highestMultiplicator;
    _tareWeightOfAddTs += _mainCascade->_highestBucketMultiplicator;
    _mainCascade->_structure.back()->_isLastBucket = true;

}

std::vector<uint32_t> MultipleCascade::GetAssumptions()
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    if (_solveTareCascade)
        return {};
    return _tareAssumptions;
}

void MultipleCascade::CreateSCVectorFromCascades(std::vector<Cascade*> cascades)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    for (auto casc : cascades)
    {
        std::vector<std::vector< uint32_t >> tares;
        tares = casc->GetTareVector();

        CreateSCVectorFromTares(&tares, casc->_base);
    }
}

void MultipleCascade::CreateSCVectorFromTares(std::vector<std::vector< uint32_t >>* tares, uint32_t base)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;
    //iterate through all tares
    // i is the weight 2^i
    for (uint32_t i=0; i < tares->size(); i++)
    {
        for (uint32_t j=0; j < (*tares)[i].size(); j++)
        {
            uint64_t weight = static_cast<uint64_t>(pow(base,i));
            _tareSoftClauses.push_back(CreateSCForTare((*tares)[i][j], weight));
//            if (weight == _highestMultiplicator)
//                _tareAssumptions.push_back((*tares)[i][j] << 1);
            _tareWeight += weight;
        }
    }
}

SoftClause* MultipleCascade::CreateSCForTare(uint32_t tare, uint64_t weight)
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::vector<uint32_t> clause;
    uint32_t relaxLit = _antom->NewVariable();
    _antom->SetDontTouch(relaxLit);
//    _antom->SetDontTouch(tare);
    clause.push_back(tare << 1);
    SoftClause* softClause;
    softClause = new SoftClause(relaxLit << 1, clause, weight);
    clause.push_back(relaxLit << 1);
	_antom->AddClause(clause);
    if (_setting->verbosity > 3)
        std::cout << "SC-weight " << softClause->weight << ", tare: " << tare << ", relaxLit: " << relaxLit << std::endl;
    return softClause;
}

void MultipleCascade::AddSolutionAsUnitClauses()
{
    if (_setting->verbosity > 6)
        std::cout << __PRETTY_FUNCTION__ << std::endl;

    //frb-10-6-4.wcnf
    //int mySolution[] = {-1, -2, -3, -4, -5, 6, -7, 8, -9, -10, -11, -12, -13, 14, -15, -16, -17, -18, -19, -20, 21, -22, -23, -24, -25, -26, -27, -28, -29, 30, -31, -32, -33, -34, -35, 36, 37, -38, -39, -40, -41, -42, -43, -44, -45, 46, -47, -48, -49, 50, -51, -52, -53, -54, -55, -56, -57, -58, -59, 60};

    //frb15-9-1
    //int mySolution[] = {-1, -2, -3, -4, 5, -6, -7, -8, -9, -10, -11, 12, -13, -14, -15, -16, -17, -18, -19, -20, -21, -22, 23, -24, -25, -26, -27, -28, -29, -30, -31, -32, -33, -34, -35, 36, -37, -38, -39, -40, -41, -42, -43, 44, -45, -46, -47, -48, -49, -50, 51, -52, -53, -54, -55, -56, -57, -58, -59, 60, -61, -62, -63, -64, -65, 66, -67, -68, -69, -70, -71, -72, -73, -74, -75, -76, -77, -78, -79, -80, 81, -82, -83, -84, -85, 86, -87, -88, -89, -90, -91, -92, -93, -94, -95, -96, -97, -98, 99, -100, -101, -102, -103, -104, -105, 106, -107, -108, -109, -110, 111, -112, -113, -114, -115, -116, -117, -118, -119, -120, -121, -122, -123, -124, 125, -126, 127, -128, -129, -130, -131, -132, -133, -134, -135};

    //cat_paths_60_150_0002.txt
    //int mySolution[] = {-1, -2, -3, 4, -5, 6, -7, -8, -9, -10, -11, -12, 13, -14, -15, -16, -17, -18, 19, -20, -21, -22, -23, -24, 25, 26, -27, -28, -29, -30, -31, 32, -33, -34, 35, -36, -37, -38, -39, 40, -41, -42, 43, -44, -45, -46, -47, -48, 49, -50, 51, -52, -53, -54, -55, -56, 57, -58, -59, 60, -61, -62, 63, -64, -65, -66, -67, -68, -69, -70, -71, -72, -73, -74, -75, -76, -77, 78, -79, -80, 81, -82, 83, -84, 85, -86, -87, -88, -89, -90, -91, -92, -93, -94, 95, -96, -97, -98, -99, -100, -101, 102, -103, -104, -105, 106, -107, -108, -109, -110, 111, 112, -113, -114, -115, -116, 117, -118, 119, -120, -121, -122, -123, -124, -125, -126, -127, -128, -129, -130, -131, -132, -133, -134, -135, -136, 137, -138, -139, -140, -141, 142, -143, -144, -145, 146, -147, 148, -149, -150};

    int mySolution[] = {-1};
    for (auto value : mySolution)
    {
        uint32_t unitClauseValue = 0;
        if (value < 0)
		  {
            unitClauseValue = (((value * -1) << 1) ^ 1);
		  }
        else
		  {
            unitClauseValue = value << 1;
		  }
        assert(unitClauseValue != 0);
        std::cout << "Add " << value << " converted to " << unitClauseValue << " as unit Clause!" << std::endl;
#ifndef NDEBUG
		bool rst = _antom->AddUnit(unitClauseValue); assert(rst);
#else
		_antom->AddUnit(unitClauseValue);
#endif
    }
    return;
}

}
