
/********************************************************************************************
antom.cpp -- Copyright (c) 2014-2016, Sven Reimer

dPermission is hereby granted, free of charge, to any person obtaining a copy of this 
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

// Include antom related headers.
#include "antom.h"
#include "sorter.h"
//#include "counter.h"
#include "control.h"
#include "preprocessor.h"
#include "core.h"
#include "multiplecascade.h"
#include "cascade.h"
#include "timemeasurement.h"
#include "timevariables.h"

namespace antom
{
Antom::Antom(void) :
    AntomBase(),
    _softClauses(),
    _sorterTree(),
    _bestModel(),
    _maxSorterDepth(0),
    _maxsatResult(ANTOM_UNKNOWN),
	_optimum(NULL),
    _lastpartialsorter(false),
    _horizontalbypasses(0),
    _verticalbypasses(0),
    _comparator(0),
    _skipped(0),
    _triggervar(0),
    _lastIndex(0),
    _low(0),
    _high(0),
	_timeVariables(NULL),
	_mainCascade(NULL),
    _mainMultipleCascade(NULL),
    _greatestCommonDivisor(1),
	_formulaStructure(),
    _currentBucketForTare(0),
    _satWeight(0),
    _moreThanTwoWeights(false),
    _topWeight(0),
    _minWeight(1),
    _maxWeight(0),
    _sumOfSoftWeights(0),
    _resultUnknown(false),
	_clausesBefore(0),
    _binaryClausesBefore(0),
    _ternaryClausesBefore(0),
    _addedClauses(0),
    _addedBinaryClauses(0),
    _addedTernaryClauses(0),
    _binaryTopClauses(0),
    _ternaryTopClauses(0),
    _binaryBottomClauses(0),
    _ternaryBottomClauses(0),
    _highestBucketMultiplicator(0),
    _currentBucketMultiplicator(0),
    _estimatedSatWeight(0),
    _diffToSatWeight(0),
    _collectedAssumptions()
{
    _sorterTree.resize(1);
}

Antom::~Antom(void)
{
    for (uint32_t i = 0; i != _sorterTree.size(); ++i)
    {
        for (uint32_t j = 0; j != _sorterTree[i].size(); ++j)
        {
            delete _sorterTree[i][j];
        }
    }
    for (uint32_t i = 0; i != _softClauses.size(); ++i)
    {
        delete _softClauses[i];
    }
}

void Antom::InstanceReset(void)
{
    AntomBase::InstanceReset();
    DataReset();
}

void Antom::DataReset(void)
{
    for (uint32_t i = 0; i != _softClauses.size(); ++i)
    {
        delete _softClauses[i];
    }

    _softClauses.clear();
    for (uint32_t i = 0; i != _sorterTree.size(); ++i)
    {
        for (uint32_t j = 0; j != _sorterTree[i].size(); ++j)
        {
            delete _sorterTree[i][j];
        }
    }
    _sorterTree.clear();
    _sorterTree.resize(1);

    _bestModel.clear();

    _maxSorterDepth = 0;
    _maxsatResult = ANTOM_UNKNOWN;
    _lastpartialsorter = false;

    _horizontalbypasses = 0;
    _verticalbypasses = 0;
    _comparator = 0;
    _skipped = 0;

    _triggervar = 0;
    _lastIndex = 0;
    _low = 0;
    _high = 0;
    _antomSetting->targetOpt = -1;

    _currentBucketForTare = 0;
    _sumOfSoftWeights = 0;
    _satWeight = 0;
    _highestBucketMultiplicator = 0;
    _currentBucketMultiplicator = 0;
    _estimatedSatWeight = 0;
    _diffToSatWeight = 0;
    _antomSetting->decStrat = 0;

    _collectedAssumptions.clear();
}

// Resets SAT and MaxSAT related data
void Antom::Reset(void)
{
    AntomBase::Reset();
    DataReset();

    // TODO: add further status resets
	_antomSetting->ResetAntom();
}

void Antom::SetMaxBounds(uint32_t low, uint32_t high)
{
    assert( low <= high );
    _low = low;
    _high = high;
}

void Antom::SetLowerBound(uint32_t low)
{
    _low = low;
}

void Antom::SetHigherBound(uint32_t high)
{
    _high = high;
}

void Antom::SetOptTarget(int32_t target)
{
    _antomSetting->targetOpt = target;
}

void Antom::SetPreciseTarget(bool val)
{
    _antomSetting->preciseTarget = val;
}

void Antom::SetSearchMode(SearchMode val)
{
    _antomSetting->searchMode = val;
}

void Antom::SetPartialMode(PartialMode val)
{
    _antomSetting->partialMode = val;
}

// Returns the last used index. Does not necessary meet with "variables()"
uint32_t Antom::LastIndex(void) const
{ return _lastIndex; }

// Add a clause to the soft clause database
bool Antom::AddSoftClause(std::vector<uint32_t>& clause, uint64_t weight)
{
    std::vector < uint32_t > sclause( clause.size() );
    uint32_t i = 0;
    for( ; i < clause.size(); ++i )
    {
        sclause[i] = clause[i];
        SetDontTouch(clause[i]>>1);
    }

    if( _antomSetting->incrementalMode > 0 )
    {
        // In incremental mode every softclause is triggered by a global assumption
        if ( _globalPropertyTrigger[_stacksize] == 0 )
        {
            _globalPropertyTrigger[_stacksize] = NewVariable();
        }
        clause.push_back(_globalPropertyTrigger[_stacksize]<<1);
    }

    uint32_t trigger( NewVariable() );
    SetDontTouch(trigger);
    clause.push_back(trigger<<1);

    // Create a new variable which is added to each clause in sorter network
    // Used in incremental mode where the sorter network and all conflict clauses resulting
    // from the network has to be deleted after every step
    SoftClause* sclaus = new SoftClause (trigger<<1, sclause, weight);
    _softClauses.push_back( sclaus );

    _sumOfSoftWeights += weight;

    return AddClause(clause);
}

bool Antom::AddWeightedSoftClause(std::vector<uint32_t>& clause, uint64_t weight)
{
    SoftClause* sclaus = new SoftClause (0, clause, weight);
    _softClauses.push_back( sclaus );

    return true;
}

bool Antom::AddSoftClauseVector()
{
    _sumOfSoftWeights = 0;
    for( uint32_t i = 0; i < _softClauses.size(); ++i )
    {
        std::vector < uint32_t > clause( _softClauses[i]->clause.size() );
        _sumOfSoftWeights += _softClauses[i]->weight;

        assert(_softClauses[i]->clause.size() != 0 );
        for( uint32_t j = 0; j < _softClauses[i]->clause.size(); ++j )
        {
            clause[j] = _softClauses[i]->clause[j];
            SetDontTouch(_softClauses[i]->clause[j]>>1);
        }
        if( _antomSetting->incrementalMode > 0 )
        {
            // In incremental mode every softclause is triggered by a global assumption
            if ( _globalPropertyTrigger[_stacksize] == 0 )
            {
                _globalPropertyTrigger[_stacksize] = NewVariable();
            }
            clause.push_back(_globalPropertyTrigger[_stacksize]<<1);
        }

        uint32_t trigger( NewVariable() );
        SetDontTouch(trigger);
        clause.push_back(trigger<<1);

        // Create a new variable which is added to each clause in sorter network
        // Used in incremental mode where the sorter network and all conflict clauses resulting
        // from the network has to be deleted after every step
        _softClauses[i]->relaxationLit = (trigger<<1);

        if (!AddClause(clause))
        {
            return false;
        }
    }
    return true;
}

// Proceed all soft clauses
// Count positive and negative occurences of all variables in soft clauses
// Set decision strategy for the variables such that the polarity with more occurences is prefered
// If there is no larger occurence value for a polarity, set the decision strategy to "0"
// If "pos" is false the polarity with less occurences is prefered
void Antom::SetStrategyforSoftVariables(bool pos)
{
    std::vector< uint32_t > occur( ((_lastIndex<<1)+2), 0 );

    for( uint32_t i = 0; i != _softClauses.size(); ++i )
    {

        std::vector <uint32_t > lits( _softClauses[i]->clause);
        for ( uint32_t pos = 0; pos != lits.size(); ++pos )
        { ++occur[lits[pos]]; }
    }

    for( uint32_t v = 1; v <= _lastIndex; ++v )
    {
        uint32_t poslit( v<<1 );
        uint32_t neglit( (v<<1)^1 );
        if( occur[poslit] > occur[neglit] )
        {
            if( pos )
            { SetDecisionStrategyForVariable(3,v); }
            else
            { SetDecisionStrategyForVariable(2,v); }
        }
        else if ( occur[poslit] > occur[neglit] )
        {
            if( pos )
            { SetDecisionStrategyForVariable(2,v); }
            else
            { SetDecisionStrategyForVariable(3,v); }
        }
        else
        {
            SetDecisionStrategyForVariable(0,v);
        }
    }
}

void Antom::SetIncrementalMode(uint32_t val)
{ _antomSetting->incrementalMode = val; }

// De-/activates constant detection with SAT
void Antom::SetSatConst(uint32_t val)
{ _antomSetting->satconst = val; }
  
// Check for constants in the candidates list by applying the SAT solver
// Attention! Use with care, process might be costly
// If quickcheck is set, only assumptions are deduced, instead of solving the instance
// -> less powerful, but much faster
bool Antom::FindConstantsWithSat(std::vector< uint32_t >& candidates, bool quickcheck)
{
    assert( _antomSetting->satconst > 0);
    std::cout << __func__ << "currently not fully implemented" << std::endl;
    // Need to distinguish between lastmodel and bestmodel for maxsat
    assert(false);
    exit(0);

    if( _antomSetting->verbosity > 2 )
    { std::cout << "c constant check with SAT..." << std::endl; }

    std::vector< uint32_t > seenvalues(Variables()+1,0);
    std::vector< DecisionStrategy > tempstrategy(Variables()+1,CACHEANDTOGGLE);
    std::vector< bool > checkvariable(Variables()+1,true);

    TrivialAssignment();
    // Fill variables to proceed with current variables
    uint32_t w = 0;
    for( uint32_t i = 0; i != candidates.size(); ++i )
    {
        assert( candidates[i] != 0 );
        assert( candidates[i] < Model().size() );

        if( Model()[candidates[i]] == 0 )
        {
            candidates[w] = candidates[i];
            ++w;
        }
        tempstrategy[candidates[i]] = _core[_sID]->_modeDSForVariables[candidates[i]];
    }
    candidates.resize(w);

    uint32_t satcalls(0);
    std::vector< uint32_t > assumptions;

    bool quickresult(true);
    uint32_t result = Solve(assumptions);
    ++satcalls;

    if( result != ANTOM_SAT )
    {
        // TODO: handle timeout and unsat case
        assert(false);
    }

    _lastModel = Model();

    if( _antomSetting->application == WEIGHTEDMAXSAT )
    {
        uint64_t satWeight = CountSatisfiedSoftClauses(NULL,_lastModel);
        if( satWeight > _satWeight )
        { _satWeight = satWeight; }
    }

    // add initial model
    for( uint32_t i = 0; i != candidates.size(); ++i )
    {
        assert( _lastModel[candidates[i]] != 0 );
        seenvalues[candidates[i]] = _lastModel[candidates[i]];
        // force solver to take opposite polarity
        // 2: always false, 3: always true
        DecisionStrategy strategy( (_lastModel[candidates[i]]&1)?ALWAYSTRUE:ALWAYSFALSE);
        SetDecisionStrategyForVariable(strategy, candidates[i]);
    }

    // Now proceed all candidates

    uint32_t constants(0);
    uint32_t size( candidates.size());
    for( uint32_t i = 0; i != size; ++i )
    {
        uint32_t nextvar( candidates[i] );
        if ( !checkvariable[nextvar] )
        { continue; }
        bool polarity = !(seenvalues[nextvar]&1);

        assert( assumptions.empty());
        // push opposite of last seen value as assumption
        assumptions.push_back( (nextvar << 1)^polarity);

        if( quickcheck )
        {
            quickresult = DeduceAssumptions(assumptions);
        }
        else
        {
            result = Solve(assumptions);
        }
        ++satcalls;

        assumptions.pop_back();

        // "nextvar" is constant
        if ( (result == ANTOM_UNSAT) || !quickresult  )
        {

            TrivialAssignment();
            _lastModel = Model();

            if ( _lastModel[nextvar] == 0 )
            {
                ++constants;
                SetDecisionStrategyForVariable(tempstrategy[nextvar],nextvar);
                if ( !AddUnit( (nextvar << 1)^!polarity ) )
                { return false; }
            }

            for( uint32_t j = i+1; j < size; ++j )
            {
                uint32_t checkvar( candidates[j] );
                // remove already assigned variables from list
                if( _lastModel[checkvar] != 0 )
                {
                    SetDecisionStrategyForVariable(tempstrategy[checkvar],checkvar);
                    candidates[j--] = candidates[--size];
                    candidates.pop_back();
                }
            }
        }
        else if ( (result == ANTOM_SAT) || quickcheck )
        {
            assert( quickresult );
            _lastModel = Model();

            if( _antomSetting->application == WEIGHTEDMAXSAT )
            {
                uint64_t satWeight = CountSatisfiedSoftClauses(NULL,_lastModel);
                if( satWeight > _satWeight )
                { _satWeight = satWeight; }
            }

            // remove candidates which cannot be constant
            for ( uint32_t j = i; j != size; ++j )
            {
                uint32_t checkvar( candidates[j] );
                uint32_t seenvalue( _lastModel[checkvar] );
                uint32_t lastseenvalue( seenvalues[checkvar] );
                assert( quickcheck || (seenvalue != 0 ) );

                if( (seenvalue != 0) && (lastseenvalue != seenvalue) )
                {
                    SetDecisionStrategyForVariable(tempstrategy[checkvar],checkvar);
                    checkvariable[checkvar] = false;
                }
            }
        }
        else
        {
            assert( false );
        }
    }

    // redo decision strategies
    for( uint32_t i = 0; i != candidates.size(); ++i )
    {
        _core[_sID]->SetDecisionStrategyForVariable(tempstrategy[candidates[i]],candidates[i]);
    }

    _core[_sID]->_statistics.constantVariablesBySAT += constants;
    if( _antomSetting->verbosity > 2 )
    {
        std::cout << "c found " << constants << " constants with ";
        if( quickcheck )
        {
            std::cout << "1 SAT-solver call and " << (satcalls-1) << " assumption deductions ";
        }
        else
        {
            std::cout << satcalls << " SAT-solver calls ";
        }
        std::cout << "(checking " << w << " variables)" << std::endl;
    }

    // Do not propagate units...
    return true;
}

uint32_t Antom::MaxSolveLocal(int64_t& optimum)
{
    std::vector< uint32_t > externalAssumptions;
    return MaxSolveLocal( externalAssumptions, optimum  );
}

// Solves the current (partial) MaxSAT formula. Returns SAT/UNSAT and modifies "optimum", representing the minimum number of
// unsatisfied clauses. The modes of operation are:
// mode = 0 --> unsatisfiability-based (multi-threaded: internal portfolio),
// mode = 1 --> satisfiability-based (multi-threaded: internal portfolio).
// mode = 2 --> binary search (multi-threaded: internal portfolio).
uint32_t Antom::MaxSolveLocal(const std::vector< uint32_t >& externalAssumptions, int64_t& optimum)
{
    if( _core[_sID]->EmptyClause() )
    { return ANTOM_UNSAT; }

    struct rusage resources;
    getrusage(RUSAGE_SELF, &resources);
    double timeS  = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
    _control->SetStartTime(timeS);
    _control->SetSumTime(true);

    // Initialize triggervar
    if( _triggervar == 0 )
    {
        _triggervar = NewVariable();
        SetDontTouch(_triggervar);
        AddUnit(_triggervar<<1);
    }

    // There should be no additional assumption in non-incremental mode
    // exception: incomplete and partial mode
    assert( externalAssumptions.empty() || ( _antomSetting->incrementalMode != 0 ) || (_antomSetting->partialMode != NONE)) ;

    if( _antomSetting->incrementalMode > 0 )
    {
        for( uint32_t i = 0; i != externalAssumptions.size(); ++i )
        {
            SetDontTouch(externalAssumptions[i]>>1);
        }
    }

    uint32_t numberofsoftclauses((uint32_t)_softClauses.size());

    _collectedAssumptions = externalAssumptions;

    // No softclauses... nothing to do
    if( numberofsoftclauses == 0 )
    {
        uint32_t res = Solve(_collectedAssumptions);
        _control->SetSumTime(false);
        optimum = 0;
        return res;
    }

    if ( (_antomSetting->partialMode==NONE) && _antomSetting->verbosity > 1 )
    {
        std::cout << "c softclauses: " << numberofsoftclauses << std::endl
                  << "c hardclauses: " << (Clauses()-numberofsoftclauses) << std::endl;
    }

    // Build sorter if not in partial mode
    // Sorters are built seperatly in partialmode
    if( _antomSetting->partialMode==NONE )
    {

        // Do preprocessing?
        if( _antomSetting->doPreprocessing != NOPREPRO )
        {
            // (Incrementally) preprocess formula without sorter
            if ( Preprocess( PREPROCESS ) != ANTOM_UNKNOWN )
            {
                if( _antomSetting->incrementalMode > 0 )
                {
                    InvalidateSoftClauses();
                }
                _control->SetSumTime(false);
                return ANTOM_UNSAT;
            }

            if ( _control->GetTimeOutReached() || _control->GetMemOutReached() )
            {
                ExitTimeout();
                return ANTOM_UNKNOWN;
            }
        }

        // Add the sorting network to the clause database.
        // Collect trigger clauses of soft variables
        assert( _sorterTree[0].empty() );

        Sorter* sorter = new Sorter((uint32_t)_softClauses.size(), this);

        for( uint32_t i = 0; i != _softClauses.size(); ++i )
        {
            sorter->AddSoftClauseToSorter(_softClauses[i]);
        }

        _sorterTree[0].push_back(sorter);

        CreateNextSorter();
    }

    SetDecisionStrategiesForMaxSAT();

    if( _antomSetting->doPreprocessing == INCREMENTAL )
    {
        // (Incrementally) preprocess formula including sorter
        if ( Preprocess( INCREMENTAL ) != ANTOM_UNKNOWN )
        {
            if( _antomSetting->incrementalMode > 0 )
            {
                InvalidateSoftClauses();
            }
            _control->SetSumTime(false);
            return ANTOM_UNSAT;
        }

        if ( _control->GetTimeOutReached() || _control->GetMemOutReached() )
        {
            ExitTimeout();
            return ANTOM_UNKNOWN;
        }
    }

    Sorter* sorter = GetNextSorter();
    uint32_t sortersize = sorter->size();

    // Handle user given bounds
    if( ( _low > 0 ) && _low <= sortersize )
    {
        uint32_t lowerlit((sorter->GetOutputs()[_low-1]<<1)^1);
        bool rst = AddBound(lowerlit);

        if( !rst )
        {
            _control->SetSumTime(false);
            return ANTOM_UNSAT;
        }
    }
    if( ( _high > 0 ) && ( _high <= sortersize ) )
    {
        uint32_t higherlit(sorter->GetOutputs()[sortersize-_high]<<1);

#ifndef NDEBUG
        bool rst = AddBound(higherlit); assert( rst );
#else
        AddBound(higherlit);
#endif
    }

    if( _antomSetting->verbosity > 3 )
    {
        TrivialAssignment();
        sorter->Print();
    }

    // Activate property in incremental mode
    if( ( _antomSetting->incrementalMode > 0 ) )
    {
        for( uint32_t i = 0; i <= _stacksize; ++ i )
        {
            assert( _globalPropertyTrigger[i] != 0 );
            _collectedAssumptions.push_back( (_globalPropertyTrigger[i]<<1)^1 );
        }
    }

    // Initialization.
    uint32_t result(ANTOM_UNKNOWN);
    uint32_t pos(0);

    // What about "mode"?
    switch (_antomSetting->searchMode)
    {

    // UNSAT based
    case UNSATBASED:

        result = UnsatBasedSearch( sorter, pos );
        break;

        // SAT-based
    case SATBASED:

        result = SatBasedSearch( sorter, pos );
        break;

        // Binary search
    case BINARYBASED:

        result = BinaryBasedSearch( sorter, pos );
        break;
    }

    // Add final unit in incomplete mode
    if (_antomSetting->incompleteMode && pos > 0)
    {
#ifndef NDEBUG
        bool rst = AddUnit((sorter->GetOutputs()[pos-1] << 1) ^ 1); assert(rst);
#else
        AddUnit((sorter->GetOutputs()[pos-1] << 1) ^ 1);
#endif		
    }

    // found a valid solution for precise target
    // now calculate value as close as possible to target
    if (_antomSetting->partialMode == NONE && _antomSetting->preciseTarget)
    {
        int32_t gap = _antomSetting->targetOpt - static_cast<int>(sorter->Weight() - pos);
        assert(gap >= 0);
        result = GetPreciseTarget(sorter, pos, gap);
    }

    // _minWeight is one, unless it is called from MaxSolveWeightedMaxSAT
    // with reduced Weight factor of _minWeight
    optimum = (sorter->Weight() - pos) * _minWeight;
    if( _antomSetting->incrementalMode == 1 )
    {
        InvalidateSoftClauses();
    }

    // Remind lower and upper bounds for this sorter
    // Needed for partial modes
    if ( (uint32_t)optimum < _softClauses.size() )
    {
        sorter->SetMinSatisfied(sorter->NumberOfSoftclauses() - pos);
    }

    // Reset collected assumptions for incremental mode
    _collectedAssumptions = externalAssumptions;

    // Return "result".
    _control->SetSumTime(false);

    if (result != ANTOM_UNKNOWN)
    {
        sorter->SetProceeded(true);
    }

    if (_antomSetting->partialMode == NONE)
    {
        _control->SetExtendedResult(_maxsatResult);
    }

    return result;
}

// sat based search for the optimal position, where outputs[pos] is sat and outputs[pos+1] is unsat
uint32_t Antom::SatBasedSearch(Sorter* sorter, uint32_t& pos)
{
    // Initialization.
    uint32_t result = ANTOM_UNSAT;
    uint32_t currentresult = ANTOM_UNKNOWN;

    struct rusage resources;

#ifndef NDEBUG
    uint32_t lastpos = 0;
#endif

    uint32_t maxopt = sorter->Weight();

    // Solve the MaxSAT problem, using a satisfiability-based approach.
    while (true)
    {

        if( _antomSetting->maxInprocess && (_antomSetting->partialMode==NONE) )
        {
            // Do inprocessing
            currentresult = Preprocess(INPROCESS);

            // Have we found the optimum?
            if (currentresult != ANTOM_UNKNOWN)
            { break; }
        }

        // Solve the CNF.
        currentresult = Solve(_collectedAssumptions);

        // Have we found the optimum?
        if (currentresult == ANTOM_UNSAT)
        {
            if (_antomSetting->partialMode == NONE)
            {
                _core[_sID]->SetModel(_lastModel);
                _maxsatResult = ANTOM_SAT;
            }
            break;
        }
        // Reached timeout?
        else if ( currentresult == ANTOM_UNKNOWN )
        {
            if (_antomSetting->partialMode == NONE && _maxsatResult == ANTOM_UNKNOWN_WITH_PRE_RESULT)
            {
                _core[_sID]->SetModel(_lastModel);
            }
            return ANTOM_UNKNOWN;
        }
        else
        {
            result = currentresult;
        }

        _lastModel = Model();

        // Try to satisfy more soft clauses by chance.
        pos = CountSatisfiedSoftClauses( sorter, _lastModel );
        _maxsatResult = ANTOM_UNKNOWN_WITH_PRE_RESULT;

#ifndef NDEBUG
        CheckGates();
#endif

        //std::cout << "lastpos: " << lastpos << " pos: " << pos << std::endl;
#ifndef NDEBUG
        assert( lastpos < pos || pos == 0 );
        lastpos = pos;
#endif

        assert(pos <= maxopt);
        uint32_t localoptimum = maxopt-pos;
        // Output our current best result.
        if ((_antomSetting->verbosity > 0) && (_lastpartialsorter || _antomSetting->partialMode==NONE))
        {
            // _minWeight is one, unless it is called from MaxSolveWeightedMaxSAT
            // with reduced Weight factor of _minWeight
            std::cout << "o " << localoptimum * _minWeight << std::endl;

			// DEBUG
			//if (localoptimum * _minWeight <= 323)
			//			  { exit(0);}
            if( _antomSetting->verbosity>1)
            {
                getrusage(RUSAGE_SELF, &resources);
                double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
                std::cout << "c " << (timeC-_control->GetStartTime()) << "s" <<std::endl;
            }
        }

        // Have we found the very special case with all relaxation variables set to FALSE?
        if (localoptimum == 0)
        { break; }

        // We found our demanded target => return result
        if (_antomSetting->targetOpt >= static_cast<int>(localoptimum))
        {
            if ((_antomSetting->partialMode == NONE) )
            {
                result = ANTOM_UNKNOWN;
                break;
            }

            else
            {
                unsigned int globalopt = UpdateSorterOptimum();
                if ( _antomSetting->targetOpt >= static_cast<int>(_sumOfSoftWeights - globalopt))
                {
                    break;
                }
            }

        }

        uint32_t litbound( (sorter->GetOutputs()[pos] << 1) ^ 1 );

        if ( !AddBound(litbound) )
        {
            _core[_sID]->SetModel(_lastModel);
            break;
        }
    }

    if( _antomSetting->incrementalMode == 1 )
    {
        InvalidateSoftClauses();
    }

    // Return "result".
    return result;
}

// unsat based search for the optimal position, where outputs[pos] is sat and outputs[pos+1] is unsat
uint32_t Antom::UnsatBasedSearch(Sorter* sorter, uint32_t& pos)
{
    // Initialization.
    uint32_t currentresult = ANTOM_UNSAT;
    uint32_t maxopt = sorter->Weight();

    pos = maxopt;

    if( _low > 0 )
    { pos = maxopt-_low-1;  }

    // Now, solve the MaxSAT problem, using a unsatisfiability-based approach.
    while (true)
    {
        // Decrement "pos".
        --pos;

        if( _antomSetting->maxInprocess && (_antomSetting->partialMode==NONE) )
        {
            currentresult = Preprocess(INPROCESS);

            // Have we found the optimum?
            if (currentresult == ANTOM_UNSAT)
            { break; }
        }

        // Output.
        if (_antomSetting->verbosity > 0)
        { std::cout << "c checking o = " << (maxopt - pos - 1) << "..." << std::endl; }

        _collectedAssumptions.push_back( (sorter->GetOutputs()[pos]<<1)^1 );

        // Solve the CNF, taking our assumption into account.
        currentresult = Solve(_collectedAssumptions);

        // Have we found the optimum?
        if (currentresult == ANTOM_SAT)
        {
            // Update "optimum".
            break;
        }
        // Reached timeout?
        else if (currentresult == ANTOM_UNKNOWN)
        {
            ExitTimeout();
            return ANTOM_UNKNOWN;
        }

        if (( _antomSetting->incrementalMode > 0 ) || ( (_antomSetting->partialMode!=NONE) && !_antomSetting->incompleteMode ))
        {
            assert( !_collectedAssumptions.empty() );
            _collectedAssumptions[_collectedAssumptions.size()-1] ^= 1;
        }
        else
        {
            uint32_t assumption( _collectedAssumptions.back() );
            // Remove from collected assumption
            _collectedAssumptions.pop_back();
            // Flip the current assumption and add it as a unit clause to the clause database.
            if (!AddUnit(assumption^1) )
            { assert(currentresult == ANTOM_UNSAT); break; }
        }

        // Did we reached the last element?
        if( pos == 0 )
        {
            assert( currentresult == ANTOM_SAT );
            return currentresult;
        }
    }
    // set pos to the last unsatisfiable value
    ++pos;
    return currentresult;
}

uint32_t Antom::BinaryBasedSearch(Sorter* sorter, uint32_t& pos)
{
    uint32_t maxopt = sorter->Weight();

    uint32_t min = 0;
    uint32_t max = maxopt;

    struct rusage resources;

    // Initialization.
    uint32_t currentresult = Solve(_collectedAssumptions);

    // Have we found the optimum?
    if (currentresult == ANTOM_UNSAT)
    {
        min=max;
        _control->SetSumTime(false);
        return ANTOM_UNSAT;
    }
    // Reached timeout?
    else if ( currentresult == ANTOM_UNKNOWN )
    {
        ExitTimeout();
        return ANTOM_UNKNOWN;
    }

    // a lazy quick and dirty way to preserve the best model
    _lastModel = Model();

    min = CountSatisfiedSoftClauses( sorter, _lastModel );

    uint32_t t = maxopt;
    bool res = true;
    // search for upper bound

    _collectedAssumptions.push_back( (sorter->GetOutputs()[maxopt-1]<<1)^1 );
    do {
        --t;
        _collectedAssumptions.pop_back();
        _collectedAssumptions.push_back( (sorter->GetOutputs()[t]<<1)^1 );
        res = DeduceAssumptions(_collectedAssumptions);
    } while ( !res );

    _collectedAssumptions.pop_back();
    max = t+1;

    pos = (min+max-1) >> 1;

    // Output our current best result.
    if ( _antomSetting->verbosity > 0 )
    {
        std::cout << "o " << (maxopt-min) << std::endl;
        if ( _antomSetting->verbosity > 1 )
        {
            getrusage(RUSAGE_SELF, &resources);
            double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
            std::cout << "c " << (timeC-_control->GetStartTime()) << "s" <<std::endl;
            std::cout << "c known bounds: (" << (maxopt-max) << " / " << (maxopt-min) << ")" << std::endl;
        }
    }


    if ( min > 0 )
    {
        if ( !AddBound( (sorter->GetOutputs()[min-1]<<1)^1 ) )
        {
            max=min;
        }
    }
    if ( max < maxopt )
    {
#ifndef NDEBUG
        bool rst = AddBound(sorter->GetOutputs()[max]<<1);
        assert( rst );
#else
        AddBound(sorter->GetOutputs()[max]<<1);
#endif
    }

    while ( min < max )
    {
        if ( _antomSetting->verbosity > 1 )
        { std::cout << "c TRY with " << (maxopt-pos-1) << std::endl; }

        _collectedAssumptions.push_back( (sorter->GetOutputs()[pos]<<1)^1);

        currentresult = Solve(_collectedAssumptions);

        if ( currentresult == ANTOM_SAT) // SAT
        {
            _lastModel = Model();
            if ( _antomSetting->verbosity > 1 )
            { std::cout << "c Found with " << (maxopt-pos-1) << std::endl; }

            uint32_t nxtmin = CountSatisfiedSoftClauses(sorter, _lastModel);
            assert( nxtmin > min );
            min = nxtmin;

            // Output our current best result.
            if ( _antomSetting->verbosity > 0 )
            {
                std::cout << "o " << (maxopt-min) << std::endl;
                if ( _antomSetting->verbosity > 1)
                {
                    getrusage(RUSAGE_SELF, &resources);
                    double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
                    std::cout << "c " << (timeC-_control->GetStartTime()) << "s" <<std::endl;
                }
            }
        }
        else if ( currentresult == ANTOM_UNSAT ) // UNSAT
        {
            if ( _antomSetting->verbosity > 1 )
            { std::cout << "c Failed with " << (maxopt-pos-1) << std::endl; }

            max = pos;

            assert( !_collectedAssumptions.empty() );
            _collectedAssumptions[_collectedAssumptions.size()-1] ^= 1;

        }
        else if ( currentresult == ANTOM_UNKNOWN ) // Reached timeout?
        {
            ExitTimeout();
            if( _antomSetting->incrementalMode > 0 )
            {
                for( uint32_t i = 0; i != _collectedAssumptions.size(); ++i )
                { SetDontTouch((_collectedAssumptions[i]>>1),false); }
            }
            return ANTOM_UNKNOWN;
        }
        else
        { assert(false); }

        if ( ( _antomSetting->incrementalMode == 0 ) && ( (_antomSetting->partialMode== NONE) || _antomSetting->incompleteMode ) )
        {
#ifndef NDEBUG
            bool rst(false); rst = AddUnit( _collectedAssumptions.back() ); assert( rst );
#else
            AddUnit( _collectedAssumptions.back() );
#endif
            _collectedAssumptions.pop_back();
        }

        // Do inprocessing
        if( _antomSetting->maxInprocess && (_antomSetting->partialMode==0)  )
        {
            currentresult = Preprocess(INPROCESS);

            // instance unsatisfiable with added assumption?
            // then the currentoptimum is the best result
            if( currentresult != ANTOM_UNKNOWN )
            { break; }
        }

        pos = (min+max-1) >> 1;

        if ( _antomSetting->verbosity > 1 )
        { std::cout << "c known bounds: (" << (maxopt-max) << " / " << (maxopt-min) << ")" << std::endl; }

    }

    if ( maxopt>pos )
    {
        std::vector< uint32_t > model_assumption;
        model_assumption.reserve(_lastIndex+1);
        // solve with the best found model, so that the antom-model is set correctly
        // just add the model for the original clauses.
        // Due to our "by chance" method the virtual best model is not complete
        for ( uint32_t m = 1; m <= _lastIndex; ++m )
        {
            if( _lastModel[m] != 0 )
            { model_assumption.push_back(_lastModel[m]);}
        }

        // this is the quick and dirty way! TODO
        currentresult = Solve(model_assumption);

        assert( currentresult == ANTOM_SAT );
    }
    ++pos;
    return currentresult;
}

uint32_t Antom::GetPreciseTarget(Sorter* sorter, uint32_t& pos, int32_t gap)
{
    // TODO: not really working yet
    // Check how to avoid by chance soft clauses
    assert(false);
    exit(0);
    if (gap==0)
    {
        return ANTOM_SAT;
    }

    // forbid last optimum
    uint32_t currentResult(ANTOM_UNKNOWN);
    int32_t lastGap(gap);
    do
    {
        if (pos == 0)
        {
            break;
        }
        uint32_t lit = sorter->GetOutputs()[pos-1]<<1;
        std::cout << "add unit " << helper::Lit(lit) << std::endl;
#ifndef NDEBUG
        bool rst = AddUnit(lit); assert(rst);
#else
        AddUnit(lit);
#endif
        lastGap = gap;
        currentResult = Solve();

        std::cout << "currentresult: " << currentResult << std::endl;
        // Timeout Result
        if (currentResult == ANTOM_UNKNOWN)
        {
            _core[_sID]->SetModel(_lastModel);
            _maxsatResult = ANTOM_UNKNOWN_WITH_PRE_RESULT;
            return ANTOM_UNKNOWN;
        }
        // There may be no solution less than found optimum -> we are done
        else if(currentResult == ANTOM_UNSAT)
        {
            break;
        }
        _lastModel = Model();
        pos = CountSatisfiedSoftClauses(sorter, _lastModel);
        std::cout << "pos: " << pos << std::endl;
        gap = _antomSetting->targetOpt - static_cast<int>(sorter->Weight() - pos);
        std::cout << "lastgap: " << lastGap << " newgap: " << gap << std::endl;
    }
    while (gap > 0);

    // readjust model if the lastGap is closer than the last result
    if (lastGap < -gap)
    {
        _core[_sID]->SetModel(_lastModel);
    }

    _maxsatResult = ANTOM_SAT;

    return ANTOM_SAT;
}

// Returns false if added bounded lead to unsat
bool Antom::AddBound(uint32_t lit)
{
    if( ( _antomSetting->incrementalMode > 0 ) || ( !_lastpartialsorter && (_antomSetting->partialMode!=NONE) ) )
    {
        //std::cout << __func__ << " assumption " << helper::Lit(lit) << std::endl;
        _collectedAssumptions.push_back( lit );
    }
    else
    {
        //std::cout << __func__ << " unit " << helper::Lit(lit) << std::endl;
        assert(!_antomSetting->incompleteMode || _lastpartialsorter);
        // if the resulting problem is unsatisfied, we have found our optimum!
        if (!AddUnit( lit ) )
        {
            return false;
        }
    }
    return true;
}


uint32_t Antom::MaxSolve(int64_t& optimum)
{
    std::vector< uint32_t > externalAssumptions;
    return MaxSolve( externalAssumptions, optimum );
}

uint32_t Antom::MaxSolve(const std::vector< uint32_t >& externalAssumptions, int64_t& optimum)
{
  _antomSetting->application = MAXSAT;
    assert(!_antomSetting->preciseTarget || _antomSetting->targetOpt >= 0);
    assert(!_antomSetting->preciseTarget || _antomSetting->searchMode == SATBASED);

    if( _antomSetting->partialMode == NONE)
    {
        return MaxSolveLocal(externalAssumptions, optimum);
    }

    // Consistency checks.
    assert(!_antomSetting->incompleteMode || ( _antomSetting->searchMode == SATBASED ));

    if( _core[_sID]->EmptyClause() )
    {
        if( _antomSetting->incrementalMode > 0 )
        {
            InvalidateSoftClauses();
        }
        return ANTOM_UNSAT;
    }

    struct rusage resources;
    getrusage(RUSAGE_SELF, &resources);
    double timeS  = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
    _control->SetStartTime(timeS);
    _control->SetSumTime(true);

    // Set incrementalmode to "incremental + partial"
    if( _antomSetting->incrementalMode == 1 )
    { _antomSetting->incrementalMode = 2; }

    for( uint32_t i = 0; i != externalAssumptions.size(); ++i )
    { SetDontTouch(externalAssumptions[i]>>1); }

    uint32_t numberofsoftclauses((uint32_t)_softClauses.size());

    if( numberofsoftclauses == 0 )
    {
        uint32_t res = Solve(externalAssumptions);
        _control->SetSumTime(false);
        if( _antomSetting->incrementalMode > 0 )
        {
            InvalidateSoftClauses();
        }
        return res;
    }

    // Do preprocessing?
    if( _antomSetting->doPreprocessing != NOPREPRO )
    {
        // (Incrementally) preprocess formula without sorter
        if ( Preprocess( PREPROCESS ) != ANTOM_UNKNOWN )
        {
            if( _antomSetting->incrementalMode > 0 )
            {
                InvalidateSoftClauses();
            }
            _control->SetSumTime(false);
            return ANTOM_UNSAT;
        }

        if ( _control->GetTimeOutReached() || _control->GetMemOutReached() )
        {
            ExitTimeout();
            return ANTOM_UNKNOWN;
        }
    }

    uint32_t currentresult(ANTOM_UNKNOWN);
    uint32_t result(ANTOM_UNKNOWN);
    uint32_t depth(0);
    int32_t currentoptimum(_sumOfSoftWeights);
    optimum = currentoptimum+1;

    if( _triggervar == 0 )
    {
        _triggervar = NewVariable();
        SetDontTouch(_triggervar);
        AddUnit(_triggervar<<1);
    }

    if ( _antomSetting->verbosity > 1 )
    {
        std::cout << "c softclauses: " << numberofsoftclauses << std::endl
                  << "c hardclauses: " << (Clauses()-numberofsoftclauses) << std::endl;
    }

    // calc and set all trivially conflicting softclauses
    CheckAllConflictingSoftclauses();
    // Constant check
    CalcSplitWidth();

    uint32_t prepart(0);
    // build datastructure for parts
    while( prepart*_antomSetting->splittedWidth < numberofsoftclauses )
    {
        uint32_t s ( prepart*_antomSetting->splittedWidth );
        uint32_t maxsize = s+_antomSetting->splittedWidth;

        if( numberofsoftclauses < maxsize )
        { maxsize = numberofsoftclauses; }

        Sorter* sorter = new Sorter(_antomSetting->splittedWidth, this);

        for( ; s < maxsize; ++s )
        {
            sorter->AddSoftClauseToSorter(_softClauses[s]);
        }
        _sorterTree[0].push_back(sorter);
        ++prepart;
    }

    if( _lastIndex == 0 )
    {
        _lastIndex = Variables();
    }

    // Add all original variables as candidates
    std::vector< uint32_t > satconstcandidates(_lastIndex,0);
    for( uint32_t v = 1; v <= _lastIndex; ++v )
    { satconstcandidates[v-1] = v; }

    if ( _antomSetting->verbosity > 1 )
    { std::cout << "c start solving depth " << depth << " width: " << _antomSetting->splittedWidth << std::endl; }

    Sorter* currentsorter = NULL;
    do
    {
        currentsorter = MergeNextSorter();
        currentsorter->SetProceedNext(true);

        // Reached last sorter -> we can treat this "unpartially"
        if ( _sorterTree[_maxSorterDepth][0]->NumberOfSoftclauses() ==  numberofsoftclauses )
        {
            _lastpartialsorter = true;
        }

        std::vector< uint32_t > localassumptions = externalAssumptions;
        CollectRelaxationLits( localassumptions );

        // Add lower bounds of satisfied soft clauses as assumption for this run
        uint32_t minsat = SetLowerPartialParts(localassumptions);

        if ( minsat > (_softClauses.size()-currentsorter->NumberOfSoftclauses()) )
        {
            uint32_t pos = minsat-(_softClauses.size()-currentsorter->NumberOfSoftclauses());

            uint32_t lit((currentsorter->GetOutputs()[pos-1]<<1)^1);
#ifndef NDEBUG
            bool rst = AddUnit(lit); assert(rst);
#else
            AddUnit(lit);
#endif			
        }

        // solve (merged) instance
        int64_t partialoptimum = -1;

        if( _antomSetting->verbosity > 1 )
        {
            uint32_t sum(0);
            for( uint32_t i = 0; i != _sorterTree.size(); ++i )
            {
                std::cout << "c " << i << ": ";
                for( uint32_t j = 0; j != _sorterTree[i].size(); ++j )
                {
                    std::cout << _sorterTree[i][j]->NumberOfSoftclauses() << ", ";
                    sum += (uint32_t)_sorterTree[i][j]->NumberOfSoftclauses();
                }
                std::cout << " (" << sum << ")" << std::endl;
            }
            assert( sum <= _softClauses.size() );
        }


        if( ( _antomSetting->satconst > 0 ) && (currentsorter->NumberOfSoftclauses() < numberofsoftclauses ) )
        {
#ifndef NDEBUG
            bool rst = FindConstantsWithSat(satconstcandidates, _antomSetting->satconst==1); assert(rst);
#else
            FindConstantsWithSat(satconstcandidates, true);
#endif
        }

        if( _antomSetting->maxInprocess )
        { currentresult = Preprocess(INPROCESS); }
        else
        { currentresult = ANTOM_UNKNOWN; }

        currentresult = MaxSolveLocal( localassumptions, partialoptimum );
        _control->SetSumTime(true);

        if( currentresult == ANTOM_SAT )
        {
            result = ANTOM_SAT;
            _maxsatResult = ANTOM_UNKNOWN_WITH_PRE_RESULT;

            currentoptimum = _sumOfSoftWeights - UpdateSorterOptimum();
            assert( currentoptimum >= 0 );

            if( currentoptimum < optimum )
            {
                _bestModel = _lastModel;
                optimum = currentoptimum;
                if ( _antomSetting->verbosity > 0 )
                {
                    std::cout << "o " << optimum << std::endl;
                    if( _antomSetting->verbosity > 1 )
                    {
                        getrusage(RUSAGE_SELF, &resources);
                        double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
                        std::cout << "c " << (timeC-timeS) << "s" << std::endl;
                    }
                }
            }

            if (_antomSetting->verbosity > 1)
            { std::cout << "c partial optimum: " << currentoptimum << std::endl; }
        }
        else if ( currentresult == ANTOM_UNSAT )
        {
            _maxsatResult = ANTOM_UNSAT;
            _control->SetExtendedResult(ANTOM_UNSAT);
            _control->SetSumTime(false);
            if ( _antomSetting->incrementalMode > 0 )
            {
                InvalidateSoftClauses();
            }
            return ANTOM_UNSAT;
        }
        // Timeout
        else if ( currentresult == ANTOM_UNKNOWN )
        {
            if (_maxsatResult == ANTOM_UNKNOWN_WITH_PRE_RESULT)
            {
                // bestmodel might not have the correct size,
                // in this case lastmodel is also the best model
                if(_bestModel.size() < _lastModel.size() )
                {
                    _core[_sID]->SetModel(_lastModel);
                }
                else
                {
                    _core[_sID]->SetModel(_bestModel);
                }
            }

            _control->SetExtendedResult(_maxsatResult);
            ExitTimeout();
            return ANTOM_UNKNOWN;
        }

        if ( _antomSetting->incompleteMode && _sorterTree[0].empty() )
        { break; }

        // We have found our demanded optimum
        if (_antomSetting->targetOpt >= static_cast<int>(optimum))
        {
            assert(_maxsatResult == ANTOM_UNKNOWN_WITH_PRE_RESULT);
            _core[_sID]->SetModel(_bestModel);
            _control->SetSumTime(false);
            if( _antomSetting->incrementalMode > 0 )
            { InvalidateSoftClauses(); }
            _control->SetExtendedResult(_maxsatResult);
            return ANTOM_UNKNOWN;
        }

        currentsorter->SetProceedNext(false);
        assert( !_sorterTree[_maxSorterDepth].empty() );
        // Do until one part is left and all parts are considered
    } while ( currentsorter->NumberOfSoftclauses() < numberofsoftclauses );

    if ( _antomSetting->verbosity > 1 )
    {
        std::cout << "c comparator: " << _comparator << " skipped: " << _skipped << std::endl;
    }

    _core[_sID]->SetModel(_bestModel);
    _control->SetSumTime(false);
    if( _antomSetting->incrementalMode > 0 )
    { InvalidateSoftClauses(); }

    _maxsatResult = ANTOM_SAT;
    _control->SetExtendedResult(_maxsatResult);
    return result;
}

uint32_t Antom::SetLowerPartialParts(std::vector< uint32_t >& localassumptions) const
{
    uint32_t sum(0);
    // proceed current considered trigger parts
    for( uint32_t i = 1; i != _sorterTree.size(); ++i )
    {
        for( uint32_t j = 0; j != _sorterTree[i].size(); ++j )
        {
            uint32_t minsat( _sorterTree[i][j]->GetMinSatisfied() );
            if( minsat > 0 )
            {
                uint32_t lit( (_sorterTree[i][j]->GetOutputs()[minsat-1]<<1)^1 );
                localassumptions.push_back(lit);
                sum += minsat;
            }
        }
    }
    assert( sum <= _softClauses.size() );
    return sum;
}

int32_t Antom::UpdateSorterOptimum(void)
{
    int32_t result(0);
    // proceed current considered trigger parts
    for( uint32_t i = 0; i != _sorterTree.size(); ++i )
    {
        for( uint32_t j = 0; j != _sorterTree[i].size(); ++j )
        {
            uint32_t pos = CountSatisfiedSoftClauses( _sorterTree[i][j], _lastModel );
            result += pos;
        }
    }
    return result;
}

// If sorter == NULL, we want to count for all soft clauses
uint64_t Antom::CountSatisfiedSoftClauses(Sorter* sorter, const std::vector<uint32_t>& model)
{

    std::vector< SoftClause* > softclauses;
    if( sorter == NULL )
    {
        softclauses = _softClauses;
    }
    else
    {
        softclauses = sorter->GetSoftClauses();
    }

    uint64_t result = CountSatisfiedSoftClauses(softclauses, model);

    if( sorter != NULL )
    {
        sorter->SetMinSatisfied(result);
    }
    return result;
}


uint64_t Antom::CountSatisfiedSoftClauses(std::vector< SoftClause* > softclauses, const std::vector<uint32_t>& model)
{

//std::cout << __func__ << std::endl;
//    if (_antomSetting->verbosity > 4)
//    {
//        std::cout << __PRETTY_FUNCTION__ << std::endl;
//        std::cout << "SC.size: " << softclauses.size() << std::endl;
//        std::cout << "Clauses.size: " << Clauses() << std::endl;
//    }

uint64_t result(0);
// Proceed all soft clauses
for( uint32_t i = 0; i != softclauses.size(); ++i )
{
    uint32_t relaxlit = softclauses[i]->relaxationLit;

	//std::cout << helper::Lit(relaxlit) << " " << helper::Lit(model[relaxlit>>1]) << "; " ;
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
                result += softclauses[i]->weight;
                softclauses[i]->lastassignment = relaxlit^1;
                //std::cout << "++" << result;
                break;
            }
        }
        if ( pos == clause.size() )
        {
            softclauses[i]->lastassignment = relaxlit;
        }
    }
    else if ( model[relaxlit>>1] != 0 )
    {
        assert( model[relaxlit>>1] == (relaxlit^1));
        softclauses[i]->lastassignment = relaxlit^1;
        result += softclauses[i]->weight;
        //std::cout << " ++relax sat " << result;
    }
    //std::cout << std::endl;
}
return result;

}

void Antom::CalcSplitWidth(void)
{
    uint32_t parts(1);
    uint32_t splitwidth((uint32_t)_softClauses.size());

    if( _antomSetting->maxWidth != 0 )
    {
        // First get largest power of 2 number below "_maxWidth"
        uint32_t factor( 1 );
        while ( ( factor * _antomSetting->maxWidth ) < splitwidth )
        { factor <<= 1; }

        // If there is any left over split into "(largest power of 2 number) - 1 + rest" parts
        _antomSetting->splittedWidth = splitwidth/factor;

        if( (splitwidth%factor) != 0 )
        { ++_antomSetting->splittedWidth; }
    }

    if ( _antomSetting->maxParts != 0 )
    {
        while ( parts > _antomSetting->maxParts )
        { splitwidth = splitwidth<<1; parts = parts>>1; }
        _antomSetting->splittedWidth = splitwidth;
    }
}

void Antom::SetSplittedWidth(uint32_t width)
{ assert( width > 0 ); _antomSetting->splittedWidth = width; }

void Antom::SetMaxWidth(uint32_t width)
{ _antomSetting->maxWidth = width; }

void Antom::SetMaxParts(uint32_t parts)
{ _antomSetting->maxParts = parts; }

void Antom::SetSortSoftClauses(uint32_t val)
{ assert( val < 4 ); _antomSetting->sortSoftClauses = val; }

void Antom::SetSkip(bool val)
{ _antomSetting->doSkipping = val; }
  
void Antom::SetCSC(uint32_t val)
{ _antomSetting->setCSC = val; }
  
void Antom::SetRelaxationLits(bool val)
{ _antomSetting->setRelaxationLits = val; }
  
void Antom::SetIncompleteMode(bool val)
{
    _antomSetting->incompleteMode = val;
    if ( val )
    { _antomSetting->partialMode = BREADTHFIRST; }
}

void Antom::SetGridMode(uint32_t val)
{ assert( val < 4 ); _antomSetting->bypassGrid = val; }
void Antom::SetHorizontalWidth(uint32_t val)
{ assert( val > 0 ); _antomSetting->horizontalWidth = val; }

void Antom::SetNetworktype(SorterType val)
{ _antomSetting->networkType = val; }

void Antom::SetMaxInprocessing(bool val)
{ _antomSetting->maxInprocess = val; }

void Antom::CheckAllConflictingSoftclauses(void)
{
    if( _antomSetting->setCSC != 2 )
    { return; }

    size_t size = _softClauses.size();

    if( _antomSetting->verbosity > 1 )
    { std::cout << "c " << __func__ << " " << size << std::endl; }

    uint32_t bypasses(0);
    uint32_t unsats(0);

    for( size_t i = 0; i != size; ++i )
    {
        std::vector< uint32_t > assumptions;
        uint32_t lit( _softClauses[i]->relaxationLit );
        assumptions.push_back( lit^1 );

        // Deduce activation of current soft clause
        bool rst = DeduceAssumptions(assumptions);
        _softClauses[i]->contra = 0;

        // Instance is not solveable if current soft clause is activated
        if( !rst )
        {
            rst = AddUnit( lit ); assert( rst );
            ++unsats;
        }
        else
        {
            const std::vector< uint32_t >& solvermodel = Model();

            // Check model
            for( uint32_t j = 0; j < size; ++j )
            {
                uint32_t checklit( _softClauses[j]->relaxationLit );
                assert( (checklit>>1) < solvermodel.size() );
                // Contradiction... j'th soft clause is deactivated if i'th soft clause is activated
                if( solvermodel[checklit>>1] == checklit )
                {
                    ++bypasses;
                    ++_softClauses[i]->contra;
                    std::vector < uint32_t > clause( 2, 0 );
                    clause[0] = lit;
                    clause[1] = checklit;
                    rst = AddClause(clause); assert(rst);
                }
                // The i'th and j'th soft clause can be activated at once
                else if( ( i != j ) && ( solvermodel[checklit>>1] == (checklit^1) ) )
                {
                    ++bypasses;
                    std::vector < uint32_t > clause( 2, 0 );
                    clause[0] = lit;
                    clause[1] = (checklit)^1;
                    rst = AddClause(clause); assert(rst);
                }
            }
        }
    }

    if( _antomSetting->verbosity > 1 )
    {
        std::cout << "c find overall " << bypasses << " bypasses and " << unsats << " unsats" << std::endl; }

    if( _antomSetting->sortSoftClauses == 1 )
    { std::sort( _softClauses.begin(), _softClauses.end(), SoftClauseContraSorter() ); }
    else if( _antomSetting->sortSoftClauses == 2 )
    {	std::sort( _softClauses.rbegin(), _softClauses.rend(), SoftClauseContraSorter() ); }
    else if ( _antomSetting->sortSoftClauses == 3 )
    {	std::random_shuffle(_softClauses.begin(),_softClauses.end());  }
}

void Antom::CheckAllWeightedConflictingSoftclauses(void)
{
    if( _antomSetting->setCSC != 2 )
    { return; }

    std::size_t size = _softClauses.size();

    if( _antomSetting->verbosity > 1 )
    { std::cout << "c " << __func__ << " " << size << std::endl; }

    uint32_t bypasses(0);
    uint32_t unsats(0);

    for( std::size_t i = 0; i != size; ++i )
    {
        std::vector< uint32_t > assumptions;
        uint32_t lit( _softClauses[i]->relaxationLit );
        assumptions.push_back( lit^1 );

        // Deduce activation of current soft clause
        bool rst = DeduceAssumptions(assumptions);
        _softClauses[i]->contra = 0;

        // Instance is not solveable if current soft clause is activated
        if( !rst )
        {
            rst = AddUnit( lit ); assert( rst );
            ++unsats;
        }
        else
        {
            const std::vector< uint32_t >& solvermodel = Model();

            // Check model
            for( uint32_t j = 0; j < size; ++j )
            {
                uint32_t checklit( _softClauses[j]->relaxationLit );
                assert( (checklit>>1) < solvermodel.size() );
                // Contradiction... j'th soft clause is deactivated if i'th soft clause is activated
                if( solvermodel[checklit>>1] == checklit )
                {
                    ++bypasses;
                    ++_softClauses[i]->contra;
                    std::vector < uint32_t > clause( 2, 0 );
                    clause[0] = lit;
                    clause[1] = checklit;
                    rst = AddClause(clause); assert(rst);
                }
                // The i'th and j'th soft clause can be activated at once
                else if( ( i != j ) && ( solvermodel[checklit>>1] == (checklit^1) ) )
                {
                    ++bypasses;
                    std::vector < uint32_t > clause( 2, 0 );
                    clause[0] = lit;
                    clause[1] = (checklit)^1;
                    rst = AddClause(clause); assert(rst);
                }
            }
        }
    }

    if( _antomSetting->verbosity > 1 )
    {
        std::cout << "c find overall " << bypasses << " bypasses and " << unsats << " unsats" << std::endl; }

    if( _antomSetting->sortSoftClauses == 1 )
    { std::sort( _softClauses.begin(), _softClauses.end(), SoftClauseContraSorter() ); }
    else if( _antomSetting->sortSoftClauses == 2 )
    {	std::sort( _softClauses.rbegin(), _softClauses.rend(), SoftClauseContraSorter() ); }
    else if ( _antomSetting->sortSoftClauses == 3 )
    {	std::random_shuffle(_softClauses.begin(),_softClauses.end());  }
}



// Adds last assignment of soft triggers to "assumptions"
void Antom::CollectRelaxationLits(std::vector< uint32_t >& assumptions)
{
    if( !_antomSetting->setRelaxationLits )
    { return; }

    // Add last assignment of soft triggers to assumptions

    // proceed proceeded parts
    for( uint32_t i = 0; i < _sorterTree.size(); ++i )
    {
        for(uint32_t j = 0; j != _sorterTree[i].size(); ++j )
        {
            Sorter* proceedSorter = _sorterTree[i][j];
            // Skip current sorter
            if( proceedSorter->GetProceedNext())
            { continue; }

            for( uint32_t j = 0; j != proceedSorter->NumberOfSoftclauses(); ++j )
            {
                SoftClause* sclause = proceedSorter->GetSoftClauses()[j];
                if ( sclause->lastassignment != 0 )
                {
                    assumptions.push_back( sclause->lastassignment );
                }
                else
                {
                    assumptions.push_back( sclause->relaxationLit );
                }
            }
        }
    }
}


// Deletes the soft clause datastructure and all clauses containing these clauses
void Antom::InvalidateSoftClauses(void)
{
    assert( _antomSetting->incrementalMode > 0 );

    // Set triggers to true, i.e. invalid softclauses
    for ( uint32_t i = 0; i != _softClauses.size(); ++i )
    {
        assert( _lastIndex < (_softClauses[i]->relaxationLit>>1) );

        // Remove "don't touch" status for soft clause variables
        for( uint32_t pos = 0; pos != _softClauses[i]->clause.size(); ++pos )
        { SetDontTouch(_softClauses[i]->clause[pos]>>1, false); }
    }

    // clear soft clause data structure
    _softClauses.clear();
    _sorterTree.clear();

    assert( _lastIndex < _globalPropertyTrigger[_stacksize] );
    // reset trigger
    _globalPropertyTrigger[_stacksize] = 0;

    // Delete variables and clauses from last time step
    ClearVariables(_lastIndex+1, _core[_sID]->Variables());

    // Reset max variable index to "_lastIndex", the datastructure capacity is conserved
    SetMaxIndex(_lastIndex);

    for (uint32_t t = 0; t < _antomSetting->threads; ++t)
    {
        _preprocessor[t]->SetPreVarIndex(_lastIndex+1);
        _preprocessor[t]->SetPreClauseIndex();
    }

    // Checks, whether all clauses with variable index > "_lastIndex" are deleted
    assert( _core[_sID]->CheckMaxIndex(_lastIndex) );
}

// creates sorter for depth
Sorter* Antom::CreateNextSorter(uint32_t depth, uint32_t targetdepth)
{
    assert( !_sorterTree[depth].empty() );
    // clear model
    TrivialAssignment();

    Sorter* sorter = _sorterTree[depth].back();
    _sorterTree[depth].pop_back();

    sorter->EncodeSorter(targetdepth);

    if( _maxSorterDepth < targetdepth )
    {
        _maxSorterDepth = targetdepth;
        _sorterTree.resize(targetdepth+1);
    }
    _sorterTree[targetdepth].push_back( sorter );

    if( _antomSetting->verbosity > 2 )
    {
        std::cout << "c comparator: " << _comparator << " skipped: " << _skipped << std::endl;
        std::cout << "c horizontal bypasses: " << _horizontalbypasses << " vertical bypasses: " << _verticalbypasses << std::endl;
    }
    return _sorterTree[targetdepth].back();
}

Sorter* Antom::MergeNextSorter(void)
{
    // Go throw sorters, starting with highest depth
    // If two sorters of the same depth are available -> merge them
    assert( _sorterTree.size() > _maxSorterDepth );
    uint32_t depth = _maxSorterDepth;
    uint32_t depth2 = 0;
    switch( _antomSetting->partialMode )
    {
    case DEPTHFIRST :
        for( ; depth > 0; --depth )
        {
            if( _sorterTree[depth].size() > 1 )
            {
                break;
            }
        }

        break;
    case BREADTHFIRST :
        depth = 0;
        for( ; depth <= _maxSorterDepth; ++depth )
        {
            if( _sorterTree[depth].size() > 0 )
            {
                break;
            }
        }

        break;
    default:
        assert(false);
        break;
    }

    // All sorters are built, no merge on higher levels possible
    // -> We have to merge over different levels
    if( (depth == 0 && _sorterTree[0].empty()) || ( (depth > 0) && _sorterTree[depth].size() == 1 ) )
    {
        for( ; depth <= _maxSorterDepth; ++depth )
        {
            if( _sorterTree[depth].size() > 0 )
            {
                if( depth2 != 0 )
                { break; }
                depth2 = depth;
            }
        }
    }
    else
    {
        depth2 = depth;
    }

    // Reched level 0
    // -> We have to create the next sorter for level 1
    if ( depth == 0 )
    {
        return CreateNextSorter();
    }

    Sorter* tempSorter1 = _sorterTree[depth].back();
    _sorterTree[depth].pop_back();
    Sorter* tempSorter2 = _sorterTree[depth2].back();
    _sorterTree[depth2].pop_back();

    tempSorter2->MergeWithSorter(*tempSorter1);

    delete tempSorter1;

    // Pushed merged list on next level
    if( depth == _maxSorterDepth )
    {
        ++_maxSorterDepth;
        _sorterTree.resize(_maxSorterDepth+1);
    }
    _sorterTree[depth+1].push_back(tempSorter2);
    return _sorterTree[depth+1].back();
}

Sorter* Antom::GetNextSorter(void)
{
    switch( _antomSetting->partialMode )
    {
    case NONE :
        assert( _sorterTree[1].size() == 1 );
        assert( _maxSorterDepth == 1 );
        return _sorterTree[1].back();

        break;
    case DEPTHFIRST :
        for (uint32_t i = 1; i != _sorterTree.size(); ++i )
        {
            if( !_sorterTree[i].empty() )
            {
                assert( !_sorterTree[i].back()->GetProceeded() );
                return _sorterTree[i].back();
            }
        }
        assert(false);
        break;
    case BREADTHFIRST :
        for (uint32_t i = _maxSorterDepth; i > 0; --i )
        {
            if( !_sorterTree[i].empty() )
            {
                assert( !_sorterTree[i].back()->GetProceeded() );
                return _sorterTree[i].back();
            }
        }
        assert(false);
        break;
    default :
        std::cout << "c unknown partial mode " << _antomSetting->partialMode << " => exit" << std::endl;
        exit(0);
        break;
    }
    return NULL;
}

void Antom::ExitTimeout(void)
{
    assert( (_antomSetting->cpuLimit > 0.0 ) || ( _antomSetting->memLimit > 0 ) || _antomSetting->targetOpt >= 0);

    if( _antomSetting->incrementalMode > 0 )
    {
        InvalidateSoftClauses();
    }
    _control->SetSumTime(false);
}

void Antom::SetDecisionStrategiesForMaxSAT(void)
{
    if (_antomSetting->decStrat <= 1)
    {
        // set outputs of last sorter to first group
        assert( !_sorterTree.empty());
        if( !_sorterTree[_maxSorterDepth-1].empty() )
        {
            assert( _sorterTree[_maxSorterDepth-1][0] != NULL );
            for (uint32_t i = 0; i < _sorterTree[_maxSorterDepth-1][0]->size(); ++i)
            {
                uint32_t outputvar( _sorterTree[_maxSorterDepth-1][0]->GetOutputs()[i] );
                SetDecisionStrategyForVariable(2,outputvar);
            }
        }

        if (_antomSetting->decStrat == 0)
        {
            // set all relaxationliterals to second group
            for (uint32_t i = 0; i < _softClauses.size(); ++i)
            {
                uint32_t relaxvar(_softClauses[i]->relaxationLit >> 1);
                SetDecisionStrategyForVariable(2,relaxvar);
            }

            // set tare variables to third group
            for (int32_t j = _maxSorterDepth - 2; j >= 0; --j)
            {
                assert(_sorterTree[j][0] != NULL );
                if ( _sorterTree[j].empty() )
                { continue; }

                for( uint32_t k = 0; k != _sorterTree[j][0]->GetTares().size(); ++k )
                {
                    uint32_t tarevar(_sorterTree[j][0]->GetTares()[k]);
                    SetDecisionStrategyForVariable(3,tarevar);
                }
            }
        }
    }
}

  void Antom::PrintStatistics() const
  {
	// Output.
	std::cout << "c #ID fastest thread.....: " << SolvingThread()            << std::endl
			  << "c #variables.............: " << Variables()                << std::endl
			  << "c    #used...............: " << UsedVariables()            << std::endl
			  << "c #clauses...............: " << Clauses()                  << std::endl
			  << "c    #binary.............: " << CurrentBinaryClauses()     << std::endl
			  << "c    #ternary............: " << CurrentTernaryClauses()    << std::endl
			  << "c    #nary...............: " << CurrentNaryClauses()       << std::endl
			  << "c #literals..............: " << Literals()                 << std::endl
			  << "c #decisions.............: " << Decisions()                << std::endl
			  << "c #bcp operations........: " << Bcps()                     << std::endl
			  << "c #implications..........: " << Implications()             << std::endl
			  << "c #conflicts.............: " << Conflicts()                << std::endl
			  << "c #restarts..............: " << Restarts()                 << std::endl
			  << "c #blocked restarts......: " << BlockedRestarts()          << std::endl
			  << "c #simplifications.......: " << Simplifications()          << std::endl
			  << "c #synchronizations......: " << Synchronizations()         << std::endl
			  << "c #lhbr clauses..........: " << Lhbr()                     << std::endl
			  << "c #learnt unit clauses...: " << LearntUnitClauses()        << std::endl
			  << "c #learnt binary clauses.: " << LearntBinaryClauses()      << std::endl
			  << "c #learnt ternary clauses: " << LearntTernaryClauses()     << std::endl
			  << "c #minimized literals....: " << MinimizedLiterals()        << std::endl
			  << "c average lbd............: " << AvgLBD()                   << std::endl
			  << "c average cc length......: " << AvgCCLength()              << std::endl
			  << "c average dec. level.....: " << AvgDL()                    << std::endl
			  << "c average lev. cleared...: " << AvgDLclearedCA()           << std::endl
			  << "c average vars unassigned: " << AvgVarsUnassignedCA()      << std::endl
			  << "c #inprocessings.........: " << Inprocessings()            << std::endl;
	if( _antomSetting->doPreprocessing || ( ( _antomSetting->doInprocessing || _antomSetting->maxInprocess) && (Inprocessings() > 0) ) )
	  {
		std::cout << "c pre/inpro stats---------" << std::endl
				  << "c #binary constants......: " << BinaryConstants() << std::endl
				  << "c #binary equivalances...: " << BinaryEquivalences()<< std::endl
				  << "c #upla constants........: " << UplaConstants() << std::endl
				  << "c #upla equivalances.....: " << UplaEquivalences() << std::endl
				  << "c #constants by SAT check: " << SatConstants() << std::endl
				  << "c #variable eliminations.: " << VariableEliminations() << std::endl
				  << "c literal reduction elim.: " << LiteralEliminations() << std::endl
				  << "c #blocked clauses.......: " << BlockedClauses() << std::endl
				  << "c #hidden tautologies....: " << HiddenTautologies() << std::endl
				  << "c #hidden subsumptions...: " << HiddenSubsumptions() << std::endl
				  << "c #monotone variables....: " << MonotoneVariables() << std::endl
				  << "c #dc variables..........: " << DcVariables() << std::endl
				  << "c #subsumed clauses......: " << SubsumedClauses() << std::endl
				  << "c #selfsubsumed literals.: " << SelfsubsumedLiterals() << std::endl
				  << "c #bva variables.........: " << BvaVariables() << std::endl
				  << "c literal reduction bva..: " << BvaLiterals() << std::endl
				  << "c #vivify subsumptions...: " << VivifySubsumptions() << std::endl
				  << "c #vivify units..........: " << VivifyUnits() << std::endl
				  << "c #vivify literal diff...: " << VivifyDiff() << std::endl
				  << "c #unit propagations.....: " << UnitPropagations() << std::endl
				  << "c runtime upla...........: " << RuntimeUPLA() << std::endl
				  << "c runtime subsumptions...: " << RuntimeSubsumption() << std::endl
				  << "c runtime varEliminations: " << RuntimeVarElim() << std::endl
				  << "c runtime bce............: " << RuntimeBCE() << std::endl
				  << "c runtime hte............: " << RuntimeHTE() << std::endl
				  << "c runtime bva............: " << RuntimeBVA() << std::endl
				  << "c runtime vivification...: " << RuntimeVivify() << std::endl
				  << "c pre/inpro time.........: " << RuntimePrepro() << std::endl
				  << "c ------------------------" << std::endl;
	}
  }

// Weigthed maxsat
// ____________________________________________________________________________________________________________________
uint32_t Antom::MaxSolveWeightedPartial(int64_t& optimum)
{
  std::vector<uint32_t> externalAssumptions;

  if (_antomSetting->networkType==WARNERS)
      return MaxSolveWarnersWeightedPartial(externalAssumptions, optimum);
  else
      return MaxSolveWeightedPartial(externalAssumptions, optimum);
}

// ____________________________________________________________________________________________________________________
uint32_t Antom::MaxSolveWarnersWeightedPartial(const std::vector<uint32_t> & /* externalAssumptions*/, int64_t &optimum)
{
    int comp=1;
    _optimum = &optimum;
    _antomSetting->application = WEIGHTEDMAXSAT;

    _timeVariables = new TimeVariables();
    struct rusage resources;
    getrusage(RUSAGE_SELF, &resources);
    double timeS = static_cast<double>(resources.ru_utime.tv_sec + 1.e-6 ) * static_cast<double>( resources.ru_utime.tv_usec );
    //double timeC = timeS;
    _control->SetStartTime(timeS);
    _control->SetSumTime(true);
    uint32_t currentresult(ANTOM_UNKNOWN);

    std::cout << "c clauses before.........: " << Clauses() << std::endl;
    std::cout << "c #softclauses...........: " << _softClauses.size() << std::endl;
    std::cout << "c variables before.......: " << Variables() << std::endl;
    std::cout << "c #sum of SoftWeights....: " << _sumOfSoftWeights << std::endl;

    currentresult = Solve();
    if (currentresult != ANTOM_SAT)
    {
        return currentresult;
    }

    TimeMeasurement timeSolvedFirst(&_timeVariables->solvedFirst);

    CalculateOverallOptimum(_satWeight, true);
    int64_t lastOptimum = *_optimum;

    std::cout << "c SATWeight solved first.: " << _satWeight << std::endl;
    if (_satWeight == _sumOfSoftWeights)
        return ANTOM_SAT;

//    vec<Lit> lits;
//	  vec<Lit> linkingVar;
    _clausesBefore = Clauses();
    _binaryClausesBefore = CurrentBinaryClauses();
    _ternaryClausesBefore = CurrentTernaryClauses();
    uint32_t naryClausesBefore = CurrentNaryClauses();
    uint32_t variablesBefore = Variables();

    std::vector<uint32_t> lits;
    std::vector<uint32_t> linkingLit;
    GenWarners0(comp, lits, linkingLit);

    std::cout << "c #clauses of coding.....: " << Clauses() - _clausesBefore << std::endl;
    std::cout << "c    #binary of coding...: " << CurrentBinaryClauses() - _binaryClausesBefore << std::endl;
    std::cout << "c    #ternary of coding..: " << CurrentTernaryClauses() - _ternaryClausesBefore << std::endl;
    std::cout << "c    #nary of coding.....: " << CurrentNaryClauses() - naryClausesBefore << std::endl;
    std::cout << "c #variables of coding...: " << Variables() - variablesBefore << std::endl;
    std::cout << "c linkingLit.size()......: " << linkingLit.size() << std::endl;


    std::vector<long long int> cc; // cardinality constraints

    while ( true ) {

        currentresult = Solve();
        if (currentresult == ANTOM_UNSAT)
        {
            currentresult = ANTOM_SAT;
            break;
        }
        else if ( currentresult == ANTOM_UNKNOWN )
        {
            ExitTimeout();
            return ANTOM_UNKNOWN;
        }

        CalculateOverallOptimum(_satWeight, true);
        if (lastOptimum < *_optimum)
            return ANTOM_UNKNOWN;

        if (_satWeight == _sumOfSoftWeights)
        {
            currentresult = ANTOM_SAT;
            break;
        }
        int64_t opti = *_optimum;

//        lessthan(linkingLit, opti, cc, lits);
        assert(opti > 0);

        // Warners encoding
        std::vector<long long int> cls;
        cls.clear();

          opti--;
          if (opti%2 == 0)
              cls.push_back(1);

          opti = opti/2;
          int cnt = 1;
          long long int pos = 0x0002LL;
          while (opti > 0)
          {
            if (opti%2 == 0)
            {
                cls.push_back(pos);
            }
            //    else if (cls.size() == 0) cls.push(pos);
            else
            {
                for(size_t i = 0; i < cls.size(); i++)
                    cls[i] = cls[i] | pos;
            }
            pos = pos << 1;
            opti = opti/2;
            cnt++;
          }
          for(size_t i = cnt; i < linkingLit.size(); i++)
          {
            cls.push_back(pos);
            pos = pos << 1;
          }
          for(size_t i = 0; i < cls.size(); i++)
          {
            long long int x = cls[i];
            bool found = false;
            for(size_t j = 0; j < cc.size(); j++)
            {
              if (x == cc[j])
              {
                found = true;
                break;
              }
            }
            if (!found)
            {
              cc.push_back(x); // koshi 2013.10.04
              lits.clear();
              int j = 0;
              while (x > 0)
              {
                // is this a test if the last bit is true?
                if ((x & 0x0001L) == 0x0001L)
                {
                  lits.push_back(linkingLit[j]^1);
                }
                x = x >> 1;
                j++;
              }
              AddClause(lits);
            }
          }

        lastOptimum = *_optimum;
    }

    std::cout << "c #solver calls..........: " << _satSolverCalls << std::endl;
    optimum = optimum * _greatestCommonDivisor;
    return ANTOM_SAT;




}

//void genWarners0
//              (
//              vec<long long int>& weights,    <- weights of _softclauses, don't needed
//              vec<Lit>& blockings,            <- RelaxLits of _softclauses, don't needed
//              long long int max,              <- sum of weights
//              long long int k,                <- o-value, (*_optimum [* _greatestCommonDivisor])
//              int comp,                       <- in competition mode comp=1; Variants of SAT-encodings for Cardinality Constraints for warn 0,1,2,10,11
//              Solver& S,                      <- _antom, don't needed
//              vec<Lit>& lits,                 <- softclausevector
//              vec<Lit>& linkingVar            <- softclausevector
//              )

//          genWarners0(weights, blockings, max, k, comp, S, lits, linkingVar)
void Antom::GenWarners0(int comp, std::vector<uint32_t>& lits, std::vector<uint32_t>& linkingLit)
{
  // koshi 20140109
  std::cout << "c Warners' encoding for Cardinality Constraints" << std::endl;


  int logk = static_cast<int>( ceil( log2( *_optimum) ) );
//  std::cout << logk << " calculated c++ like" << std::endl;

//  getting log_2 of optimum
//  logk = 1;
//  uint64_t k = *_optimum;
//  // rightshift of k
//  while ((k >>= 1) > 0)
//      logk++;
//  std::cout << logk << " calculated as before" << std::endl;

//  Lit zero = mkLit(S.newVar());
//  lits.clear();
//  lits.push_back(~zero);
//  AddClause(lits);

  uint32_t zero = NewVariable() << 1;
  //  SetDontTouch(zero);

  // TOBI: WHY IN HAVENS SAKE SHOULD I GENERATE A NEW VARIABLE AND ADD IT DIRECTLY AS UNIT CLAUSE...
  // IN QMAXSAT IT IS A SPECIAL ADD CLAUSE - COMMENT FROM QMAXSAT:
  //Add a clause to the solver without making superflous internal copy. Will change the passed vector 'ps'.
  AddUnit( zero^1 );

  GenWarners(_softClauses, _sumOfSoftWeights, logk, comp, zero, lits, linkingLit);
//  genWarners(weights,blockings, max,logk, comp, S, zero,lits,linkingVar);
}


// koshi 2013.03.25
// Parallel counter
// koshi 2013.04.16, 2013.05.23
//void genWarners(
//          vec<long long int>& weights,    <- all softWeights -> in _softclauses, don't needed
//          vec<Lit>& blockings,            <- relaxing literals
//          long long int max,              <- first time max is sum of softWeights; max is answer in main;
//          int k,                          <- o-value, (*_optimum * _greatestCommonDivisor)
//          int comp,                       <- in competition mode comp=1; Variants of SAT-encodings for Cardinality Constraints for warn 0,1,2,10,11
//          Solver& S,                      <- _antom, don't needed
//          const Lit zero,                 <- first literal added in genWarners0
//          vec<Lit>& lits,                 <- _probably softclausevector...
//          vec<Lit>& linkingVar
//          ) {
/**
 * @brief Antom::GenWarners - Code from QMaxSAT
 *                      Parallel counter.
 *                      Is a divide and conquer algorithm.
 * @param softCl        current soft clause vector.
 * @param sumOfSW       sum of softweights of soft clause vector softCl.
 * @param logk          log2 of *_optimum.
 * @param comp          computation mode, 1 is as in qmaxSAT Auto Mode
 * @param zero          literal generated in GenWarners0, and added as negated unitClause.
 * @param lits          vector of literals.
 * @param linkingVar    vector of linkingLiterals ( RelaxLits and zero pushed into it ).
 *
 */
void Antom::GenWarners( std::vector< SoftClause* >  softCl,
                        uint64_t /* sumOfSW */,
                        int logk,
                        int comp,
                        uint32_t zero,
                        std::vector<uint32_t>& lits,
                        std::vector<uint32_t>& linkingLit)
{

    linkingLit.clear();

//      linkingVar.clear();
//    TOBI: dvar is only used to give the solver the initial assignment. Not used like this in antom!
//    bool dvar = (comp == 11) ? false : true;

//      if (weights.size() == 1) {
    if (softCl.size() == 1)
    {
//        long long int weight = weights[0];
        uint64_t weight = softCl[0]->weight;
//        vec<bool> pn;
        std::vector<bool> pn;
//        pn.clear();
        pn.clear();

        // create binary weight representation
        while (weight > 0) {
          if (weight%2 == 0)
              pn.push_back(false);
          else
              pn.push_back(true);
          weight /= 2;
        }
        for(size_t i = 0; i < pn.size(); i++) {
          if (pn[i])
//              linkingVar.push_back(blockings[0]);
              linkingLit.push_back(softCl[0]->relaxationLit);
          else
              linkingLit.push_back(zero);
        }
        pn.clear();
    }
    else if (softCl.size() > 1)
    {
//        long long int weightL = 0;
//        long long int weightR = 0;
        uint64_t weightL = 0;
        uint64_t weightR = 0;

//      vector of weights
//        vec<long long int> weightsL, weightsR;

        // vector of literals
//        vec<Lit> blockingsL, blockingsR;
//        std::vector<uint32_t> blockingsL, blockingsR;
        std::vector< SoftClause* > softClausesL, softClausesR;

//        long long int half = max/2;
//        uint64_t half = sumOfSW / 2;
//        std::cout << "Max: " << sumOfSW << "  Half: " << half << std::endl;

        std::for_each(softCl.begin(), softCl.begin() + softCl.size()/2, [&](SoftClause* SC)
        {
            weightL += SC->weight;
            softClausesL.push_back(SC);
        });

        std::for_each(softCl.begin() + softCl.size()/2, softCl.end(), [&](SoftClause* SC)
        {
            weightR += SC->weight;
            softClausesR.push_back(SC);
        });

//        assert(weightL+weightR == sumOfSW);

        std::vector<uint32_t> alpha;
        std::vector<uint32_t> beta;

        uint32_t sum = NewVariable() << 1;
        uint32_t carry = NewVariable() << 1;

        GenWarners(softClausesL, weightL, logk, comp, zero, lits, alpha);
        GenWarners(softClausesR, weightR, logk, comp, zero, lits, beta);

        bool lessthan = (alpha.size() < beta.size());
        std::vector<uint32_t> &smalls = lessthan ? alpha : beta;
        std::vector<uint32_t> &larges = lessthan ? beta : alpha;

        assert(smalls.size() <= larges.size());

        GenWarnersHalf(smalls[0], larges[0], carry, sum, comp, lits);

        linkingLit.push_back(sum);

        size_t i = 1;
        uint32_t carryN;
        for(; i < smalls.size(); i++) {
          sum = NewVariable() << 1;
          carryN = NewVariable() << 1;
          GenWarnersFull(smalls[i], larges[i], carry, carryN, sum, comp, lits);
          linkingLit.push_back(sum);
          carry = carryN;
        }
        for(; i < larges.size(); i++) {
//          sum = mkLit(S.newVar(true,dvar));
//          carryN = mkLit(S.newVar(true,dvar));
          sum = NewVariable() << 1;
          carryN = NewVariable() << 1;
          GenWarnersHalf(larges[i], carry, carryN, sum, comp, lits);
          linkingLit.push_back(sum);
          carry = carryN;
        }
        linkingLit.push_back(carry);
        alpha.clear();
        beta.clear();
    }
    size_t lsize = linkingLit.size();
    for (size_t i = logk; i < lsize; i++) { // koshi 2013.05.27
        //    printf("shrink: k = %d, lsize = %d\n",k,lsize);
        lits.clear();
        lits.push_back(linkingLit[i]^1);
//        AddClause(lits);
        AddClause(lits);
    }
//      for (int i = logk; i < lsize; i++) linkingLit.shrink(1); // koshi 2013.05.27
      // koshi 2013.05.27
    for (size_t i = logk; i < lsize; i++)
        linkingLit.pop_back();
}


/*
  Cardinality Constraints:
  Joost P. Warners, "A linear-time transformation of linear inequalities
  into conjunctive normal form",
  Information Processing Letters 68 (1998) 63-69
 */

// koshi 2013.04.16
//void genWarnersHalf(Lit& a, Lit& b, Lit& carry, Lit& sum, int comp,
//		       Solver& S, vec<Lit>& lits) {

void Antom::GenWarnersHalf( uint32_t& a,
                            uint32_t& b,
                            uint32_t& carry,
                            uint32_t& sum,
                            int comp,
                            std::vector<uint32_t>& lits) {

  // carry
  lits.clear();
  lits.push_back(a^1);
  lits.push_back(b^1);
  lits.push_back(carry);
  AddClause(lits);

  // sum
  lits.clear();
  lits.push_back(a);
  lits.push_back(b^1);
  lits.push_back(sum);
  AddClause(lits);

  lits.clear();
  lits.push_back(a^1);
  lits.push_back(b);
  lits.push_back(sum);
  AddClause(lits);

  //
  if (comp == 1 || comp == 2) {
    lits.clear();
    lits.push_back(carry);
    lits.push_back(sum);
    lits.push_back(a^1);
    AddClause(lits);

    lits.clear();
    lits.push_back(carry);
    lits.push_back(sum);
    lits.push_back(b^1);
    AddClause(lits);
  }
  if (comp == 2) {
    lits.clear();
    lits.push_back(carry^1);
    lits.push_back(sum^1);
    AddClause(lits);

    lits.clear();
    lits.push_back(carry^1);
    lits.push_back(sum);
    lits.push_back(a);
    AddClause(lits);

    lits.clear();
    lits.push_back(carry^1);
    lits.push_back(sum);
    lits.push_back(b);
    AddClause(lits);
  }
  // koshi 2013.05.31
  if (comp == 10 || comp == 11) { // [Warners 1996]
    // carry
    lits.clear();
    lits.push_back(a);
    lits.push_back(carry^1);
    AddClause(lits);

    lits.clear();
    lits.push_back(b);
    lits.push_back(carry^1);
    AddClause(lits);
    // sum
    lits.clear();
    lits.push_back(a^1);
    lits.push_back(b^1);
    lits.push_back(sum^1);
    AddClause(lits);

    lits.clear();
    lits.push_back(a);
    lits.push_back(b);
    lits.push_back(sum^1);
    AddClause(lits);
  }
}

// koshi 2013.04.16
void Antom::GenWarnersFull(     uint32_t& a,
                                uint32_t& b,
                                uint32_t& c,
                                uint32_t& carry,
                                uint32_t& sum,
                                int comp,
                                std::vector<uint32_t>& lits)
{
  // carry
  lits.clear();
  lits.push_back(a^1); lits.push_back(b^1); lits.push_back(carry); AddClause(lits);
  lits.clear();
  lits.push_back(a^1); lits.push_back(c^1); lits.push_back(carry); AddClause(lits);
  lits.clear();
  lits.push_back(b^1); lits.push_back(c^1); lits.push_back(carry); AddClause(lits);
  // sum
  lits.clear();
  lits.push_back(a); lits.push_back(b); lits.push_back(c^1); lits.push_back(sum);
  AddClause(lits);
  lits.clear();
  lits.push_back(a); lits.push_back(b^1); lits.push_back(c); lits.push_back(sum);
  AddClause(lits);
  lits.clear();
  lits.push_back(a^1); lits.push_back(b); lits.push_back(c); lits.push_back(sum);
  AddClause(lits);
  lits.clear();
  lits.push_back(a^1); lits.push_back(b^1); lits.push_back(c^1); lits.push_back(sum);
  AddClause(lits);
  if (comp == 1 || comp == 2) {
    lits.clear();
    lits.push_back(carry); lits.push_back(sum); lits.push_back(a^1); AddClause(lits);
    lits.clear();
    lits.push_back(carry); lits.push_back(sum); lits.push_back(b^1); AddClause(lits);
    lits.clear();
    lits.push_back(carry); lits.push_back(sum); lits.push_back(c^1); AddClause(lits);
  }
  if (comp == 2) {
    lits.clear();
    lits.push_back(carry^1); lits.push_back(sum^1); lits.push_back(a); AddClause(lits);
    lits.clear();
    lits.push_back(carry^1); lits.push_back(sum^1); lits.push_back(b); AddClause(lits);
    lits.clear();
    lits.push_back(carry^1); lits.push_back(sum^1); lits.push_back(c); AddClause(lits);
  }
  // koshi 2013.05.31
  if (comp == 10 || comp == 11) {// [Warners 1996]
    // carry
    lits.clear();
    lits.push_back(a); lits.push_back(b); lits.push_back(carry^1); AddClause(lits);
    lits.clear();
    lits.push_back(a); lits.push_back(c); lits.push_back(carry^1); AddClause(lits);
    lits.clear();
    lits.push_back(b); lits.push_back(c); lits.push_back(carry^1); AddClause(lits);
    // sum
    lits.clear();
    lits.push_back(a); lits.push_back(b); lits.push_back(c); lits.push_back(sum^1);
    AddClause(lits);
    lits.clear();
    lits.push_back(a^1); lits.push_back(b^1); lits.push_back(c); lits.push_back(sum^1);
    AddClause(lits);
    lits.clear();
    lits.push_back(a^1); lits.push_back(b); lits.push_back(c^1); lits.push_back(sum^1);
    AddClause(lits);
    lits.clear();
    lits.push_back(a); lits.push_back(b^1); lits.push_back(c^1); lits.push_back(sum^1);
    AddClause(lits);
  }
}

// ____________________________________________________________________________________________________________________
uint32_t Antom::MaxSolveWeightedPartial(const std::vector<uint32_t>& /* externalAssumptions */, int64_t& optimum)
{
    // to change optimum value with another function called from cascade and bucket
    _optimum = &optimum;
    _antomSetting->application = WEIGHTEDMAXSAT;

    _timeVariables = new TimeVariables();
    struct rusage resources;
    getrusage(RUSAGE_SELF, &resources);
    double timeS = static_cast<double>(resources.ru_utime.tv_sec + 1.e-6 ) * static_cast<double>( resources.ru_utime.tv_usec );
    //double timeC = timeS;
    _control->SetStartTime(timeS);
    _control->SetSumTime(true);
    uint32_t currentresult(ANTOM_UNKNOWN);

    //uint32_t averageBucketEntries(0);

#ifndef NDEBUG
    uint64_t softWeights = 0;
    std::for_each(_softClauses.begin(), _softClauses.end(), [&](SoftClause* SC) { softWeights += SC->weight; });
    assert(softWeights == _sumOfSoftWeights);
#endif

    _clausesBefore = Clauses();
    _binaryClausesBefore = CurrentBinaryClauses();
    _ternaryClausesBefore = CurrentTernaryClauses();
    std::cout << "c clauses before.........: " << _clausesBefore << std::endl;
    std::cout << "c #softclauses...........: " << _softClauses.size() << std::endl;
    if (!_antomSetting->solveAtFirst)
        std::cout << "c #sum of SoftWeights....: " << _sumOfSoftWeights << std::endl;

    if (_antomSetting->onlyByTares || _antomSetting->mcDivideStrategy != SOLVEINNORMALCASCADEMODE)
    {
        _antomSetting->encodeStrategy = ENCODEONLYIFNEEDED;
    }

    if (_antomSetting->solveAtFirst)
    {
        TimeMeasurement timeSolvedFirst(&_timeVariables->solvedFirst);
        currentresult = Solve();
        if (currentresult == ANTOM_UNSAT)
        {
            return ANTOM_UNSAT;
        }
        else if (currentresult == ANTOM_UNKNOWN)
        {
            return ANTOM_UNKNOWN;
        }

        CalculateOverallOptimum(_satWeight, true);
        std::cout << "c #sum of SoftWeights....: " << _sumOfSoftWeights << std::endl;
        std::cout << "c SATWeight solved first.: " << _satWeight << std::endl;
    }
    if (_satWeight == _sumOfSoftWeights)
        return ANTOM_SAT;

    // calc and set all trivially conflicting softclauses
    // TOASK: maybe better seperately for each bucket - to know boundaries of bucket!
    // together with mode of solving bucket parts to get max
    CheckAllWeightedConflictingSoftclauses();

    // Do preprocessing?
    if ( _antomSetting->doPreprocessing != NOPREPRO )
    {
        // (Incrementally) preprocess formula without sorter
        // TOASK: Probably not working with new structure!
        if (Preprocess( PREPROCESS ) != ANTOM_UNKNOWN)
        {
            if( _antomSetting->incrementalMode > 0 )
            {
                InvalidateSoftClauses();
            }
            _control->SetSumTime(false);
            return ANTOM_UNSAT;
        }

        if ( _control->GetTimeOutReached() || _control->GetMemOutReached() )
        {
            ExitTimeout();
            return ANTOM_UNKNOWN;
        }
    }

    if (_antomSetting->mcDivideStrategy != SOLVEINNORMALCASCADEMODE)
    {

        _antomSetting->encodeStrategy = ENCODEONLYIFNEEDED;
        _mainMultipleCascade = new MultipleCascade(this, _antomSetting->onlyByTares, _antomSetting->tareCascadeOnlyByTares, _antomSetting->cascadeDivider, _antomSetting->maxBucketSize, _antomSetting->nOfCasc, _antomSetting->interimResult, _sumOfSoftWeights);
//        std::cout << "Cascade Divider: " << _cascadeDivider << std::endl;
        if(_mainMultipleCascade->DivideAndConnectCascades(_antomSetting->mcDivideStrategy, &_softClauses))
        {
            currentresult = _mainMultipleCascade->Solve();

            if (currentresult != ANTOM_UNKNOWN)
                currentresult = ANTOM_SAT;
        }
    } else
    {
        _antomSetting->cascadeDivider = 0;
    }

    // no else if, because if _cascadeDivider returns only one cascade then _cascadeDivider == 0
    if (_antomSetting->cascadeDivider == 0)
    {
        _mainMultipleCascade = NULL;
        _mainCascade = new Cascade(this, nullptr, _antomSetting->onlyByTares);

        _mainCascade->Fill(&_softClauses, _antomSetting->partitionStrategy, _antomSetting->encodeStrategy);

    //if (_antomSetting->verbosity > 2 && _antomSetting->base > 2)
    //{
        //std::cout << "How often is Softclause in Bucket, because of base " << _baseMode << std::endl;
        //DumpBucketStructure();
    //}

    //Answered: why at this position - is it possible to move it after encoding?
	//ANSWER: No
    std::vector< uint32_t > satconstcandidates;
    if( _antomSetting->satconst > 0 )
    {
        for ( uint32_t v = 1; v <= Variables(); ++v )
        {
            satconstcandidates.push_back(v);
        }
#ifndef NDEBUG
        bool rst = FindConstantsWithSat(satconstcandidates, _antomSetting->satconst==1); assert(rst);
#else
        FindConstantsWithSat(satconstcandidates, _antomSetting->satconst==1);
#endif
    }

    if (_antomSetting->cascadeDivider == 0)
    {
        if (_antomSetting->encodeStrategy == ENCODEALL)
        {
            _clausesBefore = Clauses();
            _binaryClausesBefore = CurrentBinaryClauses();
            _ternaryClausesBefore = CurrentTernaryClauses();
        }
        if (! _mainCascade->Encode())
        {
            ExitTimeout();
            return ANTOM_UNKNOWN;
        }
        if (_antomSetting->encodeStrategy == ENCODEALL)
        {
            _addedClauses = Clauses() - _clausesBefore;
            _addedBinaryClauses= CurrentBinaryClauses() - _binaryClausesBefore;
            _addedTernaryClauses= CurrentTernaryClauses() - _ternaryClausesBefore;
        }
    }

//    if( _dopreprocessing == INCREMENTAL )
//    {
//        // (Incrementally) preprocess formula including sorter
//        if (Preprocess( INCREMENTAL ) != ANTOM_UNKNOWN)
//        {
//            if( _antomSetting->incrementalMode > 0 )
//            {
//                InvalidateSoftClauses();
//            }
//            _control->SetSumTime(false);
//            return ANTOM_UNSAT;
//        }
//        if ( _control->GetTimeOutReached() || _control->GetMemOutReached() )
//        {
//            ExitTimeout();
//            return ANTOM_UNKNOWN;
//        }
//    }


    // TOBI: SetDecisionStrategiesForCascadeModel!!!!
    //SetDecisionStrategiesForMaxSAT();
        currentresult = _mainCascade->Solve();
    }

    _timeVariables->DumpVariables();
    std::cout << "c #clauses of coding.....: " << _addedClauses << std::endl;
    std::cout << "c    #binary of coding...: " << _addedBinaryClauses << std::endl;
    std::cout << "c    #ternary of coding..: " << _addedTernaryClauses << std::endl;
    std::cout << "c #solver calls..........: " << _satSolverCalls << std::endl;


    if ( currentresult == ANTOM_UNKNOWN )
        {
            ExitTimeout();
            return ANTOM_UNKNOWN;
        }

    // if there is more than one cascade, we have to calculate
    // the OverallOptimum again to sum up the weights.
    //CalculateOverallOptimum(_satWeight, true);
    optimum = optimum * _greatestCommonDivisor;
    return currentresult;
}


StructureInfo Antom::AnalyzeandConvertStructure()
{
    StructureInfo structure = AnalyzeStructure();
    switch ( structure )
    {
    case ISWEIGHTEDMAXSAT:
        break;
    case ISSAT:
        return ISSAT;
    case ISMAXSAT:
        if (_minWeight > 1)
        {
            DivideAllSoftClausesByFactor(_minWeight);
            std::cout << "DivideByDivisormW: " << _minWeight << std::endl;
        }
        break;
    case CONVERTTOMAXSAT:
        ConvertFormulaToMaxSAT(_maxWeight);
        // TOASK: better - TOTELL
        // This line makes it faster with a high factor!!
        // But the value of _minWeight has to be multiplied to optimum!
        if (_minWeight > 1)
        {
            DivideAllSoftClausesByFactor(_minWeight);
            std::cout << "DivideByDivisormW: " << _minWeight << std::endl;
        }

        structure = ISMAXSAT;
        break;
    case antom::DIVIDEWEIGHTSBYDIVISOR:
        DivideAllSoftClausesByFactor(_greatestCommonDivisor);
        structure = ISWEIGHTEDMAXSAT;
        break;
    }
    // Add information to SoftClauses
    AddSoftClauseVector();
    return structure;
}

StructureInfo Antom::AnalyzeStructure()
{
    //assert(_sumOfSoftWeights > _topWeight);
    bool maxIsMin = (_maxWeight == _minWeight);
    bool withHC = Clauses() != 0;
    bool minIsSet = _minWeight != ((uint64_t) - 1);
    bool maxIsSet = _maxWeight != 0;
    // is SAT formula
    bool onlyHC = !minIsSet && !maxIsSet && withHC;
    // is MaxSAT formula -- NOT PARTIAL
    bool onlySCwithOneWeight = !_moreThanTwoWeights && maxIsMin && !withHC;
    // is partial MaxSAT formula - Set _minWeight to one
    bool oneHCWeightAndSCWeight = maxIsMin && withHC;
    // if sum (_minWeight) < _maxWeight --> is partial MaxSAT formula
    bool onlySCwithTwoWeights = !_moreThanTwoWeights && !maxIsMin && minIsSet && maxIsSet && !withHC;
    // Calculate the greatest common divisor of all softclauseweights.
    uint64_t greatestCommonDivisor = _minWeight;

    for (uint32_t ind = 0; ind < _softClauses.size(); ++ind)
    {
        if (greatestCommonDivisor == 1)
            break;
        greatestCommonDivisor = GreatestCommonDivisor(greatestCommonDivisor, _softClauses[ind]->weight);
    }
    _greatestCommonDivisor = greatestCommonDivisor;


    if (_antomSetting->verbosity > 2)
    {
        std::cout << std::setw(30) << "_topWeight " << _topWeight << std::endl;
        std::cout << std::setw(30) << "_minWeight " << _minWeight << std::endl;
        std::cout << std::setw(30) << "_maxWeight " << _maxWeight << std::endl;
        std::cout << std::setw(30) << "_moreThanTwoWeights " << _moreThanTwoWeights << std::endl;
        std::cout << std::setw(30) << "maxIsMin " << maxIsMin << std::endl;
        std::cout << std::setw(30) << "withHC " << withHC << std::endl;
        std::cout << std::setw(30) << "minIsSet " << minIsSet << std::endl;
        std::cout << std::setw(30) << "maxIsSet " << maxIsSet << std::endl;
        std::cout << std::setw(30) << "onlyHC " << onlyHC << std::endl;
        std::cout << std::setw(30) << "onlySCwithOneWeight " << onlySCwithOneWeight << std::endl;
        std::cout << std::setw(30) << "oneHCWeightAndSCWeight " << oneHCWeightAndSCWeight << std::endl;
        std::cout << std::setw(30) << "onlySCwithTwoWeights " << onlySCwithTwoWeights << std::endl;
        std::cout << std::setw(30) << "parser -> _sumOfSoftWeights " << _sumOfSoftWeights << std::endl;
    }

    if(onlyHC)
    {
        return ISSAT;
    }
    if (onlySCwithTwoWeights)
    {
        uint32_t sumOfMinWeights(0);

        // the sum of all _minWeights has to be smaller than _maxWeight
        // then it can be solved like MaxSAT formula.
        for (uint32_t ind = 0; ind < _softClauses.size(); ++ind)
        {
            if (_softClauses[ind]->weight != _minWeight)
                continue;

            sumOfMinWeights += _softClauses[ind]->weight;
        }

        if (sumOfMinWeights < _maxWeight)
        {
            _sumOfSoftWeights = sumOfMinWeights;
            return CONVERTTOMAXSAT;
        }
    }
    if (oneHCWeightAndSCWeight || onlySCwithOneWeight || oneHCWeightAndSCWeight)
    {
        return ISMAXSAT;
    }
    if (greatestCommonDivisor > 1)
    {
        _greatestCommonDivisor = greatestCommonDivisor;
        std::cout << "c greatest common divisor: " << greatestCommonDivisor << std::endl;
        return antom::DIVIDEWEIGHTSBYDIVISOR;
    }
    return ISWEIGHTEDMAXSAT;
}

uint64_t Antom::GreatestCommonDivisor(uint64_t a, uint64_t b)
{
    uint64_t temp;
    while(b > 0) {
        temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

void Antom::ConvertFormulaToMaxSAT(uint64_t maxWeight)
{
    // configure all _maxWeights as hardclauses
    std::vector< SoftClause* > newSoftClauseVector;
    for (uint32_t ind = 0; ind < _softClauses.size(); ++ind)
    {
        if (_softClauses[ind]->weight != maxWeight)
        {
            newSoftClauseVector.push_back(_softClauses[ind]);
            continue;
        }

        AddClause(_softClauses[ind]->clause);
    }

    _softClauses = newSoftClauseVector;
}

void Antom::DivideAllSoftClausesByFactor(uint64_t factor)
{
    for (uint32_t ind = 0; ind < _softClauses.size(); ++ind)
    {
        _softClauses[ind]->weight /= factor;
    }
    _sumOfSoftWeights /= factor;
}

uint64_t Antom::CalculateOverallOptimum(uint64_t satWeight, bool countAgain)
{
    if (_antomSetting->verbosity > 3)
        std::cout << __PRETTY_FUNCTION__ << std::endl;


    struct rusage resources;
    // Cascade and Bucket counts the satisfied SC only for their SC's
    // for an overall result count again.
    if (countAgain)
    {
        // maybe _maxSatWeight - because it could be a local optima of this Solver call
        satWeight = CountSatisfiedSoftClauses(NULL, Model());
    }

    if (*_optimum > static_cast<int64_t>(_sumOfSoftWeights - satWeight) || *_optimum == -1)
    {
        if (_antomSetting->verbosity > 0)
        {
            getrusage(RUSAGE_SELF, &resources);
            double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
            std::cout << "c " << (timeC - _control->GetStartTime()) << "s" << std::endl;
        }
        // better actualize satWeight once at the end.
        _satWeight = satWeight;
        *_optimum = static_cast<int64_t>(_sumOfSoftWeights - satWeight);
        std::cout << "o " << *_optimum * _greatestCommonDivisor << std::endl;
        _lastModel = Model();

        if (_antomSetting->verbosity > 2)
            std::cout << std::setw(50) << "Calculated Global SATWeight: " << satWeight << std::endl;

    } else if (_antomSetting->verbosity > 2)
    {
        if (_antomSetting->verbosity > 3)
        {
            getrusage(RUSAGE_SELF, &resources);
            double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;
            std::cout << "c " << (timeC - _control->GetStartTime()) << "s" << std::endl;
        }
        std::cout << "o " << *_optimum * _greatestCommonDivisor << std::endl;
        std::cout << std::setw(50) << "Calculated Global SATWeight: " << satWeight << std::endl;
    }

    return _satWeight;
}

void Antom::SetEncode01Mode(bool val, bool val2)
{
    _antomSetting->encode01 = val;
    _antomSetting->lastPos1 = val2;
}

void Antom::SetDecStratMode(uint32_t val)
{
    assert( val < 3 );
    _antomSetting->decStrat = val;
}


void Antom::SetMoreThanTwoWeights(bool val)
{
    _moreThanTwoWeights = val;
}

void Antom::SetTopWeight(uint64_t val)
{
    assert( val > 0 );
    _topWeight = val;
}

void Antom::SetMinWeight(uint64_t val)
{
    // case only HC
    //assert( val < (uint64_t) - 1 );
    _minWeight = val;
}

void Antom::SetMaxWeight(uint64_t val)
{
    // case only HC
    //assert( val > 0 );
    _maxWeight = val;
}

void Antom::SetSumOfSoftWeights(uint64_t val)
{
    // case only HC
    //assert( val > 0 );
    _sumOfSoftWeights = val;
}

void Antom::SetBaseMode(uint32_t val)
{
    assert( val > 1 );
	_antomSetting->base = val;
}

void Antom::SetPartitionStrategy(PartitionStrategy val)
{
    //assert( val < 3 );
    _antomSetting->partitionStrategy = val;
}

void Antom::SetHeuristic(uint32_t val)
{
  assert(val < 11);
  _antomSetting->groupHeuristic = val;
}

void Antom::SetPercentOff(uint32_t val)
{
  assert(val <= 100);
  _antomSetting->percentOff = val;
}
  
  void Antom::SetPercentOffreinsert(bool val)
{
  _antomSetting->percentOffReinsert = val;
}

void Antom::SetCascadeDivider(uint32_t val)
{
  _antomSetting->cascadeDivider = val;
}

void Antom::SetEqualWeight(uint32_t val)
{
  _antomSetting->equalWeight = val;
}

void Antom::SetSolveAtFirst(bool val)
{
  _antomSetting->solveAtFirst = val;
}

void Antom::SetEncodeStrategy(EncodeStrategy val)
{
  _antomSetting->encodeStrategy = val;
}

void Antom::SetCreateGraphFile(std::string val)
{
  _antomSetting->createGraphFile = val;
}

void Antom::SetOnlyByTares(bool val)
{
  _antomSetting->onlyByTares = val;
}

void Antom::SetSolveTareCascadeOnlyByTares(bool val)
{
  _antomSetting->tareCascadeOnlyByTares = val;
}

void Antom::SetMultipleCascadeDivideStrategy(MultipleCascadeDivideStrategy val)
{
  _antomSetting->mcDivideStrategy = val;
}

void Antom::SetSepHiWeight(bool val)
{
  _antomSetting->sepHiWeight = val;
}

void Antom::SetFeatureTest(bool val)
{
  _antomSetting->featureTest = val;
}

void Antom::SetWeightPlusOne(bool val)
{
  _antomSetting->weightPlusOne = val;
}

void Antom::SetInterimResult(InterimResult val)
{
  _antomSetting->interimResult = val;
}

void Antom::SetMaxBucketSize(uint32_t val)
{
  _antomSetting->maxBucketSize = val;
}

void Antom::SetNOfCascades(uint32_t val)
{
  _antomSetting->nOfCasc = val;
}


}
