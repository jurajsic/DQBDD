/********************************************************************************************
control.hpp -- Copyright (c) 2013-2017, Tobias Schubert, Sven Reimer

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

#ifndef ANTOM_CONTROL_H_
#define ANTOM_CONTROL_H_

// Include standard headers.
#include <cassert>
#include <cstddef>
#include <vector>
#include <sys/resource.h>
#include <fstream>
#include <unistd.h>
#include "settings.h"


namespace antom
{
  // The "ThreadData" class.
  struct ThreadData
  {

  public:

    // Constructor.
    ThreadData(void) : 
      decisionStack(NULL),
      dsStart(0),
      dsEnd(0)
    { }

    // Decision stack related data.
    uint32_t* decisionStack;
    uint32_t dsStart; 
    uint32_t dsEnd; 

  private: 
    
    // Copy constructor.
    ThreadData (const ThreadData&);

    // Assignment operator.
    ThreadData& operator = (const ThreadData&);

  };

  // The "Control" class.
  class Control
  {

  public:
    // Constructor.
  explicit Control(Settings* setting) : 
	_done(false),
	  _setting(setting),
      _threadData(1),
      _unitClauses(1),
	  _timeS(0.0),
	  _timeOutReached(false),
	  _memOutReached(false),
	  _sumTime(false),
	  _timeSetted(false),
	  _virtualMemory(0.0),
	  _residentSetMemory(0.0),
	  _result(ANTOM_UNKNOWN)
    {
      // Consistency check.
      assert(_setting != nullptr);   

      // Initialize "_threadData".
      for (uint32_t t = 0; t < _setting->threads; ++t)
		{ _threadData[t] = new ThreadData; }
    }

    // Destructor.
    ~Control(void)
    {
      // Delete "_threadData".
      for (uint32_t t = 0; t < _setting->threads; ++t)
		{ delete _threadData[t]; }
    }

	void InstanceReset(void)
	{
	  _done = false;
	  _timeOutReached = false;
	  _memOutReached = false;
	  _sumTime = false;
	  _virtualMemory = 0.0;
	  _timeSetted = false;
	  _residentSetMemory = 0.0;
	  _result = ANTOM_UNKNOWN;
	}

	void Reset(void)
	{
	  InstanceReset();
	  _timeS = 0.0;
	  _timeSetted = false;
	}

    // Returns whether the search process can be stopped or not.
    bool Done(void) const 
	{ return _done; }

    // Sets "_done" to TRUE.
    void SetDone(void) 
	{ _done = true; }

    // Sets "_done" to FALSE.
    void ResetDone (void) 
	{ _done = false; }

	void SetStartTime(double startTime)
	{ 
	  if( !_sumTime && !_timeSetted ) 
		{ _timeS = startTime; _timeSetted = true; } 
	}

	double GetStartTime(void) const 
	{ return _timeS; }

	void SetSumTime(bool val)
	{ _sumTime = val; }

	bool GetSumTime(void) const 
	{ return _sumTime; }

	void SetTimeSetted(bool val)
	{ _timeSetted = val; }

	bool GetTimeSetted(void) const 
	{ return _timeSetted; }

	bool ReachedTimeout(void) 
	{ 
	  if( _setting->cpuLimit == 0.0 )
		{ return false; }

	  struct rusage resources;

	  getrusage(RUSAGE_SELF, &resources);
	  double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double) resources.ru_utime.tv_usec;

	  _timeOutReached = ( (timeC - _timeS ) > _setting->cpuLimit );
	  if( _timeOutReached && _setting->verbosity > 0 )
		{
		  std::cout << "c reached timeout of " << (timeC - _timeS ) << "s. Limit was " << _setting->cpuLimit << "s" << std::endl;
		}
	  return _timeOutReached; 
	}

	bool GetTimeOutReached(void) const
	{
	  return _timeOutReached;
	}

	bool GetMemOutReached(void) const
	{
	  return _memOutReached;
	}

	double GetResidentSetMemory(void) const
	{
	  return _residentSetMemory;
	}

	double GetVirtualMemory(void) const
	{
	  return _virtualMemory;
	}

	bool ReachedMemout(void)
	{ 
	  if( _setting->memLimit == 0.0 )
		{ return false; }

	  ProcessMemUsage();
	  _memOutReached = ( _residentSetMemory > _setting->memLimit );
	  if( _memOutReached && _setting->verbosity > 0 )
		{
		  std::cout << "c reached memory limit " << _residentSetMemory << " MB. Limit was " << _setting->memLimit << " MB" << std::endl;
		}
	  return _memOutReached; 
	}

	bool ReachedLimits(void)
	{
	  return (ReachedTimeout() || ReachedMemout() );
	}

	void SetExtendedResult(uint32_t result)
	{
	  _result = result;
	}
	
	// Get extended solver result 
	uint32_t GetExtendedResult(void) const
	{ return _result; }


    // Used by the threads to import unit clauses found by other threads.

    // Used by the threads to export all of their assignments made on decision level 0.
    void ExportUnitClauses(uint32_t* decisionStack, uint32_t start, uint32_t end, uint32_t id)
    { 
      assert(id < _setting->threads);
      _threadData[id]->decisionStack = decisionStack; 
      _threadData[id]->dsStart = start;
      _threadData[id]->dsEnd = end; 
    }

    uint32_t* ImportUnitClauses(uint32_t id)
    {
      assert(id < _setting->threads); 
      _unitClauses[id].clear(); 
      for (uint32_t t = 0; t < _setting->threads; ++t)
		{
		  if (t != id)
			{
			  for (uint32_t l = _threadData[t]->dsStart; l < _threadData[t]->dsEnd;  ++l)
				{ _unitClauses[id].push_back(_threadData[t]->decisionStack[l]); }
			}
		}
      _unitClauses[id].push_back(0); 
      return &_unitClauses[id][0]; 
    }

  private:

	void ProcessMemUsage(void)
	{
	  _virtualMemory = 0.0;
	  _residentSetMemory = 0.0;

	  struct rusage resources;
	  getrusage(RUSAGE_SELF, &resources);
	  _residentSetMemory = static_cast<double>(resources.ru_maxrss/1024);
	}

    // A flag indicating whether the search process can be stopped or not.
    bool _done;

	Settings* _setting;
  
    // For each thread we store a pointer to some decision stack related data.
    std::vector<ThreadData*> _threadData; 

    // For each thread we maintain a vector, storing all unit clauses found by other threads. 
    std::vector<std::vector<uint32_t> > _unitClauses; 

	// Start time of the SAT core
	double _timeS;

	bool _timeOutReached;
	bool _memOutReached;

	// Do we have to accumulate run time?
	bool _sumTime;

	// Starttime already set?
	bool _timeSetted;

	// Memory usage
	double _virtualMemory;
	double _residentSetMemory;

	// The SAT solver result
	uint32_t _result;
  };
}

#endif
