/********************************************************************************************
circuit.h -- Copyright (c) 2016, Sven Reimer

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

#ifndef ANTOM_CIRCUIT_H_
#define ANTOM_CIRCUIT_H_

#include <vector>
#include <string>
#include <map>

#include "helper.h"
#include "antombase.h"

namespace antom
{
  // Datastructure for managing circuit infrastructure
  class Circuit
  {
  public:
	// Constructor
	Circuit(void);
	~Circuit(void);

	void AddGate(const std::vector<uint32_t>& inputLits, uint32_t outputLit, antom::GateType type);

	void ToDotFile(std::string filename) const;

  private:
	struct CircuitNode
	{
	CircuitNode(antom::GateType t, uint32_t outLit) :
	  type(t),
		inputLits(),
		outputLit(outLit),
		inputGates(),
		outputGates()
	  {}

	  std::string GateTypeToString;
	  
	  antom::GateType type;
	  std::vector<uint32_t> inputLits;
	  uint32_t outputLit;
	  
	  std::vector<uint32_t> inputGates;
	  std::vector<uint32_t> outputGates;
	};

	// maps outputvar to nodeID
	std::map<uint32_t, uint32_t> _usedOutputVars;
	std::vector<CircuitNode*> _nodes;
  };

}

#endif
