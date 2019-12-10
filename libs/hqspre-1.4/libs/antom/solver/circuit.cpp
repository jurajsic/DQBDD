
/********************************************************************************************
circuit.cpp -- Copyright (c) 2016, Sven Reimer

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

// Include antom related headers.
#include "circuit.h"
#include "helper.h"
#include <fstream>

namespace antom
{
  Circuit::Circuit(void) :
	_usedOutputVars(),
	_nodes()
  {
	// Create dummy node representing all inputs to the circuit
	CircuitNode* inputNode = new CircuitNode(antom::INPUT, 0);
	_nodes.push_back(inputNode);
  }

  Circuit::~Circuit(void)
  {
	for (CircuitNode* node : _nodes)
	  {
		delete node;
	  }
  }

  void Circuit::AddGate(const std::vector<uint32_t>& inputLits, uint32_t outputLit, antom::GateType type)
  {
	// outputlit should not be used
	assert(_usedOutputVars.find(outputLit>>1) == _usedOutputVars.end());

	CircuitNode* node = new CircuitNode(type, outputLit);

	for (uint32_t input : inputLits)
	  {
		// Search whether we have already used the input as an output
		auto result = _usedOutputVars.find(input>>1);
		if (result != _usedOutputVars.end())
		  {
			// We have already used this input as an ouput
			node->inputGates.push_back(result->second);
			// Add this gate as outputgate of found gate
			_nodes[result->second]->outputGates.push_back(_nodes.size());
		  }
		else
		  {
			node->inputGates.push_back(0);
		  }
		node->inputLits.push_back(input);
	  }
	_usedOutputVars[outputLit>>1] = _nodes.size();
	_nodes.push_back(node);
  }

  void Circuit::ToDotFile(std::string filename) const
  {
	std::cout << "Exporting to outfile : " << filename << std::endl;;
	std::string type;
	std::ofstream myfile;
	myfile.open (filename.c_str());
	myfile <<  "digraph \"" << filename << "\" {\n";
	myfile <<  "rankdir=\"LR\"\n";

	const CircuitNode* currentNode = NULL;

	// Define all nodes first
	for(uint32_t i = 0; i < _nodes.size(); i++)
	  {
		currentNode = _nodes[i];
		
		type = GateTypeToString(currentNode->type);

		// TODO:
		myfile  <<  "\"" << i
				<< "\" [label=\"" << i << "| "<< type;

		myfile << "|out: " << currentNode->outputLit;

		myfile		<< "\", shape=record ";

		if (currentNode->type == INPUT)
		  {
			myfile		<< ",style=filled, fillcolor=\"#CCCC00\"";
		  }
		else if (currentNode->outputGates.size() == 1)
		  {
			myfile		<< ",style=filled, fillcolor=\"#00CC00\"";
		  }
		myfile		<< "]\n";

	  }

	// Define all connections between nodes
	for(uint32_t i = 0; i < _nodes.size(); i++)
	  {
		currentNode = _nodes[i];
		
		for(uint32_t j = 0; j < currentNode->outputGates.size(); ++j)
		  {
			myfile <<  "\n\"" << i << "\" -> \"" << currentNode->outputGates[j] << "\"";
		  }
	  }

	myfile <<  "}\n";
	myfile.close();
  }
}
