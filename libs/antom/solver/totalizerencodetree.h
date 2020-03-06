/********************************************************************************************
totalizerencodetree.h -- Copyright (c) 2017, Tobias Paxian

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

#ifndef TOTALIZERENCODETREE_H
#define TOTALIZERENCODETREE_H

#include <vector>
#include <iostream>
#include <fstream>
#include "sorter.h"

namespace antom
{

struct TotalizerEncodeTree
{
    TotalizerEncodeTree(uint32_t size) :
  _encodedOutputs(size, 0),
	_size(size),
	_depth(0),
	_howOftenUsed(0),
	_maxPos(0),
	_allOutputsEncoded(false),
	_hasBeenBucketBefore(false),
	_onesEncoded(false),
	_everyNthOutput(1),
	_child1(NULL),
	_child2(NULL)
    {
        //std::cout << _size << std::endl;
    }
    // the inputs are the _encodedOutputs of the children.
    std::vector<uint32_t> _encodedOutputs;
    uint32_t _size;
    uint32_t _depth;
    uint32_t _howOftenUsed;
    uint32_t _maxPos;

    bool _allOutputsEncoded;
    bool _hasBeenBucketBefore;
    bool _onesEncoded;

    // as standard every output counts.
    uint32_t _everyNthOutput;

    TotalizerEncodeTree* _child1;
    TotalizerEncodeTree* _child2;

    void CutVectorAbove(uint32_t maxPos)
    {
        assert(_size == _encodedOutputs.size());
        std::cout << std::endl << "Size before Vector Cut ABOVE EncodeTree: " << _encodedOutputs.size() << "  SIZE: " << _size << std::endl;
        _maxPos = maxPos;
        while(_size > maxPos)
        {
            _encodedOutputs.pop_back();
            _size--;
        }
        std::cout << "Size after Vector Cut ABOVE EncodeTree: " << _encodedOutputs.size() << "  SIZE: " << _size << std::endl << std::endl;
    }

    void CutVectorBelow(uint32_t minPos)
    {
        assert(_size == _encodedOutputs.size());

        std::cout << std::endl << "Size before Vector Cut BELOW EncodeTree: " << _encodedOutputs.size() << "  SIZE: " << _size << std::endl;
        if(_size > minPos)
            _encodedOutputs.erase(_encodedOutputs.begin(), _encodedOutputs.begin()+minPos);

        _size = _encodedOutputs.size();

        std::cout << "Size after Vector Cut BELOW EncodeTree: " << _encodedOutputs.size() << "  SIZE: " << _size << std::endl << std::endl;
    }


    uint32_t ReturnOutputEncodeIfNecessary(uint32_t index, Sorter* sorter, bool encodeOnlyOnes = false)
    {
//        std::cout << std::endl << "NodeSize: " << _size;
//        if(_everyNthOutput > 1)
//            std::cout << "/" << _everyNthOutput;
//        else if (_hasBeenBucketBefore)
//            std::cout << "*" << _howOftenUsed;
//        std::cout << "  Encoded Outputs:  (  ";
//        for (uint32_t i = 0; i < _encodedOutputs.size(); i++)
//            std::cout << _encodedOutputs[i] << "  ";
//        std::cout << ")" << std::endl;

//        std::cout << "requested Index: " << index;
//        std::cout << "   _everyNthOutput; " << _everyNthOutput << "   _size: " << _size << "   index: " << index << std::endl;

        index = (_everyNthOutput * (index + 1)) - 1;

        //TOBI: WHY??
        //std::cout << "size: " << _size << "  index: " << index << "  encodedOutputs.size(): " << _encodedOutputs.size() << std::endl;
        assert(_size > index);
        assert(_size = _encodedOutputs.size());

        if (encodeOnlyOnes && !_onesEncoded && _size != 1)
        {
            _encodedOutputs[index] = sorter->TotalizerEncodeOnes(this, index, _encodedOutputs[index]);
            _onesEncoded = true;
        }
        else if (_encodedOutputs[index] != 0)
        {
            return _encodedOutputs[index];
        }
        else
        {
            _encodedOutputs[index] = sorter->TotalizerEncodeOutput(this, index);
        }




//        std::cout << "                                 One Index more Encoded!! NodeLabel:            " << _size;
//        if(_everyNthOutput > 1)
//            std::cout << "/" << _everyNthOutput;
//        else if (_hasBeenBucketBefore)
//            std::cout << "*" << _howOftenUsed;
//        std::cout << "  Encoded Outputs:  (  ";
//        for (uint32_t i = 0; i < _encodedOutputs.size(); i++)
//            std::cout << _encodedOutputs[i] << "  ";
//        std::cout << ")" << std::endl;

        return _encodedOutputs[index];
    }

    /**
     * @brief ActualizeValues
     *          After connecting two subbuckets, the values have to been actualized.
     */
    void ActualizeValues() {
        _size = _child1->_size / _child1->_everyNthOutput;
        _size += _child2->_size / _child2->_everyNthOutput;

//        std::cout << "child1: " << _child1->_size / _child1->_everyNthOutput << std::endl;
//        std::cout << "child2: " << _child2->_size / _child2->_everyNthOutput << std::endl;
//        std::cout << "ACTUALIZE VALUES SIZE: " << _size << std::endl;

        _depth = (_child1->_depth > _child2->_depth) ? _child1->_depth + 1 : _child2->_depth + 1;
        _encodedOutputs.resize(_size, 0);
    }

    // create empty tree for given size()
    /**
     * @brief CreateOutputTreeReturnMaxDepth
     *                    create empty tree; leaves contain input values.
     * @param lo          Points to first element of inputVector to look at.
     *                    A hard lower bound can be given!
     * @param hi          Points to last element + 1 of inputVector to look at.
     *                    A hard upper bound can be given!
     * @param inputVector Pointer to whole given inputVector.
     * @return            Depth of tree, root has greatest depth.
     */
    uint32_t CreateOutputTreeReturnMaxDepth(uint32_t lo, uint32_t hi, std::vector< uint32_t >* inputVector)
    {
        if( (hi-lo) > 1)
        {
            assert( hi > lo );
            uint32_t m = ((hi-lo)>>1);
            _child1 = new TotalizerEncodeTree( m );
            uint32_t depth1 = _child1->CreateOutputTreeReturnMaxDepth(lo, lo+m, inputVector);
            _child2 = new TotalizerEncodeTree( _size - m );
            uint32_t depth2 = _child2->CreateOutputTreeReturnMaxDepth(lo+m, hi, inputVector);
            _depth = (depth1 > depth2) ? depth1 : depth2;
        } else {
            assert( hi - lo == 1 );
            assert( hi > lo );
            _encodedOutputs[0] = (*inputVector)[lo];
            _allOutputsEncoded = true;
            _depth = 0;
        }
        return _depth + 1;
    }

//    uint32_t AttachTwoChildrenReturnMaxDepth(TotalizerEncodeTree* firstChild, TotalizerEncodeTree* secondChild)
//    {
//        assert(firstChild->_size == 0 || secondChild->_size == 0);

//        _child1 = firstChild;
//        _child2 = secondChild;
//        _depth = (firstChild->_depth > secondChild->_depth) ? firstChild->_depth + 1 : secondChild->_depth + 1;

//        return _depth;
//    }


    /**
     * @brief DumpOutputTree
     *          The output is an Trivial Graph format output - *.tgf!!
     *          To visualize it, use a program like yed.
     *          At first print out all existing node numbers with their label.
     *          Followed by the # symbol.
     *          Then cout all connections between the nodes.
     */
    void DumpOutputTree(std::string filename, bool labelWithOutputs = false)
    {
        std::ofstream out(filename);
        std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
        std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

        DumpNodes(1, labelWithOutputs);
        std::cout << "#" << std::endl;
        DumpConnections(1);
        std::cout.rdbuf(coutbuf); //reset to standard output again
        std::cout << filename << " file written." << std::endl;
    }

    /**
     * @brief DumpNodes
     *          Cout all existing node numbers with their label.
     * @param nodeNumber
     *          Call with next not yet given node number.
     * @return  Highest node subNode number.
     */
    uint32_t DumpNodes(uint32_t nodeNumber, bool labelWithOutputs)
    {

        if (_size == 1)
        {
            std::cout << nodeNumber << " " << _encodedOutputs[0] << std::endl;
            return nodeNumber;
        }
        std::cout << nodeNumber << " ";
        if(labelWithOutputs)
        {        std::cout << "(";
            for (uint32_t i = 0; i < _encodedOutputs.size(); i++)
            {
                if (i != _encodedOutputs.size() - 1)
                    std::cout << _encodedOutputs[i] << ",";
                else
                    std::cout << _encodedOutputs[i] << ")";
            }
        } else
        {
            std::cout << _size;
            if(_everyNthOutput > 1)
                std::cout << "/" << _everyNthOutput;
            else if (_hasBeenBucketBefore)
                std::cout << "*" << _howOftenUsed;
        }
        std::cout << std::endl;
        uint32_t nodeNumber1 = _child1->DumpNodes(nodeNumber + 1, labelWithOutputs);
        uint32_t nodeNumber2 = _child2->DumpNodes(nodeNumber1 + 1, labelWithOutputs);
        return nodeNumber2;
    }

    /**
     * @brief DumpConnections
     *          Cout all possible connections between nodes.
     * @param nodeNumber
     *          Call with next not yet given node number.
     * @return  Highest node subNode number.
     */
    uint32_t DumpConnections(uint32_t nodeNumber)
    {
        if (_size == 1)
            return nodeNumber;

        std::cout << nodeNumber << " " << nodeNumber + 1 << std::endl;
        uint32_t nodeNumber1 = _child1->DumpConnections(nodeNumber + 1);
        std::cout << nodeNumber << " " << nodeNumber1 + 1 << std::endl;
        uint32_t nodeNumber2 = _child2->DumpConnections(nodeNumber1 + 1);
        return nodeNumber2;
    }

private:
  	// Copy constructor.
  TotalizerEncodeTree (const TotalizerEncodeTree&) = default;

    // Assignment operator.
  TotalizerEncodeTree& operator = (const TotalizerEncodeTree&) = default;

};


}


#endif // TOTALIZERENCODETREE_H
