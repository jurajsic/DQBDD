/********************************************************************************************
varheap.h -- Copyright (c) 2013-2017, Tobias Schubert, Sven Reimer

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

#ifndef ANTOM_VARHEAP_H_
#define ANTOM_VARHEAP_H_

// Include standard headers.
#include <cassert>
#include <vector>

namespace antom
{
  // The "VarHeap" class.
  template < template < typename > class Order, typename T>
  class VarHeap
  {

  public:

    // Constructor.
  explicit VarHeap(const Order<T>& ord) : 
	_heap(),
      _position(),
      _size(0),
      _variables(0),
	  _comp(ord),
      _status()
		{} 

	~VarHeap (void) {}

    // Updates all data structures to be able to handle "var" variables.
    void resize(uint32_t var) 
    {
      // Resize "_position".
      _position.resize(var + 1, -1);

      // Resize "_heap".
      _heap.resize(var, 0); 

      // Update "_variables".
      _variables = var;
    }

    // Clears the entire heap.
    void clear(void)
    {
      for (uint32_t m = 0; m < _size; ++m)
		{ 
		  _position[_heap[m]] = -1; 
		}
      _size = 0; 
    }

    // Saves the current status.
    void saveStatus(void)
    {
      // Store all status variables/flags.
      _status.size = _size;
      _status.variables = _variables; 

      // Resize all vectors.
      _status.position.resize(_variables + 1, -1); 
      _status.heap.resize(_variables, 0); 

      // Store all vectors. 
      for (uint32_t v = 1; v <= _variables; ++v)
		{ _status.position[v] = _position[v]; }
      for (uint32_t p = 0; p < _size; ++p)
		{ _status.heap[p] = _heap[p]; assert(_position[_heap[p]] != -1); }
    }

    // Restores the status saved before by "saveStatus()". 
    void restoreStatus(void)
    {
      // Restore all status variables/flags.
      _size = _status.size;
      _variables = _status.variables; 

      // Resize all vectors.
      _position.resize(_variables + 1, -1); 
      _heap.resize(_variables, 0); 

      // Restore all vectors. 
      for (uint32_t v = 1; v <= _variables; ++v)
		{ _position[v] = _status.position[v]; }
      for (uint32_t p = 0; p < _size; ++p)
		{ _heap[p] = _status.heap[p]; }
    }

    // Deletes the status saved before by "saveStatus()". 
    void deleteStatus(void)
    {
      std::vector<uint32_t>().swap(_status.heap);
      std::vector<uint32_t>().swap(_status.position);
    }

    // Returns whether the heap is empty or not.
    bool empty(void) const 
	{ return (_size == 0); }

    // Returns the current size of "_heap".
    uint32_t size(void) const 
	{ return _size; }

    // Returns whether "var" is an element of "_heap".
    bool inHeap(uint32_t var) const 
    {
      // "var" has to be less or equal "_variables".
      assert(var <= _variables);

      // Return whether "var" is part of the heap.
      return (_position[var] > -1); 
    }

    // Updates the position of "var" within the heap.
    void update(uint32_t var)
    {
      // "var" has to be less or equal "_variables".
      assert(var <= _variables);

      // "var" has to be an element of "_heap".
      assert(_position[var] > -1);

      // Ensure that the heap property holds. Assumes, that "var's" 
      // activity has been incremented by the SAT solving core. 
      shiftUpwards(_position[var]); 
    }

    // Inserts variable "var" into "_heap". 
    void insert(uint32_t var) 
    {
      // "var" has to be less or equal "_variables".
      assert(var <= _variables);

      // If "var" is an element of "_heap", we've got a problem.
      assert(_position[var] == -1); 

      // Update "_position".
      _position[var] = _size;  
      
      // Add "var" to "_heap".
      _heap[_size] = var;
      
      // Increment "_size".
      ++_size; 
	
      // Ensure that the heap property holds.
      shiftUpwards(_position[var]);
    }

    // Returns the root variable.
    uint32_t top(void)
    {
      // If "_heap" is empty, we've got a problem.
      assert(_size > 0); 
      
      // Get the root variable.
      uint32_t var = _heap[0]; 

      // Decrement "_size".
      --_size;
      
      // Overwrite "_heap[0]" with the last element of "_heap".
      _heap[0] = _heap[_size];

      // Update "_position".
      _position[var] = -1;
   
      // If we removed the last element from the heap, we can skip the following operations.
      if (_size > 0)
		{
		  // Update "_position".
		  _position[_heap[0]] = 0; 
	  
		  // Ensure that the heap property holds. 
		  shiftDownwards(0);
		}
	  
      // Return "var".
      return var;
    }

    // Removes "var" from the heap.
    void remove(uint32_t var)
    {
      // If "_heap" is empty, we've got a problem.
      assert(_size > 0); 
      
      // Initialization.
      int pos(_position[var]); 

      // If "var" is not part of the heap, we've got a problem.
      assert(pos > -1); 

      // Decrement "_size".
      --_size;
      
      // Overwrite "_heap[pos]" with the last element of "_heap".
      _heap[pos] = _heap[_size];

      // Update "_position".
      _position[var] = -1;
   
      // If we removed the right-most element from the heap, we can skip the following operations.
      if ((uint32_t) pos != _size)
		{
		  // Update "_position".
		  _position[_heap[pos]] = pos; 
	  
		  // Ensure that the heap property holds. 
		  shiftDownwards(pos);
		}

      // Consistency check.
      assert(_position[var] == -1);
    }

  private:

    // Returns the position of the "father" of the element stored on position "pos".
    uint32_t father(uint32_t pos) const 
	{ return ((pos - 1) >> 1); }

    // Returns the position of the left "son" of the element stored on position "pos".
    uint32_t left(uint32_t pos) const 
	{ return ((pos << 1) + 1); }

    // Returns the position of the right "son" of the element stored on position "pos".
    uint32_t right (uint32_t pos) const 
	{ return ((pos + 1) << 1); }

    // Ensures the heap property by shifting the element on position "pos" of "_heap" upwards.
    void shiftUpwards(uint32_t pos)
    {
      // Get the variable stored on position "pos".
      uint32_t var = _heap[pos]; 
      
      // Determine the correct position of "var" within "_heap".
	  while (pos > 0 && _comp(var,_heap[father(pos)]))
		{
		  _heap[pos]             = _heap[father(pos)];
		  _position[_heap[pos]] = pos;
		  pos                     = father(pos);
		}
      
      // Store "var" at position "pos".
      _heap[pos]     = var;
      _position[var] = pos;
    }

    // Ensures the heap property by shifting the element on position "pos" of "_heap" downwards.
    void shiftDownwards(uint32_t pos)
    {      
      // Get the variable stored on position "pos".
      uint32_t var = _heap[pos]; 

      // Determine the correct position of "var" within "_heap".
      while (left(pos) < _size)
		{
		  uint32_t r = right(pos);
		  uint32_t child = r < _size && _comp(_heap[r],_heap[left(pos)]) ? r : left(pos);
	  	
		  if (!_comp(_heap[child],var))
			{ break; }
	  
		  _heap[pos]             = _heap[child];
		  _position[_heap[pos]] = pos;
		  pos                     = child; 
		}
      
      // Store "var" at position "pos".
      _heap[pos]     = var;
      _position[var] = pos;
    }

	
    // The heap.
    std::vector<uint32_t> _heap;

    // The position of a particular variable within "_heap".
    std::vector<int32_t> _position;

    // The current size of "_heap".
    uint32_t _size; 

    // The maximum number of variables for which memory has been reserved.
    uint32_t _variables;

	// Ordering class 
	Order<T> _comp;

    // A struct to be able to store the current state of the variable ordering.
    struct HeapStatus
    {
	HeapStatus(void) :
	  heap(),
		position(),
		size(0),
		variables(0)
	  {}
	  
	  ~HeapStatus(void) = default;
	  
      std::vector<uint32_t> heap;
      std::vector<uint32_t> position;
      uint32_t size; 
      uint32_t variables;
    };

	HeapStatus _status;

  };
}

#endif
