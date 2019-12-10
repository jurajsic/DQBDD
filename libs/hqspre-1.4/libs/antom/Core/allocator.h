/********************************************************************************************
allocator.h -- Copyright (c) 2016-2017, Sven Reimer

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

#ifndef ANTOM_ALLOCATOR_H_
#define ANTOM_ALLOCATOR_H_

#include <errno.h>

namespace antom {

  class OutOfMemoryException {};
  static inline void* Xrealloc(void* ptr, size_t size)
  {
    void* mem = realloc(ptr, size);
    if (mem == NULL && errno == ENOMEM) 
	  {
		throw OutOfMemoryException();
	  }
	else 
	  {
		return mem;
	  }
  }


  // Adopted from minisat
  // RegionAllocator -- a simple class for allocating memory for clauses:
  
  template<class T>
  class RegionAllocator
  {
  public:

	typedef uint32_t Ref;
	// Undefined reference for dummy references in case memory was freed
    enum { Ref_Undef = UINT32_MAX };
	// Reference size for statistics
    enum { Unit_Size = sizeof(uint32_t) };

	explicit RegionAllocator(uint32_t initialCap = 1024*1024) :
	  _memory(NULL), 
	  _size(0),
	  _capacity(0), 
	  _wasted(0)
	{
	  //std::cout << "init: " << __func__ << " with " << initialCap << std::endl;
	  SetCapacity(initialCap);
	}

	~RegionAllocator()
	{
	  if (_memory != NULL)
		{
		  ::free(_memory);
		}
	}

	uint32_t size(void) const 
	{ return _size; }
	uint32_t Wasted(void) const 
	{ return _wasted; }
	void Free(uint32_t size) 
	{ _wasted += size; }

	Ref Alloc(uint32_t size)
	{ 
	  assert(size > 0);
	  SetCapacity(_size+size);

	  uint32_t prevSize = _size;
	  _size += size;
    
	  // Overflow
	  if (_size < prevSize)
		{
		  exit(0);
		}

	  return prevSize;
	}

	void AllocRegion(uint32_t size)
	{
	  _capacity = size;
	  _wasted = 0;
	  _size = 0;
	  _memory = (T*)Xrealloc(_memory, sizeof(T)*_capacity);
	}


	T& operator[](Ref r) 
	{ assert(r < _size); return _memory[r]; }
	const T& operator[](Ref r) const 
	{ assert(r < _size); return _memory[r]; }

	T* GetAdress(Ref r) 
	{ 
	  assert(r < _size ); 
	  return &_memory[r]; 
	}
	T* GetAdress(Ref r) const
	{ 
	  assert(r < _size); 
	  return &_memory[r]; 
	}

	Ref GetReference(const T* t)  
	{
	  assert((void*)t >= (void*)&_memory[0] && (void*)t < (void*)&_memory[_size-1]);
	  return  (Ref)(t - &_memory[0]); 
	}

	void CopyTo(RegionAllocator& to) 
	{
	  if (to._memory != NULL) 
		{
		  ::free(to._memory);
		}

	  to._memory = _memory;
	  to._size = _size;
	  to._capacity = _capacity;
	  to._wasted = _wasted;
	}

	void MoveTo(RegionAllocator& to) 
	{
	  CopyTo(to);

	  _memory = NULL;
	  _capacity = 0;
	  _size = 0;
	  _wasted = 0;
	}

	void Reset(void)
	{
	  _size = 0;
	  _wasted = 0;
	}

  uint32_t GetSize(void) const
  {
	return _size;
  }

  uint64_t GetCapacity(void) const
  {
	return _capacity;
  }

  private:
  // Copy constructor.
  RegionAllocator (const RegionAllocator&) = default;

  // Assignment operator.
  RegionAllocator& operator = (const RegionAllocator&) = default;
	
	void SetCapacity(uint32_t min_cap)
	{
	  if (_capacity >= min_cap) return;
	  //std::cout << __func__ << " " << min_cap << " current: " << _capacity << std::endl;

	  uint32_t prev_cap = _capacity;
	  while (_capacity < min_cap)
		{
		  // See minisat Region allocator:
		  // NOTE: Multiply by a factor (13/8) without causing overflow, then add 2 and make the
		  // result even by clearing the least significant bit. The resulting sequence of capacities
		  // is carefully chosen to hit a maximum capacity that is close to the '2^32-1' limit when
		  // using 'uint32_t' as indices so that as much as possible of this space can be used.
		  uint64_t delta = ((_capacity >> 1) + (_capacity >> 3) + 2) & ~1;
		  _capacity += delta;

		  // Overflow
		  if (_capacity <= prev_cap)
			{
			  exit(0);
			}
		}

	  //std::cout << "capacity: " << _capacity << " alloc: " << sizeof(T)*_capacity << std::endl;
	  assert(_capacity > 0);
	  _memory = (T*)Xrealloc(_memory, sizeof(T)*_capacity);
	}

	T* _memory;
	uint32_t _size;
	uint32_t _capacity;
	uint32_t _wasted;
  };
}

#endif
