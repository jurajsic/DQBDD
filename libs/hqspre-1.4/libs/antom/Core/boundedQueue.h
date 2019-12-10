/********************************************************************************************
boundedQueue.h -- Copyright (c) 2018, Sven Reimer

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

#ifndef ANTOM_BOUNDEDQUEUE_H_
#define ANTOM_BOUNDEDQUEUE_H_


namespace antom {

  // based on glucose idea of bounded queue for maintaining the sum of the last x entries
  template <class T>
  class BoundedQueue 
  {
  public:

  BoundedQueue() :
	_first(0),
	  _last(0),
	  _sumOfQueue(0),
	  _maxsize(0),
	  _queuesize(0),
	  _elements()
		{}

	void InitSize(size_t size)
	{
	  _elements.resize(size,0);
	  _first = 0;
	  _maxsize = size;
	  _queuesize = 0;
	  _last = 0;
	}

	void Push (T element)
	{
	  //std::cout << __func__ << " " << element << std::endl;
	  // full queue -> replace last added element
	  if (_queuesize == _maxsize)
		{
		  assert(_first == _last);
		  _sumOfQueue -= _elements[_last];
		  if ((++_last) == _maxsize)
			{
			  _last = 0;
			}
		}
	  else // queue not full yet
		{
		  ++_queuesize;
		}

	  _sumOfQueue += element;
	  _elements[_first] = element;
	  if ((++_first) == _maxsize)
		{
		  _first = 0;
		  _last = 0;
		}
	}

	T Top() const
	{
	  assert(_queuesize>0);
	  return _elements[_last];
	}

	void Pop()
	{
	  _sumOfQueue -= _elements[_last];
	  --_queuesize;
	  if ( (++_last) == _maxsize)
		{
		  _last = 0;
		}
	}

	uint64_t GetSum() const
	{ return _sumOfQueue; }

	uint32_t GetAvg() const
	{
	  assert(_queuesize > 0);
	  return static_cast<uint32_t>(_sumOfQueue/static_cast<uint64_t>(_queuesize));
	}

	size_t MaxSize() const
	{ return _maxsize; }

	size_t Size() const
	{ return _queuesize; }

	// BoundedQueue is valid if it is full
	bool IsValid() const
	{
	  return (_queuesize == _maxsize);
	}

	void FastClear()
	{
	  _first = 0;
	  _last = 0;
	  _queuesize = 0;
	  _sumOfQueue = 0;
	}
		
	
  private:
	uint32_t _first;
	uint32_t _last;
	uint64_t _sumOfQueue;
	size_t _maxsize;
	size_t _queuesize;
	std::vector<T> _elements;
  };	

}

#endif
