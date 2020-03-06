/********************************************************************************************
softClauseNodes.h -- Copyright (c) 2016, Tobias Paxian

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

#ifndef SOFTCLAUSENODES_H
#define SOFTCLAUSENODES_H

#include <vector>
#include <iostream>
#include "sorter.h"
#include "bucket.h"

namespace antom
{
    struct SoftClause;
    struct Bucket;

    /**
     * @brief The SoftClauseNodes struct - Datastructure for managing grouping by sorting
     */
    struct SoftClauseNodes
    {
       /**
       * @brief SoftClauseNode with Softclause.
       */
      SoftClauseNodes(SoftClause* SoftClause, uint32_t base) :
	  weight(SoftClause->weight),
		_base(base),
		occursHowOftenInBucket(),
		highestBucket(0),
		inHowManyBuckets(0),
		howManyOccurrencesWithoutBuckets(0),
		hasSoftClause(true),
		hasSubBuckets(false),
		softClause(SoftClause),
		subBuckets()
      {
          fillBucketOccurences();
      };
      
      /**
       * @brief SoftClauseNode with subBucket.
       */
      SoftClauseNodes(Bucket* bucket, uint64_t weight, uint32_t base) :
	  weight(weight),
		_base(base),
		occursHowOftenInBucket(),
		highestBucket(0),
		inHowManyBuckets(0),
		howManyOccurrencesWithoutBuckets(0),
		hasSoftClause(false),
		hasSubBuckets(true),
		softClause(NULL),
		subBuckets()
      {
          subBuckets.push_back(bucket);
          fillBucketOccurences();
          //subBuckets[0]->_howOftenUsed = inHowManyBuckets;
          subBuckets[0]->_howOftenUsed = HowOftenUsed();
      };

      // Destructor
      ~SoftClauseNodes(void) {};



      uint64_t weight;
      uint32_t _base;

      std::vector<uint32_t> occursHowOftenInBucket;
      uint32_t highestBucket;
      uint32_t inHowManyBuckets;
      // The sum of all occursHowOftenInBucket entries.
      uint32_t howManyOccurrencesWithoutBuckets;


      // A SoftClauseNode contains either a SoftClause or a SubBucket.
      bool hasSoftClause;
      bool hasSubBuckets;

      // SC referring to that SoftClauseNode.
      SoftClause* softClause;

      // For now at most one Subbucket!
      std::vector< Bucket* > subBuckets;

      /**
       * @brief size returns how many SC the node contains.
       * @return #SoftClauses
       */
      uint32_t size()
      {
          // for now at max one subbucket!!
          assert(subBuckets.size() <= 1);
          if (hasSubBuckets)
              return subBuckets[0]->size();
          else
              return 1;
      }

      /**
       * @brief setWeight and sets the vars inHowManyBuckets, occursHowOftenInBucket, howManyOccurrencesWithoutBuckets, highestBucket
       * @param value - the weight
       */
      void setWeight(uint64_t value)
      {
          weight = value;
          fillBucketOccurences();
      }

      uint32_t HowOftenUsed()
      {
          uint32_t howOftenUsed(0);
          for (auto v : occursHowOftenInBucket)
              howOftenUsed += v;
          return howOftenUsed;
      }

      /**
       * @brief GetOccurrences
       * @return howManyOccurrencesWithoutBuckets * size();
       */
      uint32_t GetOccurrences()
      {
          return howManyOccurrencesWithoutBuckets * size();
      }

      /**
       * @brief SortByEntriesOccurrences
       *            to std::sort by
       *                first  - inHowManyBuckets
       *                second - GetOccurences
       * @return bool to indicate if sort order is fulfilled or not
       */
      bool static SortByEntriesOccurrences(SoftClauseNodes* s1, SoftClauseNodes* s2)
      {
          return ((s1->inHowManyBuckets > s2->inHowManyBuckets) ||
                  ((s1->inHowManyBuckets == s2->inHowManyBuckets) && (s1->GetOccurrences() > s2->GetOccurrences())));
      }

      /**
       * @brief fillBucketOccurences
       *            Recalc of all SCNode vars, necessary if new weight is set.
       */
      void fillBucketOccurences()
      {
// ONLY _base 2 WOULD BE WAY EASIER!
// weight could be an easy indicator in which bucket a SC occurs,
// but for the _bases bigger than two it is way harder to get
// directly the Buckets in which they occur.
// for _base 2 it can be done with the following bit operations.
//    occursHowOftenInBucket[position]  ==  weight & (1 << position)
//                    inHowManyBuckets  ==  __builtin_popcountll(weight)
//          if (_base == 2)
//          {
//              std::cout << "weight: " << weight << std::endl;

//              howManyOccurences = __builtin_popcountll(weight);
//              std::cout << "__builtin_popcountll(weight): " << __builtin_popcountll(weight) << std::endl;
//              inHowManyBuckets = howManyOccurences;
//              highestBucket = 63 - __builtin_clzll(weight);
//              std::cout << "64 - __builtin_clzll(weight): " << 64 - __builtin_clzll(weight) << std::endl;
//              for ( uint32_t i = 0; i <= highestBucket; i++ )
//              {
//                  occursHowOftenInBucket.push_back(((weight & (1 << i)) >= 1) ? 1 : 0);
//              }
//              std::cout << "durch" << std::endl;
//              return;
//          }
          long actualQuotient = static_cast<uint64_t>(weight);
          inHowManyBuckets = 0;
          howManyOccurrencesWithoutBuckets = 0;
          occursHowOftenInBucket.clear();

          while (actualQuotient > 0)
          {
              ldiv_t divresult = ldiv( actualQuotient, _base );
              actualQuotient = divresult.quot;
              occursHowOftenInBucket.push_back(static_cast<uint32_t>(divresult.rem));
              inHowManyBuckets += divresult.rem > 0 ? 1 : 0;
              howManyOccurrencesWithoutBuckets += divresult.rem;
          }
          highestBucket = static_cast<uint32_t>(floor(log2(weight)/log2(_base)));
      }

      /**
       * @brief dumpStructure
       *            Dumps Node Information, how many entries in which bucket!
       * @param printSize
       *            true  - the bucket entries are size*occursHowOftenInBucket
       *            false - the bucket entries are occursHowOftenInBucket
       */
      void dumpStructure(bool printSize)
      {
          std::cout << std::setw(4) << inHowManyBuckets << std::setw(8) << GetOccurrences() << std::setw(15) << weight << std::setw(3) << "|";
          for (uint32_t i = 0; i <= highestBucket; ++i)
          {
              if (occursHowOftenInBucket[i] == 0)
              {
                  std::cout << std::setw(4) << "";
                  continue;
              }
              if (hasSubBuckets && printSize)
                  std::cout << std::setw(4) << size() * occursHowOftenInBucket[i];
              else
                std::cout << std::setw(4) << occursHowOftenInBucket[i];
          }
//          if (hasSubBuckets)
//          {
//              std::cout << "  <- " << subBuckets.back()->_howOftenUsed;
//          }
          std::cout << std::endl;
      }

	private:
	  
	  // Copy constructor.
	  SoftClauseNodes (const SoftClauseNodes&) = default;
	
	  // Assignment operator.
	  SoftClauseNodes& operator = (const SoftClauseNodes&) = default;	  



    };

}

#endif // SOFTCLAUSENODES_H

