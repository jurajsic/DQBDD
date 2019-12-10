

/********************************************************************************************
lifting.cpp -- Copyright (c) 2013, Sven Reimer

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

//#include "core.hpp"

#include <deque>
#include <map>
#include <vector>

#include "antom.hpp"
#include "solverstate.hpp"
#ifndef MEMALLOCATION
	#include "preprocessor.hpp"
#endif
#include "control.hpp"

#if 0
namespace antom {

  // Performs conflict analysis for Lifting Mode
  std::vector< uint32_t > Core::analyzeLifting (uint32_t* cc )
  {
      // Increment "m_conflicts".
      ++m_statistics.totalConflicts;

      // Everything OK wrt. "cc"?
      assert(cc != NULL);

      // Update "m_incVarActivity".
      m_incVarActivity *= m_decayFactor;
 
      // Initialization.
      uint32_t elements(0);
      uint32_t pos(m_dsEndIndex - 1); 
	  uint32_t finalmlevel(m_decisionLevel);
      uint32_t uip(0); 
      uint32_t clit(0);
      uint32_t cpos(0);
      std::vector<uint32_t> conflictClause; 
      std::vector<unsigned char> seen;
      seen.resize(m_variables + 1, false); 

      // Perform conflict analysis according to the 1UIP scheme.
      do
		{
		  // Increment "m_basicOps".
		  ++m_basicOps; 

		  // Initialization.
		  clit = cc[cpos]; 

		  // Analyze the literals of the current conflicting clause.
		  do
			{
			  // Get the index of "clit".
			  uint32_t index(clit >> 1); 

			  // "clit" not checked so far and not assigned on decision level 0? 
			  // Checking the decision level on which a particular variable has been 
			  // assigned, requires that assumptions are stored on decision levels 
			  // greater 0. Otherwise we run into problems in the incremental mode. 
			  if (!seen[index] && m_level[index] > 0)
				{
				  // Update "seen".
				  seen[index] = true;

				  // Increase the activity of the variable corresponding to "clit".
				  increaseActivity(index);

				  // Do we have a forced lit?
				  if( m_forcing[index].ptr != NULL && m_forcing[index].lit == 0 )
					{

					  //Update final decision level
					  if( finalmlevel > m_level[index] )
						{ finalmlevel = m_level[index]; }
					  ++elements; 
					}
				  else
					{ conflictClause.push_back(clit); }
				}

			  // Get the next literal to be checked.
			  clit = cc[++cpos]; 
			}
		  while (clit != 0);

		caNC:

		  // Determine the next clause to be processed. 
		  while (!seen[m_decisionStack[pos] >> 1]) { --pos; }
		  uip = m_decisionStack[pos];

		  // Update the status variables.
		  cc = m_forcing[uip >> 1].ptr;
		  assert(elements == 1 || (m_forcing[uip >> 1].lit & 1) == 1 || cc != NULL);
		  --pos;
		  --elements;

		  // Next clause a binary one?
		  if (elements > 0 && (m_forcing[uip >> 1].lit & 1) == 1)
			{
			  // Initialization.
			  uint32_t index(m_forcing[uip >> 1].lit >> 2); 
			  assert(m_level[index] == m_decisionLevel);

			  // Variable corresponding to "index" not checked so far?
			  if (!seen[index])
				{ seen[index] = true; increaseActivity(index); ++elements; }
	      
			  // Get the next clause to be processed.
			  goto caNC;
			}
	  
		  // Further update the status variables.
		  cpos = 0;
		}
      while (elements > 0);

      // Flip the sign of the UIP.
      uip = uip ^ 1;

      // Initialization.
      uint32_t size(conflictClause.size());

      // Perform a simple conflict clause minimization step. See also "Towards Understanding 
      // and Harnessing the Potential of Clause Learning" by Beame, Kautz, and Sabharwal.
      for (uint32_t l = 0; l < size; ++l)
		{
		  // Get the next literal of "conflictClause".
		  uint32_t lit(conflictClause[l]); 

		  // Do we talk about an implication?
		  if ((m_forcing[lit >> 1].lit & 1) == 1 || m_forcing[lit >> 1].ptr != NULL)
			{
			  if ((m_forcing[lit >> 1].lit & 1) == 1)
				{
				  // Is it safe to remove "lit"?
				  if (seen[m_forcing[lit >> 1].lit >> 2])
					{ --size; conflictClause[l] = conflictClause[size]; --l; }
				}
			  else
				{ 
				  // Get the forcing clause of "lit".
				  uint32_t* reason(m_forcing[lit >> 1].ptr);
				  assert(reason != NULL);

				  // Check whether all literals of "reason" have been processed during conflict analysis. 
				  // In this particular case we are allowed to remove "lit" from the conflict clause.
				  do
					{
					  // Do we have to keep "lit"?
					  if (!seen[(*reason) >> 1] && m_level[(*reason) >> 1] > 0)
						{ break; }
		  
					  // Increment "reason".
					  ++reason; 
					}
				  while (*reason != 0);
	      
				  // Is it safe to remove "lit"?
				  if (*reason == 0)
					{ --size; conflictClause[l] = conflictClause[size]; --l; }
				}
			}
		}
  
      // Add the UIP to the conflictclause 
	  conflictClause.push_back(uip);

	  return conflictClause;
  }

  // lifts a solution, given by the variables in "assumptions"
  // returns the lifted solution
  // liftingmodi:
  // 0x01 : conflict driven
  // 0x02 : brute force driven
  // 0x03 : combine 1+2
  // 0x1x : add all assumptions at once (only for conflict driven approach)
  // note: needs a negated property
  // the formula with negated property must be UNSAT, with non-negated property SAT
  std::vector< uint32_t > Core::solveLifting( std::vector<uint32_t>& assumptions, uint32_t liftingmode, uint32_t sortmode ) 
  {
    std::vector< uint32_t > liftedvars;
	std::vector< uint32_t > cclause;

    // The lifting process has to be started on decision level 0.
    assert(m_decisionLevel == 0);

	std::vector< uint32_t > assumptionvector;

    if( liftingmode & 1 )
      {
#ifndef NDEBUG
		bool receivedconflict = false;
#endif

		std::vector< std::pair< uint32_t, uint32_t > > actvector;
		
		for(uint32_t a = 0; a != assumptions.size(); ++a )
		  {
			uint32_t var = assumptions[a];
			assumptionvector.push_back(var);
			uint32_t activity = m_watches[var^1].size();
		
			actvector.push_back( std::pair< uint32_t, uint32_t>(var, activity) );
		  }

		if( sortmode & 1 )
		  { std::sort( actvector.begin(), actvector.end(), sortPairBySecond<uint32_t, uint32_t>); }
		else if ( sortmode & 2 )
		  { std::sort( actvector.begin(), actvector.end(), sortPairBySecondReverse<uint32_t, uint32_t>); }

		if ( liftingmode & 8 )
		  {
			// Set all assignments specified by "assumptions".
			for (uint32_t c = 0 ; c != actvector.size(); ++c )
			  {
				// Get the next assignment.
				uint32_t lit(actvector[c].first);

				// variable has to be initialized before (TODO)
				assert( (lit >> 1) <= m_variables);

				// Push "lit" as a "fake decision" onto the decision stack.
				if (!m_assignment[lit])
				  {
					addDecision(lit); 
				  }
			  }

			uint32_t* conclause = deduce();
			assert( conclause != 0 );
#ifndef NDEBUG
			receivedconflict = true;
#endif
			cclause = analyzeLifting( conclause );
			// Backtrack to decision level 0.
			assert(m_decisionLevel > 0 );
			backtrack(0); 
			
			for( uint32_t p = 0; p != cclause.size(); ++p )
			  { liftedvars.push_back(cclause[p]^1); }
		  }
		else
		  {
			// Set all assignments specified by "assumptions".
			for (uint32_t c = 0 ; c != actvector.size(); ++c )
			  {
				// Get the next assignment.
				uint32_t lit(actvector[c].first);

				// variable has to be initialized before (TODO)
				assert( (lit >> 1) <= m_variables);

		
				// "lit" already incorrectly assigned?
				if (m_assignment[lit ^ 1])
				  {
					if( m_level[lit>>1] != 0 )
					  {
						uint32_t* forcingclause = m_forcing[lit>>1].ptr;
						if( forcingclause == NULL )
						  {
							/*TODO*/
						  }
						cclause = analyzeLifting( forcingclause );
					  }
					else
					  { 
						liftedvars.push_back(lit); 
					  }

					// Backtrack to decision level 0.
					if (m_decisionLevel > 0)
					  { backtrack(0); }
 
					for( uint32_t p = 0; p != cclause.size(); ++p )
					  { liftedvars.push_back(cclause[p]^1); }
#ifndef NDEBUG
					receivedconflict = true;
#endif
					break;
				  }
				// "lit" currently unassigned?
				else if (!m_assignment[lit])
				  { 
					// Push "lit" as a "fake decision" onto the decision stack.
					addDecision(lit); 
	    
					// What about the effects of this implication?
					uint32_t* conclause( NULL );
					if ( ( conclause = deduce() ) != NULL)
					  {
#ifndef NDEBUG
						receivedconflict = true;
#endif
						cclause = analyzeLifting( conclause );
						// Backtrack to decision level 0.
						assert(m_decisionLevel > 0 );
						backtrack(0); 

						for( uint32_t p = 0; p != cclause.size(); ++p )
						  { liftedvars.push_back(cclause[p]^1); }

						break;
					  }
				  }
			  }
		  }

#ifndef NDEBUG
		assert( receivedconflict );
#endif
		if( liftingmode & 3 )
		  { return liftedvars; }
		
      }

	std::vector< uint32_t > vars;
	std::vector< std::pair< uint32_t, uint32_t > > actvector;

	if( liftingmode == 2 )
	  { vars = assumptions; }
	else if( liftingmode == 3 )
	  { 
		vars = liftedvars; 
		liftedvars.clear();
	  }

	for(uint32_t a = 0; a != vars.size(); ++a )
	  {
		uint32_t var = vars[a];

		uint32_t activity = m_watches[var^1].size();

		actvector.push_back( std::pair< uint32_t, uint32_t>(var, activity) );
	  }

	if( sortmode & 4 )
	  { std::sort( actvector.begin(), actvector.end(), sortPairBySecond<uint32_t, uint32_t>); }
	else if ( sortmode & 8 )
	  { std::sort( actvector.begin(), actvector.end(), sortPairBySecondReverse<uint32_t, uint32_t>); }

	// Set all assignments specified by "assumptions".
	for( uint32_t a = 0; a != actvector.size(); ++a)
	  {
		// Get the next assignment.
		uint32_t lit = actvector[a].first;

		std::vector< uint32_t > currentassumptions;

		assert( (lit >> 1) <= m_variables);

		if( lit == 0 )
		  { continue; }

		for( uint32_t b = 0; b != actvector.size(); ++b )
		  {
			if( b == a )
			  {
				// Flip assumption
				currentassumptions.push_back(lit^1);
			  }
			else if( actvector[b].first != 0 )
			  {
				currentassumptions.push_back(actvector[b].first);
			  }
		  }
				
		uint32_t sat = deduceAssumptions( currentassumptions );

		// Result is unsatisfiable? 
		// Then we can lift this current flipped literal
		if( sat == 20 )
		  {
			actvector[a].first = 0;
		  }
	  }

	assert(liftedvars.empty());

	for( uint32_t v = 0; v != actvector.size(); ++v )
	  {
		if( actvector[v].first == 0 )
		  { continue; }

		liftedvars.push_back(actvector[v].first);
	  }
    return liftedvars;
  }
}
#endif
