/****************************************************************************************[Dimacs.h]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/
/* koshi 20140124
#ifndef Minisat_Dimacs_h
#define Minisat_Dimacs_h
*/
#ifndef MAXSATDIMACS_H
#define MAXSATDIMACS_H

#include <stdio.h>

/* koshi 20140210
#include "utils/ParseUtils.h"
*/
/* koshi 20170630
#include "maxsat0.2f-glucose3.0/ParseUtils.h"
*/
#include "utils/ParseUtils.h"
#include "core/SolverTypes.h"
#include "Pacose.h"
#include "Softclause.h"

/* koshi 20140124
namespace Minisat {
*/
namespace Glucose {

//=================================================================================================
// DIMACS Parser:

// koshi 20140106
template<class B, class Solver>
static bool readClause(B& in, Solver& S, vec<Lit>& lits,
                       long long int& weight, long long int top) {
    long long int parsed_lit; // koshi 20140106
    long long int var;
    bool soft = false;
    lits.clear();
    if (top == 0)
    {
        // unweighted MaxSAT
        // koshi 20140106 (ms: all clauses are 1-weighted soft)
        soft = true;
        weight = 1;
    } else
    {
        // weighted MaxSAT
        // koshi 20140106
        parsed_lit = parseInt(in);
        if ((1 <= parsed_lit && parsed_lit < top) || top == -1)
        {
            // soft clause
            // top == -1 indicates all clauses are weighted
            soft = true;
            weight = parsed_lit;
        } else if (parsed_lit != top)
        {
            // weight of hard clause must be top
            printf("Unexpected weight %c\n", *in), exit(3);
        }
    }

    for (;;){
        parsed_lit = parseInt(in);
        if (parsed_lit == 0) break;
        var = abs(parsed_lit)-1;
        //while (var >= S.nVars()) S.newVar();
        lits.push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
    }

    return soft;
}
// koshi 20140106
template<class B, class Solver>
static void parse_DIMACS_main(B& in, Solver& S)
{
//    int out_nbvar = 0;
    long long int out_top = 0;
    int out_nbsoft = 0;
    //  vec<long long int>& weights;
    //  vec<Lit>& blockings;
    vec<Lit> lits;

    long long int vars    = 0;
    long long int clauses = 0;
    long long int cnt     = 0;

    // koshi 20140106
    int uscnt = 0;
    long long int sWeight = 0; // sum of weights of soft clauses
    long long int usWeight = 0; // sum of weights of unit soft clauses
    //  weights.clear(); blockings.clear();
    long long int top = 0;

    for (;;){
        skipWhitespace(in);
        if (*in == EOF) break;
        else if (*in == 'p') {
            if (eagerMatch(in, "p cnf")){ // koshi 20140106 (Unweighted MaxSAT)
                vars    = parseInt(in);
                clauses = parseInt(in);
                // SATRACE'06 hack
                // if (clauses > 4000000)
                //     S.eliminate(true);
                top = 0; // all clauses are 1-weighted soft
            } else if (eagerMatch(in, "wcnf")) {
                vars    = parseInt(in);
                clauses = parseInt(in);
                while((*in == 9 ||*in == 32)) ++in; // koshi 20140117 skip space and tab
                // skipWhitespace(in);
                if (*in >= '0' && *in <= '9') {
                    top = parseInt(in);
                    printf("c top = %lld\n", top);
                }
                else {
                    top = -1; // weighted Max-SAT (no hard clause)
                    printf("c weighted Max-SAT (wms)\n");
                }
            } else {
                printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
            }
//            out_nbvar   = vars;
            out_top     = top;
            while (vars > S.nVars()) S.newVar();

            printf("c |  Number of variables:    %-12lld                                     |\n", vars);
            printf("c |  Number of clauses:      %-12lld                                     |\n", clauses);
            if (top == -1)
                printf("c |  all clauses are weigthed soft                                             |\n");
            else if (top == 0)
                printf("c |  all clauses are 1-weigthed soft                                           |\n");
            else
                printf("c |  Weight of hard clauses: %-12lld                                     |\n", top);
        } else if (*in == 'c' || *in == 'p')
        {
            skipLine(in);
        } else {
            // soft or hard clause
            long long int weight = 0;
            cnt++;
            // this should be always true, then out_top can be replaced by top!
            assert(out_top == top);
            if (readClause(in, S, lits, weight, out_top)) {// soft clause
                out_nbsoft++;
                sWeight += weight;
                S.AddSoftClause(lits, weight);
                if (lits.size() == 1) {// unit soft clause
                    uscnt++;
                    usWeight += weight;
                    continue;
                }
            }
            S.addClause_(lits);
        }
    }

    // koshi 20140107
    printf("c |  Number of soft clauses: %-12d                                     |\n", out_nbsoft);
    printf("c |  Number of unit soft clauses: %-12d                                |\n", uscnt);
    // koshi 2013.05.21
    printf("c |  Maximum Cardinality (by unit): %-12lld (%-12lld)                  |\n", sWeight, usWeight);

    // koshi 20140107
    if (vars != S.nVars())
        fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of variables.\n");

    if (cnt  != clauses)
    {
        fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of clauses.\n");
    }
}

// Inserts problem into solver.
//
template<class Solver>
static void parse_DIMACS(gzFile input_stream, Solver& S) {
    StreamBuffer in(input_stream);

    parse_DIMACS_main(in, S); }

//=================================================================================================
}

#endif // MAXSATDIMACS_H
