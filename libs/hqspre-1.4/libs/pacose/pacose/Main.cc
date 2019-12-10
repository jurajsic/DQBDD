/*****************************************************************************************[Main.cc]
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

#include <errno.h>

#include <signal.h>
#include <zlib.h>

#include "utils/System.h"
//#include "maxsat0.2f-glucose3.0/ParseUtils.h" // koshi 20140210
//#include "qmaxsat2016-g3/ParseUtils.h" // uemura 20161128
#include "utils/ParseUtils.h" // koshi 20170630
#include "utils/Options.h"
//#include "maxsat0.2f-glucose3.0/Dimacs.h" // koshi 20140124
//#include "qmaxsat2016-g3/Dimacs.h" // uemura 20161128
#include "pacose/MaxSATDimacs.h" // koshi 20170630
#include "core/Solver.h"
//#include "MaxSAT.h"

/* koshi 20140124
using namespace Minisat;
*/
using namespace Glucose;

/*
  koshi 20140106
  based on minisat2-070721/maxsat0.2e
 */
// koshi 20140106
#define TOTALIZER 128

//=================================================================================================


void printStats(Solver& solver)
{
    double cpu_time = cpuTime();
    double mem_used = memUsedPeak();
    //    double mem_used = 0;
    printf("c restarts              : %" PRIu64"\n", solver.starts);
    printf("c conflicts             : %-12" PRIu64"   (%.0f /sec)\n", solver.conflicts   , solver.conflicts   /cpu_time);
//    printf("c decisions             : %-12" PRIu64"   (%4.2f %% random) (%.0f /sec)\n", solver.decisions, (float)solver.rnd_decisions*100 / (float)solver.decisions, solver.decisions   /cpu_time);
    printf("c propagations          : %-12" PRIu64"   (%.0f /sec)\n", solver.propagations, solver.propagations/cpu_time);
//    printf("c conflict literals     : %-12" PRIu64"   (%4.2f %% deleted)\n", solver.tot_literals, (solver.max_literals - solver.tot_literals)*100 / (double)solver.max_literals);
    if (mem_used != 0) printf("c Memory used           : %.2f MB\n", mem_used);
    printf("c CPU time              : %g s\n", cpu_time);
}


static Solver* solver;
// Terminate by notifying the solver and back out gracefully. This is mainly to have a test-case
// for this feature of the Solver as it may take longer than an immediate call to '_exit()'.
static void SIGINT_interrupt(int signum) { solver->interrupt(); }

// Note that '_exit()' rather than 'exit()' has to be used. The reason is that 'exit()' calls
// destructors and may cause deadlocks if a malloc/free function happens to be running (these
// functions are guarded by locks for multithreaded use).
static void SIGINT_exit(int signum) {
    printf("\n"); printf("*** INTERRUPTED ***\n");
    if (solver->verbosity > 0){
        printStats(*solver);
        printf("\n"); printf("*** INTERRUPTED ***\n"); }
    _exit(1); }


//=================================================================================================
// Main:

int main(int argc, char** argv)
{

//    Solver S;
    double initial_time = cpuTime();

    Pacose maxSolver;

    try {
        setUsageHelp("USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");

    printf("c This is      MaxGlucose 2018\n");
    printf("c based on     QMaxSAT 2017\n");
    printf("c based on     Glucose 4.1\n");
    printf("c based on     MiniSat 2.2.0\n");

#if defined(__linux__)
        fpu_control_t oldcw, newcw;
        _FPU_GETCW(oldcw); newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE; _FPU_SETCW(newcw);
        printf("c WARNING: for repeatability, setting FPU to use double precision\n");
#endif
        // Extra options:
        //
        IntOption    verb   ("MAXSAT", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 0, IntRange(0, 2));
        IntOption    cpu_lim("MAXSAT", "cpu-lim","Limit on CPU time allowed in seconds.\n", INT32_MAX, IntRange(0, INT32_MAX));
        IntOption    mem_lim("MAXSAT", "mem-lim","Limit on memory usage in megabytes.\n", INT32_MAX, IntRange(0, INT32_MAX));
    // koshi 20140106
    /*
        IntOption    card   ("MAXSAT", "card",   "Type of SAT-encodings for Cardinality Constraints", 0, IntRange(0, 3));
    // 0: [Warners 1998]
    // 1: [Bailleux & Boufkhad 2003]
    // 2: [Asin et. al 2011]
    // 3: [Ogawa et. al 2013]
    // -1: auto // koshi 20140324
    */
        StringOption cardS   ("MAXSAT", "card",   "Type of SAT-encodings for Cardinality Constraints\n           warn, bail, asin, ogaw, wmto, mrwto, and auto", "auto");

        IntOption    comp   ("MAXSAT", "comp",
                 "Variants of SAT-encodings for Cardinality Constraints\n          warn -> 0,1,2,10,11,   bail -> 0,10,11,\n          asin -> 0,10,11,    ogaw -> 0,   wmto -> 0,   mrwto -> 0",
                 0, IntRange(0, 11));
    // koshi 20150629 for evaluation
        IntOption    pmodel   ("MAXSAT", "pmodel",   "Print a MaxSAT model", 1, IntRange(0, 1));

        printf("argc: %d\n", argc);

        parseOptions(argc, argv, true);

    int card;
    printf("c card = ");
    if (strcmp(cardS, "warn") == 0) {
      printf("warn, "); card = 0;
      maxSolver._settings._encoding = warn;
    }
    if (strcmp(cardS, "bail") == 0)  {
      printf("bail, "); card = 1;
//      MaxSolver._settings._encoding = bail;
    }
    if (strcmp(cardS, "asin") == 0) {
      printf("asin, "); card = 2;
//      MaxSolver._settings._encoding = asin;
    }
    if (strcmp(cardS, "ogaw") == 0) {
      printf("ogaw, "); card = 3;
//      MaxSolver._settings._encoding = ogaw;
    }
    if (strcmp(cardS, "bailw2") == 0) {
      printf("bailw2, "); card = 6;
//      MaxSolver._settings._encoding = bailw2;
    }
    if (strcmp(cardS, "wmto") == 0) {
      printf("wmto, "); card = 10;
//      MaxSolver._settings._encoding = wmto;
    }
    if (strcmp(cardS, "mrwto") == 0) {
      printf("mrwto, "); card = 11;
//      MaxSolver._settings._encoding = mrwto;
    }
    if (strcmp(cardS, "mrwto2") == 0) {
      printf("mrwto2, "); card = 12;
//      MaxSolver._settings._encoding = mrwto2;
    }
    if (strcmp(cardS, "auto") == 0) {
      printf("auto, "); card = -1;
//      maxSolver->_settings._encoding = heuristicQMaxSAT;
    }
    printf("comp = %d, pmodel = %d, verb = %d\n",
           (int) comp,(int) pmodel,(int) verb);

        printf("argc: %d\n", argc);

        maxSolver.verbosity = verb;

        // Use signal handlers that forcibly quit until the solver will be able to respond to
        // interrupts:
        signal(SIGINT, SIGINT_exit);
        signal(SIGXCPU,SIGINT_exit);

        // Set limit on CPU-time:
        if (cpu_lim != INT32_MAX){
            rlimit rl;
            getrlimit(RLIMIT_CPU, &rl);
            if (rl.rlim_max == RLIM_INFINITY || (rlim_t)cpu_lim < rl.rlim_max){
                rl.rlim_cur = cpu_lim;
                if (setrlimit(RLIMIT_CPU, &rl) == -1)
                    printf("WARNING! Could not set resource limit: CPU-time.\n");
            } }

        // Set limit on virtual memory:
        if (mem_lim != INT32_MAX){
            rlim_t new_mem_lim = (rlim_t)mem_lim * 1024*1024;
            rlimit rl;
            getrlimit(RLIMIT_AS, &rl);
            if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max){
                rl.rlim_cur = new_mem_lim;
                if (setrlimit(RLIMIT_AS, &rl) == -1)
                    printf("WARNING! Could not set resource limit: Virtual memory.\n");
            } }

        if (argc == 1)
            printf("Reading from standard input... Use '--help' for help.\n");

        gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
        if (in == NULL)
            printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);

    printf("c benchmark file.........: %s\n", argv[1]);

        if (maxSolver.verbosity > 0){
            printf("============================[ Problem Statistics ]=============================\n");
            printf("|                                                                             |\n"); }

    // koshi 20140107
//    int nbvar  = 0; // number of original variables
    /* weight of hard clause
       0 indicates ms (unweighted MaxSAT)
         i.e. all clauses are 1-weighted soft clauses
       -1 indicates wms (weighted MaxSAT)
         i.e. all clauses are weighted soft clauses
       positive value indicates pms or wpms (partial MaxSAT)
     */
//    int nbsoft = 0; // number of soft clauses
    vec<long long int> weights;
    vec<Lit> blockings;

    parse_DIMACS(in, maxSolver);

    //	printf("top = %d\n",top);
    //        parse_DIMACS(in, S);
        gzclose(in);
        FILE* res = (argc >= 3) ? fopen(argv[2], "wb") : NULL;

        if (maxSolver.verbosity > 0){
            printf("|  Number of variables:  %12d                                         |\n", maxSolver.nVars());
            printf("|  Number of clauses:    %12d                                         |\n", maxSolver.nClauses()); }

        double parsed_time = cpuTime();
        if (maxSolver.verbosity > 0){
            printf("|  Parse time:           %12.2f s                                       |\n", parsed_time - initial_time);
            printf("|                                                                             |\n"); }

        // Change to signal-handlers that will only notify the solver and allow it to terminate
        // voluntarily:
        signal(SIGINT, SIGINT_interrupt);
        signal(SIGXCPU,SIGINT_interrupt);

        if (!maxSolver.simplify()){
            if (res != NULL) fprintf(res, "UNSAT\n"), fclose(res);
            if (maxSolver.verbosity > 0){
                printf("===============================================================================\n");
                printf("Solved by unit propagation\n");
                printStats(maxSolver);
                printf("\n"); }
            printf("UNSATISFIABLE\n");
            exit(20);
        }



    // koshi 20140107
    lbool ret = maxSolver.MaxSolve(comp, card, pmodel);

    printStats(maxSolver);

#ifdef NDEBUG
        exit(ret == l_True ? 10 : ret == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
#else
        return (ret == l_True ? 10 : ret == l_False ? 20 : 0);
#endif
    } catch (OutOfMemoryException&){
        printf("===============================================================================\n");
        printf("INDETERMINATE\n");
        exit(0);
    }
}
