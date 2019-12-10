#include "Pacose.h"
#include "Settings.h"
#include "Softclause.h"

using namespace Glucose;

// koshi 13.04.05, 13.06.28, 13.07.01, 13.10.04
void lessthan(int card, vec<Lit>& linking,vec<long long int>& linkingWeight, long long int ok, long long int k,
              long long int divisor, // koshi 13.10.04
              vec<long long int>& cc, Solver& S, vec<Lit>& lits) {
    assert(k > 0);
    if (linking.size() == 0) {} else // koshi 20140124 20140129
        if (card == 1) {// Bailleux encoding (Totalizer)
            for (long long int i = k;
                 i < linking.size() && i < ok; i++) {
                lits.clear();
                lits.push(~linking[i]);
                S.addClause_(lits);
            }
        } else if (card == 2) {// Asin encoding
            for (long long int i = k-1;
                 i < linking.size() && i < ok; i++) {
                lits.clear();
                lits.push(~linking[i]);
                S.addClause_(lits);
            }

        }else if (card == 6) {//Weighted  Bailleux encoding (Weighted Totalizer) hayata 2014/12/17

            for (long long int i = 1 ; i < linking.size() ; i++) {
                long long int tmp_w = linkingWeight[i];
                if(tmp_w >= k && tmp_w < ok){
                    lits.clear();
                    lits.push(~linking[i]);
                    S.addClause_(lits);
                    //if(S.verbosity > 1)
                    //	printf("[%lld]\n" , tmp_w);
                }
                //printf("[-%d]" , var(lits[0])+1);
            }
        }
        else if (card == 3) {// Ogawa encoding (Modulo Totalizer)
            long long int upper = (k-1)/divisor;
            long long int lower = k%divisor;
            long long int oupper = ok/divisor;
            //    printf("upper = %lld, oupper = %lld\n", upper,oupper);
            if (upper < oupper)
                for (long long int i = divisor+upper+1; i < divisor+oupper+1; i++) {
                    if (linking.size() <= i) break;
                    else {
                        //  printf("linking i = %lld ",i);
                        lits.clear();
                        lits.push(~linking[i]);
                        S.addClause_(lits);
                    }
                }
            upper = k/divisor;
            lits.clear();
            lits.push(~linking[divisor+upper]);
            lits.push(~linking[lower]);
            //    printf("divisor+upper = %lld, lower = %lld\n",divisor+upper,lower);
            S.addClause_(lits);
        } else if (card == 0) {// Warners encoding
            vec<long long int> cls;
            cls.clear();

            k--;
            if (k%2 == 0) cls.push(1);
            k = k/2;
            int cnt = 1;
            long long int pos = 0x0002LL;
            while (k > 0) {
                if (k%2 == 0) cls.push(pos);
                //    else if (cls.size() == 0) cls.push(pos);
                else for(int i = 0; i < cls.size(); i++) cls[i] = cls[i] | pos;
                pos = pos << 1;
                k = k/2;
                cnt++;
            }
            for(int i = cnt; i < linking.size(); i++) {
                cls.push(pos);
                pos = pos << 1;
            }
            for(int i = 0; i < cls.size(); i++) {
                long long int x = cls[i];
                bool found = false;
                for(int j = 0; j < cc.size(); j++) {
                    if (x == cc[j]) {
                        found = true; break;
                    }
                }
                if (!found) {
                    cc.push(x); // koshi 2013.10.04
                    lits.clear();
                    int j = 0;
                    while (x > 0) {
                        if ((x & 0x0001L) == 0x0001L) {
                            lits.push(~linking[j]);
                        }
                        x = x >> 1;
                        j++;
                    }
                    S.addClause_(lits);
                }
            }
        }
}

//uemura 20161128
void lessthanMR(int /* card */, vec<vec<Lit> >& linkings, vec<vec<long long int> >& linkingWeights,
                long long int ok, long long int k, vec<long long int>& divisors,
                vec<long long int>& /* cc */, Solver& S, vec<Lit>& lits) {

    assert(k > 0);

    if (linkings[0].size() == 0) {} else {// uemura 20161112
        int ndigit = linkings.size();
        long long int tmp_k = k;
        long long int tmp_ok = ok;
        //okは前回のk

        vec<Lit> control;

        //各桁を計算して表示する
        int *sp_k = new int[ndigit];
        int *sp_ok = new int[ndigit];
        for(int i = 0;i<ndigit-1;i++){
            sp_k[i] = tmp_k % divisors[i];
            sp_ok[i] = tmp_ok % divisors[i];
            tmp_k = tmp_k / divisors[i];
            tmp_ok = tmp_ok / divisors[i];
        }
        sp_k[ndigit-1] = tmp_k;
        sp_ok[ndigit-1] = tmp_ok;

        printf("c k = %lld(",k);
        for(int i = ndigit-1;i>=0;i--){
            printf("%d/",sp_k[i]);
        }
        printf("\b)");
        for(int i=ndigit-2;i>=0;i--){
            printf("p%d = %lld ",i,divisors[i]);
        }
        printf("\n");

        /* 自分より下の桁が全て0である場合、自分の否定も問題に加える
         * そうでない場合は、自分より大きなものの否定を問題に加える。
         */
        for(int i=ndigit-1;i>=0;i--){
            int cnr_k = 0;
            int cnr_ok = 0;
            for(int k = 0;k<i;k++){
                if(sp_k[k]>0) {cnr_k=1;}
                if(sp_ok[k]>0){cnr_ok=1;}
            }
            if(cnr_k == 1) {sp_k[i]++;}
            if(cnr_ok == 1){sp_ok[i]++;}
        }

        for(int cdigit = ndigit-1; cdigit>=0; cdigit--){
            int checknextdigit = 0;
            int sp_k2 = -1;
            long long int tmp_max = 0;
            //現在の桁より下の桁がすべて0出ないことのチェック
            for(int i = cdigit-1;i>=0;i--){
                if(sp_k[i] > 0){
                    checknextdigit=1;
                }
            }
            //最上位桁の処理=======================================================================================================
            if(cdigit == ndigit-1){
                if (linkings[cdigit].size() == 0) {
                    fprintf(stderr , "ERROR : link size digit[%d] = 0\t@less thanMR\n",cdigit);
                    exit(1);
                } else {// uemura 20161112

                    for (long long int i = 1 ; i < linkings[cdigit].size() ; i++) {
                        if(linkingWeights[cdigit][i] >= sp_k[cdigit]){
                            if(linkingWeights[cdigit][i] <sp_ok[cdigit]){//tmp_ok=>sp_ok
                                lits.clear();
                                lits.push(~linkings[cdigit][i]);
                                S.addClause_(lits);
                            }
                        }else if(checknextdigit == 1){
                            if(tmp_max < linkingWeights[cdigit][i]){
                                tmp_max = linkingWeights[cdigit][i];
                                sp_k2 = i;
                            }
                        }
                    }
                    if(sp_k[cdigit] > 1 && sp_k2 > 0){
                        control.push(~linkings[cdigit][sp_k2]);
                    }
                }
            }

            else if (cdigit  >=0){
                //最上位より下の桁の処理================================================================================================-
                if (linkings[cdigit].size() == 0) {
                    fprintf(stderr , "ERROR : link size digit[%d] = 0\t@less thanMR\n",cdigit);
                    exit(1);
                } else if(sp_k[cdigit] > 0) {// uemura 20161112
                    for (long long int i = 1 ; i < linkings[cdigit].size() ; i++) {
                        if(linkingWeights[cdigit][i] >= sp_k[cdigit]){
                            lits.clear();
                            for(int ctr = 0 ; ctr < control.size() ; ++ctr){
                                lits.push(control[ctr]);
                            }
                            lits.push(~linkings[cdigit][i]);
                            S.addClause_(lits);
                        }
                        else if(checknextdigit == 1){
                            if(tmp_max < linkingWeights[cdigit][i]){
                                tmp_max = linkingWeights[cdigit][i];
                                sp_k2 = i;
                            }
                        }
                    }
                    if(sp_k2 > 0){
                        control.push(~linkings[cdigit][sp_k2]);
                    }else{
                        control.push(~linkings[cdigit][0]);
                    }
                }
            }
        }
        control.clear();
        delete []sp_k;
        delete []sp_ok;
    }
}

// koshi 2013.05.23
void Pacose::wbSort(vec<long long int>& weights, vec<Lit>& blockings,
                    vec<long long int>& sweights, vec<Lit>& sblockings) {
    sweights.clear(); sblockings.clear();
    /*
  for(int i = 0; i < weights.size(); i++) {
    int maxi = i;
    for(int j = i+1; j < weights.size(); j++) {
      if(weights[maxi] < weights[j]) maxi = j;
    }
    if (maxi != i) { // swap
      long long int tweight = weights[maxi];
      Lit tblocking = blockings[maxi];
      weights[maxi] = weights[i];
      blockings[maxi] = blockings[i];
      weights[i] = tweight;
      blockings[i] = tblocking;
    }
  }
  */
    for(int i = 0; i < weights.size(); i++) {
        sweights.push(weights[i]);
        sblockings.push(blockings[i]);
    }
}

// koshi 20140121
void Pacose::wbFilter(long long int UB, Solver& S,vec<Lit>& lits,
                      vec<long long int>& weights, vec<Lit>& blockings,
                      vec<long long int>& sweights, vec<Lit>& sblockings) {
    sweights.clear(); sblockings.clear();

    for(int i = 0; i < weights.size(); i++) {
        if (weights[i] < UB) {
            sweights.push(weights[i]);
            sblockings.push(blockings[i]);
        } else {
            lits.clear();
            lits.push(~blockings[i]);
            S.addClause(lits);
        }
    }
}

// koshi 2013.04.05, 2013.05.21, 2013.06.28, 2013.07.01, 2013.10.04
// koshi 20140121
void Pacose::genCardinals(int& card, int comp,
                          vec<long long int>& weights, vec<Lit>& blockings,
                          long long int max, long long int k,
                          long long int& divisor, // koshi 2013.10.04
                          vec<Lit>& lits, vec<Lit>& linkingVar,vec<long long int>& linkingWeight, //uemura 20161202
                          vec<long long int>& divisors, //uemura 20161128
                          vec<vec<Lit> >& linkingVars,vec<vec<long long int> >& linkingWeights) { //uemura 20161128
    assert(weights.size() == blockings.size());

    vec<long long int> sweights;
    vec<Lit> sblockings;
    //PAX: copys weights and blockings into sweits and sblockings.
    wbSort(weights,blockings, sweights,sblockings);

    //PAX: filters out all non possible weights smaller than o-Value.
    wbFilter(k,*this,lits, sweights,sblockings, weights,blockings);

    long long int sum = sumWeight(weights); // koshi 20140124
    printf("c Sum of weights = %lld\n",sum);
    printf("c A number of soft clauses remained = %d\n",blockings.size());

    if (card == -1) { // koshi 20140324 auto mode
        printf("c auto-mode for generating cardinality constraints\n");
        int logk = 0;
        int logsum = 0;
        for (long long int ok = k; ok > 0; ok = ok >> 1) logk++;
        for (long long int osum = sum; osum > 0; osum = osum >> 1) logsum++;
        printf("c logk = %d, logsum = %d\n",logk,logsum);
        if (logk+logsum < 15) {
            // Bailleux
            card = 1; comp = 0;
            printf("c Bailleux's encoding (comp=0)\n");
        } else if (k < 3) {// Warners
            card = 0; comp = 1;
            printf("c Warners' encoding (comp=1)\n");
        } else if (logsum < 17) {// Ogawa
            card = 3; comp = 0;
            printf("c Ogawa's encoding (comp=0)\n");
        } else {
            card = 0; comp = 1;
            printf("c Warners' encoding (comp=1)\n");
        }

    }

    if (weights.size() == 0) {linkingVar.clear();} else // koshi 20140124 20140129
        // koshi 2013.06.28
        if (card == 0) // Warners
            _encodings.genWarners0(weights,blockings, max,k, comp, *this, lits,linkingVar);
        else if (card == 1) // Bailleux
            _encodings.genBailleux0(weights,blockings, max,k, comp, *this, lits,linkingVar);
        else if (card == 2) // Asin
            _encodings.genAsin(weights,blockings, max,k, comp, *this, lits,linkingVar);
        else if (card == 3) // Ogawa
            _encodings.genOgawa0(card, // koshi 2013.12.24
                                 weights,blockings, max,k,divisor, comp, *this, lits,linkingVar);
        else if (card == 6) // BailleuxW2 k cardinal hayata 2015/02/06
            _encodings.genBailleuxW20(weights,blockings, max,k, comp, *this, lits,linkingVar , linkingWeight);

        else if (card == 10){//WMTO uemura 2016.11.29
            _encodings.genKWMTO0(card,weights,blockings,max,k,divisors,*this,lits,linkingVars,linkingWeights);
        }
        else if (card == 11 || card == 12){//MRWTO uemura 2016.11.29
            _encodings.genMRWTO0(card,
                                 weights,blockings ,max,k,divisors, *this, lits, linkingVars , linkingWeights );
        }
    sweights.clear(); sblockings.clear();
}

// koshi 20140106 based on minisat2-070721/maxsat0.2e
// Lit -> mkLit, addClause -> addClause_
// koshi 2013.05.23
long long int Pacose::sumWeight(const vec<long long int>& weights) const {
    long long int sum = 0;
    for (int i = 0; i < weights.size(); i++) sum += weights[i];
    return sum;
}

Pacose::Pacose():
    _settings(),
    _encodings(),
    _softClauses(),
    _nbOfOrigVars(0)
{
    verbEveryConflicts = INT32_MAX;
}


Pacose::~Pacose()
{
    for (int i = 0; i < _softClauses.size(); ++i)
    {
        delete _softClauses[i];
    }
}

void Pacose::AddSoftClause(vec<Lit> &clause, long long int weight)
{
    SoftClause *sc = nullptr;
    if (clause.size() == 1) {// unit soft clause
        sc = new SoftClause(~clause[0], clause, weight);
        _softClauses.push(sc);
    } else {
        Lit lit = mkLit(newVar());
        clause.push(lit);
        sc = new SoftClause(lit, clause, weight);
        _softClauses.push(sc);
        addClause_(clause);
    }
}

lbool Pacose::DeduceAssumptions(const vec<Lit>& externalAssumptions)
{
    // should be called before any solving process!
    assert(decisionLevel() == 0);
    if (!ok)
    {
        return l_False;
    }

    CRef confl = CRef_Undef;
    Lit next = lit_Undef;
    lbool status = l_Undef;
    externalAssumptions.copyTo(assumptions);

    for(; ;) {

        confl = propagate();

        if(confl != CRef_Undef) {
            break;
        } else {

            lastLearntClause = CRef_Undef;
            next = lit_Undef;
            while(decisionLevel() < assumptions.size()) {
                // Perform user provided assumption:
                Lit p = assumptions[decisionLevel()];
                // maybe to implement later on
                // if phaeton gives a not known variable
                // but guarantees that this var is later on
                // used!
                // while (var(p) > nVars()) newVar();
                assert(var(p) <= nVars());
                if(value(p) == l_True) {
                    // Dummy decision level:
                    newDecisionLevel();
                } else if(value(p) == l_False) {
                    // Do I need analyze final here?
                    analyzeFinal(~p, conflict);
                    return l_False;
                } else {
                    next = p;
                    break;
                }
            }

            if(next == lit_Undef) {
                break;
            }
        }
    }

    if(next == lit_Undef && confl == CRef_Undef) {
        // not unsat
        model.clear();
        // Extend & copy model:
        model.growTo(nVars());
        for(int i = 0; i < nVars(); i++) model[i] = value(i);
    } else if(confl != CRef_Undef) {
        status = l_False;
        // unsat
        if (conflict.size() == 0) {
            // unsat under any circumstances
            // explicitly assumptions independent unsat
            // on decision level 0
            ok = false;
        }
        assert(conflict.size() == 0 || ok);
    }

    cancelUntil(0);
    return l_False;

}

void Pacose::Reset()
{
    InstanceReset();

    // reset all possible settings which can be set by phaeton
    _settings.ResetCore();
}

void Pacose::InstanceReset()
{
    verbosity = 0;
    showModel = 0;
    K = opt_K;
    R = opt_R;
    sizeLBDQueue = opt_size_lbd_queue;
    sizeTrailQueue = opt_size_trail_queue;
    firstReduceDB = opt_first_reduce_db;
    incReduceDB = opt_chanseok_hack ? 0 : opt_inc_reduce_db;
    specialIncReduceDB = opt_chanseok_hack ? 0 : opt_spec_inc_reduce_db;
    lbLBDFrozenClause = opt_lb_lbd_frozen_clause;
    chanseokStrategy = opt_chanseok_hack;
    coLBDBound  = opt_chanseok_limit;
    lbSizeMinimizingClause = opt_lb_size_minimzing_clause;
    lbLBDMinimizingClause = opt_lb_lbd_minimzing_clause;
    var_decay = opt_var_decay;
    max_var_decay = opt_max_var_decay;
    clause_decay = opt_clause_decay;
    random_var_freq = opt_random_var_freq;
    random_seed = opt_random_seed;
    ccmin_mode = opt_ccmin_mode;
    phase_saving = opt_phase_saving;
    rnd_pol = false;
    rnd_init_act = opt_rnd_init_act;
    randomizeFirstDescent = false;
    garbage_frac = opt_garbage_frac;
    certifiedOutput = NULL;
    certifiedUNSAT = false;     // Not in the first parallel version
    vbyte = false;              // Not in copy constructor
    panicModeLastRemoved = 0;
    panicModeLastRemovedShared = 0;
    useUnaryWatched = false;
    promoteOneWatchedClause = true;
    solves = 0;
    starts = 0;
    decisions = 0;
    propagations = 0;
    conflicts = 0;
    conflictsRestarts = 0;
    curRestart = 1;
    glureduce = opt_glu_reduction;
    restart_inc = opt_restart_inc;
    luby_restart = opt_luby_restart;
    adaptStrategies = opt_adapt;
    luby_restart_factor = opt_luby_restart_factor;
    randomize_on_restarts = opt_randomize_phase_on_restarts;
    fixed_randomize_on_restarts = opt_fixed_randomize_phase_on_restarts;
    newDescent = 0;
    randomDescentAssignments = 0;
    forceUnsatOnNewDescent = opt_forceunsat;

    ok = true;
    cla_inc = 1;
    var_inc = 1;
    watches.clear(false);
    watchesBin.clear(false);
    unaryWatches.clear(false);
    qhead = 0;
    simpDB_assigns = -1;
    simpDB_props = 0;
    order_heap.clear(false);
    progress_estimate = 0;
    remove_satisfied = true;
    lastLearntClause = CRef_Undef;
    // Resource constraints:
    //
    conflict_budget = -1;
    propagation_budget = -1;
    asynch_interrupt = false;
    incremental = false;
    nbVarsInitialFormula = INT32_MAX;
    totalTime4Sat = 0.;
    totalTime4Unsat = 0.;
    nbSatCalls = 0;
    nbUnsatCalls = 0;

    ca.reset();

    // as long as we do not use the simp solver
    // extra clause field should be false.
    assert(ca.extra_clause_field == false);

    // Initialize  other variables
    MYFLAG = 0;

    // Initialize only first time. Useful for incremental solving (not in // version), useless otherwise
    // Kept here for simplicity
    sumLBD = 0;
    nbclausesbeforereduce = firstReduceDB;

    //copy constructor
    assigns.fill(l_Undef);
    vardata.fill(mkVarData(CRef_Undef, 0));
    activity.fill(rnd_init_act ? drand(random_seed) * 0.00001 : 0);
    seen.fill(0);
    permDiff.fill(0);
    polarity.fill(true);
    decision.fill(char());
    trail.clear(false);
    clauses.clear(false);
    learnts.clear(false);
    permanentLearnts.clear(false);
    unaryWatchedClauses.clear(false);
    lbdQueue.fastclear();
    trailQueue.fastclear();
    forceUNSAT.fill(0);
    stats.fill(0);

    // other stuff to be set re
    model.clear(false);
    conflict.clear(false);
    verbEveryConflicts = INT32_MAX;

    trail_lim.clear(false);
    assumptions.clear();
    lastDecisionLevel.clear();
    analyze_stack.clear();
    analyze_toclear.clear();
    add_tmp.clear();


}

lbool Pacose::MaxSolve(IntOption comp, int card, IntOption pmodel)
{
    _nbOfOrigVars = nVars();

    vec<Lit> blockings;
    vec<long long int> weights;
    for (int i = 0; i < _softClauses.size(); i++)
    {
        blockings.push(_softClauses[i]->relaxationLit);
        weights.push(_softClauses[i]->weight);
    }

    long long int answer = sumWeight(weights);

    /* koshi 20140404
    if (card == 3 && answer < TOTALIZER) {
      card = 1;
      printf("c Sum of weights is %lld, ", answer);
      printf("which is small, so we use normal totalizer, i.e. -card=bail\n");
    }
    */

    vec<Lit> lits;
    int lcnt = 0; // loop count
    vec<Lit> linkingVar;
    vec<long long int> linkingWeight; //uemura 20161202
    //bool mmodel[nbvar]; // koshi 2013.07.05
    bool *mmodel = new bool[_nbOfOrigVars]; //uemura 20161128
    long long int divisor = 1; // koshi 2013.10.04

    vec<long long int> ndivisor;//mrwto用の複数の基数を保存する変数 uemura 2016.12.05
    vec<vec<Lit> > linkingVarMR; //uemura 2016.12.05 for mrwto
    vec<vec<long long int> > linkingWeightMR;//uemura 2016.12.05 for mrwto

    vec<long long int> cc; // cardinality constraints
    cc.clear();
    vec<Lit> dummy;

    // koshi 20140701        lbool ret = S.solveLimited(dummy);
    lbool ret;

    while ((ret = solveLimited(dummy)) == l_True) {// koshi 20140107
        lcnt++;
        long long int answerNew = 0;

        for (int i = 0; i < blockings.size(); i++) {
            int varnum = var(blockings[i]);
            if (sign(blockings[i])) {
                if (model[varnum] == l_False) {
                    answerNew += weights[i];
                }
            } else {
                if (model[varnum] == l_True) {
                    answerNew += weights[i];
                }
            }
        }
        printf("o %lld\n",answerNew);
        //          printf("answer %lld\n",answer);
        if (lcnt > 1) assert(answerNew < answer);
        if (lcnt == 1) { // first model: generate cardinal constraints
            int ncls = nClauses();

            //genCardinals(card,comp, weights,blockings, answer,answerNew,divisor,
            //	 S, lits, linkingVar);
            // uemura 20161128
            genCardinals(card,comp, weights,blockings, answer,answerNew,divisor,
                         lits, linkingVar,linkingWeight,ndivisor,linkingVarMR , linkingWeightMR);
            //printf("c linkingVar.size() = %d\n",linkingVar.size());
            //uemura 20161129
            if (card < 10){
                printf("c linkingVar.size() = %d\n",linkingVar.size());
            }else if (card == 10 || card == 11 || card == 12){
                printf("c ");
                for(int i = 0;i<linkingVarMR.size();i++){
                    printf("linkingVar[%d].size = %d, ",i,linkingVarMR[i].size());
                }
                printf("\n");
            }
            printf("c Cardinality Constraints: %d variables and %d clauses\n",
                   nVars()-_nbOfOrigVars,nClauses()-ncls);
        }
        if (pmodel == 1) {
            for (int i = 0; i < _nbOfOrigVars; i++) {
                if (model[i]==l_True)
                    mmodel[i] = true;
                else mmodel[i] = false;
            }
        }
        if (answerNew > 0) {
            if (card == 10 || card == 11 || card == 12){
                lessthanMR(card, linkingVarMR,linkingWeightMR, answer,answerNew,ndivisor, cc, *this, lits);
            } else{
                if (card == 1 && lcnt == 1)
                    answer = linkingVar.size();
                lessthan(card, linkingVar,linkingWeight, answer,answerNew,divisor, cc, *this, lits);
            }
            answer = answerNew;
        } else {
            answer = answerNew;
            ret = l_False; // koshi 20140124
            break;
        }
    } // end of while

    // koshi 20140124
    if (ret == l_False) {
        printf((lcnt > 0) ? "s OPTIMUM FOUND\n" : "s UNSATISFIABLE\n");
    } else {
        printf("s UNKNOWN\n");
        printf("c Search is stopped by a limit (maybe time-limit)\n");
    }
    if (lcnt > 0 ) {
        if (pmodel == 1) {
            printf("v ");
            for (int i = 0; i< _nbOfOrigVars; i++) {
                printf("%s%d ", mmodel[i]?"":"-", i+1);
                if ((i+1)%20 == 0 && i+1 < _nbOfOrigVars) printf("\nv ");
            }
            printf("\n");
        }
    }

    if (lcnt > 0)
        printf("c Latest Answer = %lld by %d loops\n",answer,lcnt);

    delete[] mmodel;
    return ret;
}
