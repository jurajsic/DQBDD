
#include <errno.h>

#include <signal.h>
#include <zlib.h>

#include "utils/System.h"
#include "utils/ParseUtils.h" // koshi 20170630
#include "utils/Options.h"
#include "pacose/MaxSATDimacs.h" // koshi 20170630
#include "core/Solver.h"
#include "Encodings.h"

using namespace Glucose;

/*
  koshi 20140106
  based on minisat2-070721/maxsat0.2e
 */
// koshi 20140106
#define TOTALIZER 128

//=================================================================================================


/*
  Cardinality Constraints:
  Joost P. Warners, "A linear-time transformation of linear inequalities
  into conjunctive normal form",
  Information Processing Letters 68 (1998) 63-69
 */

// koshi 2013.04.16
bool Encodings::generateEncoding()
{
    return false;
}

void Encodings::genWarnersHalf(Lit &a, Lit &b, Lit &carry, Lit &sum, int comp,
               Solver &S, vec<Lit> &lits) {
  // carry
  lits.clear();
  lits.push(~a); lits.push(~b); lits.push(carry);  S.addClause_(lits);
  // sum
  lits.clear();
  lits.push(a); lits.push(~b); lits.push(sum);  S.addClause_(lits);
  lits.clear();
  lits.push(~a); lits.push(b); lits.push(sum);  S.addClause_(lits);
  //
  if (comp == 1 || comp == 2) {
    lits.clear();
    lits.push(carry); lits.push(sum); lits.push(~a); S.addClause_(lits);
    lits.clear();
    lits.push(carry); lits.push(sum); lits.push(~b); S.addClause_(lits);
  }
  if (comp == 2) {
    lits.clear();
    lits.push(~carry); lits.push(~sum); S.addClause_(lits);
    lits.clear();
    lits.push(~carry); lits.push(sum); lits.push(a); S.addClause_(lits);
    lits.clear();
    lits.push(~carry); lits.push(sum); lits.push(b); S.addClause_(lits);
  }
  // koshi 2013.05.31
  if (comp == 10 || comp == 11) { // [Warners 1996]
    // carry
    lits.clear(); lits.push(a); lits.push(~carry); S.addClause_(lits);
    lits.clear(); lits.push(b); lits.push(~carry); S.addClause_(lits);
    // sum
    lits.clear();
    lits.push(~a); lits.push(~b); lits.push(~sum);  S.addClause_(lits);
    lits.clear();
    lits.push(a); lits.push(b); lits.push(~sum);  S.addClause_(lits);
  }
}

// koshi 2013.04.16
void Encodings::genWarnersFull(Lit &a, Lit &b, Lit &c, Lit &carry, Lit &sum, int comp,
               Solver &S, vec<Lit> &lits) {
  // carry
  lits.clear();
  lits.push(~a); lits.push(~b); lits.push(carry); S.addClause_(lits);
  lits.clear();
  lits.push(~a); lits.push(~c); lits.push(carry); S.addClause_(lits);
  lits.clear();
  lits.push(~b); lits.push(~c); lits.push(carry); S.addClause_(lits);
  // sum
  lits.clear();
  lits.push(a); lits.push(b); lits.push(~c); lits.push(sum);
  S.addClause_(lits);
  lits.clear();
  lits.push(a); lits.push(~b); lits.push(c); lits.push(sum);
  S.addClause_(lits);
  lits.clear();
  lits.push(~a); lits.push(b); lits.push(c); lits.push(sum);
  S.addClause_(lits);
  lits.clear();
  lits.push(~a); lits.push(~b); lits.push(~c); lits.push(sum);
  S.addClause_(lits);
  if (comp == 1 || comp == 2) {
    lits.clear();
    lits.push(carry); lits.push(sum); lits.push(~a); S.addClause_(lits);
    lits.clear();
    lits.push(carry); lits.push(sum); lits.push(~b); S.addClause_(lits);
    lits.clear();
    lits.push(carry); lits.push(sum); lits.push(~c); S.addClause_(lits);
  }
  if (comp == 2) {
    lits.clear();
    lits.push(~carry); lits.push(~sum); lits.push(a); S.addClause_(lits);
    lits.clear();
    lits.push(~carry); lits.push(~sum); lits.push(b); S.addClause_(lits);
    lits.clear();
    lits.push(~carry); lits.push(~sum); lits.push(c); S.addClause_(lits);
  }
  // koshi 2013.05.31
  if (comp == 10 || comp == 11) {// [Warners 1996]
    // carry
    lits.clear();
    lits.push(a); lits.push(b); lits.push(~carry); S.addClause_(lits);
    lits.clear();
    lits.push(a); lits.push(c); lits.push(~carry); S.addClause_(lits);
    lits.clear();
    lits.push(b); lits.push(c); lits.push(~carry); S.addClause_(lits);
    // sum
    lits.clear();
    lits.push(a); lits.push(b); lits.push(c); lits.push(~sum);
    S.addClause_(lits);
    lits.clear();
    lits.push(~a); lits.push(~b); lits.push(c); lits.push(~sum);
    S.addClause_(lits);
    lits.clear();
    lits.push(~a); lits.push(b); lits.push(~c); lits.push(~sum);
    S.addClause_(lits);
    lits.clear();
    lits.push(a); lits.push(~b); lits.push(~c); lits.push(~sum);
    S.addClause_(lits);
  }
}
/*
#define wbsplit(wL,wR, ws,bs, wsL,bsL, wsR,bsR) \
  wsL.clear(); bsL.clear(); wsR.clear(); bsR.clear(); \
  for(int i = 0; i < ws.size(); i++) { \
    if (wL < wR) { \
      wsL.push(ws[i]); \
      bsL.push(bs[i]); \
      wL += ws[i]; \
    } else { \
      wsR.push(ws[i]); \
      bsR.push(bs[i]); \
      wR += ws[i]; \
    } \
  }
*/
/*
#define wbsplit(half,wL,wR, ws,bs, wsL,bsL, wsR,bsR) \
  wsL.clear(); bsL.clear(); wsR.clear(); bsR.clear(); \
  int ii = 0; \
  for(; ii < ws.size()-1; ii++) { \
    if(wL >= half) break; \
    wsL.push(ws[ii]); \
    bsL.push(bs[ii]); \
    wL += ws[ii]; \
  } \
  for(; ii < ws.size(); ii++) { \
    wsR.push(ws[ii]); \
    bsR.push(bs[ii]); \
    wR += ws[ii]; \
  }
*/
#define wbsplit(half,wL,wR, ws,bs, wsL,bsL, wsR,bsR) \
  wsL.clear(); bsL.clear(); wsR.clear(); bsR.clear(); \
  int ii = 0; \
  int wsSizeHalf = ws.size()/2; \
  for(; ii < wsSizeHalf; ii++) { \
    wsL.push(ws[ii]); \
    bsL.push(bs[ii]); \
    wL += ws[ii]; \
  } \
  for(; ii < ws.size(); ii++) { \
    wsR.push(ws[ii]); \
    bsR.push(bs[ii]); \
    wR += ws[ii]; \
  }


// koshi 2013.03.25
// Parallel counter
// koshi 2013.04.16, 2013.05.23
void Encodings::genWarners(vec<long long int>& weights, vec<Lit>& blockings,
        long long int max, int k,
        int comp, Solver& S, const Lit zero,
        vec<Lit>& lits, vec<Lit>& linkingVar) {

  linkingVar.clear();
  bool dvar = (comp == 11) ? false : true;

  if (weights.size() == 1) {
    long long int weight = weights[0];
    vec<bool> pn;
    pn.clear();
    while (weight > 0) {
      if (weight%2 == 0) pn.push(false);
      else pn.push(true);
      weight /= 2;
    }
    for(int i = 0; i < pn.size(); i++) {
      if (pn[i]) linkingVar.push(blockings[0]);
      else linkingVar.push(zero);
    }
    pn.clear();
  } else if (weights.size() > 1) {
    long long int weightL = 0; long long int weightR = 0;
    vec<long long int> weightsL, weightsR;
    vec<Lit> blockingsL, blockingsR;
    /*
    weightsL.clear(); weightsR.clear();
    blockingsL.clear(); blockingsR.clear();
    for(int i = 0; i < weights.size(); i++) {
      if (weightL < weightR) {
    weightsL.push(weights[i]);
    blockingsL.push(blockings[i]);
    weightL += weights[i];
      } else {
    weightsR.push(weights[i]);
    blockingsR.push(blockings[i]);
    weightR += weights[i];
      }
    }
    */
    const long long int half = max/2;
    (void)half; // suppress warning about unused variable
//    printf("Max: %lld  Half: %lld\n", max, half);
    wbsplit(half,weightL,weightR, weights,blockings,
        weightsL,blockingsL, weightsR,blockingsR);

    vec<Lit> alpha;
    vec<Lit> beta;
    Lit sum = mkLit(S.newVar(true,dvar));
    Lit carry = mkLit(S.newVar(true,dvar));
    genWarners(weightsL, blockingsL, weightL,k, comp, S, zero, lits,alpha);
    genWarners(weightsR, blockingsR, weightR,k, comp, S, zero, lits,beta);
    weightsL.clear(); weightsR.clear();
    blockingsL.clear(); blockingsR.clear();

    bool lessthan = (alpha.size() < beta.size());
    vec<Lit> &smalls = lessthan ? alpha : beta;
    vec<Lit> &larges = lessthan ? beta : alpha;
    assert(smalls.size() <= larges.size());

    genWarnersHalf(smalls[0],larges[0], carry,sum, comp, S,lits);
    linkingVar.push(sum);

    int i = 1;
    Lit carryN;
    for(; i < smalls.size(); i++) {
      sum = mkLit(S.newVar(true,dvar));
      carryN = mkLit(S.newVar(true,dvar));
      genWarnersFull(smalls[i],larges[i],carry, carryN,sum, comp, S,lits);
      linkingVar.push(sum);
      carry = carryN;
    }
    for(; i < larges.size(); i++) {
      sum = mkLit(S.newVar(true,dvar));
      carryN = mkLit(S.newVar(true,dvar));
      genWarnersHalf(larges[i],carry, carryN,sum, comp, S,lits);
      linkingVar.push(sum);
      carry = carryN;
    }
    linkingVar.push(carry);
    alpha.clear();beta.clear();
  }
  int lsize = linkingVar.size();
  for (int i = k; i < lsize; i++) { // koshi 2013.05.27
    //    printf("shrink: k = %d, lsize = %d\n",k,lsize);
    lits.clear();
    lits.push(~linkingVar[i]);
    S.addClause_(lits);
  }
  for (int i = k; i < lsize; i++) linkingVar.shrink(1); // koshi 2013.05.27
}

// koshi 2013.06.28
void Encodings::genWarners0(vec<long long int>& weights, vec<Lit>& blockings,
         long long int max,long long int k, int comp, Solver& S,
          vec<Lit>& lits, vec<Lit>& linkingVar) {
  // koshi 20140109
  printf("c Warners' encoding for Cardinality Constraints\n");

  int logk = 1;
  while ((k >>= 1) > 0) logk++;

  printf("c Warners' encoding for Cardinality Constraints k=%lld\n", k);
  printf("c Warners' encoding for Cardinality Constraints logk=%d\n", logk);
  Lit zero = mkLit(S.newVar());
  lits.clear();
  lits.push(~zero);
  S.addClause_(lits);
  genWarners(weights,blockings, max,logk, comp, S, zero,lits,linkingVar);
}

/*
  Cardinaltiy Constraints:
  Olivier Bailleux and Yacine Boufkhad,
  "Efficient CNF Encoding of Boolean Cardinality Constraints",
  CP 2003, LNCS 2833, pp.108-122, 2003
 */
// koshi 10.01.08
// 10.01.15 argument UB is added
void Encodings::genBailleux(vec<long long int>& weights, vec<Lit>& blockings,
         long long int total,
         Lit zero, Lit one, int comp,Solver& S,
         vec<Lit>& lits, vec<Lit>& linkingVar, long long int UB) {
  assert(weights.size() == blockings.size());

  linkingVar.clear();
  bool dvar = (comp == 11) ? false : true;

  vec<Lit> linkingAlpha;
  vec<Lit> linkingBeta;

  if (blockings.size() == 1) {// koshi 20140121
    long long int weight = weights[0];
    assert(weight < UB);
    linkingVar.push(one);
    for(int i = 0; i<weight; i++) linkingVar.push(blockings[0]);
    linkingVar.push(zero);
  } else if (blockings.size() > 1) {
    long long int weightL = 0; long long int weightR = 0;
    vec<long long int> weightsL, weightsR;
    vec<Lit> blockingsL, blockingsR;
    const long long int half = total/2;
    (void)half; // suppress warning about unused variable
    wbsplit(half, weightL,weightR, weights,blockings,
        weightsL,blockingsL, weightsR,blockingsR);

    genBailleux(weightsL,blockingsL,weightL,
        zero,one, comp,S, lits, linkingAlpha, UB);
    genBailleux(weightsR,blockingsR,weightR,
        zero,one, comp,S, lits, linkingBeta, UB);

    weightsL.clear();blockingsL.clear();
    weightsR.clear();blockingsR.clear();

    linkingVar.push(one);
    for (int i = 0; i < total && i <= UB; i++)
      linkingVar.push(mkLit(S.newVar(true,dvar)));
    linkingVar.push(zero);

    for (long long int sigma = 0; sigma <= total && sigma <= UB; sigma++) {
      for (long long int alpha = 0;
       alpha < linkingAlpha.size()-1 && alpha <= UB;
       alpha++) {
    long long int beta = sigma - alpha;
    if (0 <= beta && beta < linkingBeta.size()-1 && beta <= UB) {
      lits.clear();
      lits.push(~linkingAlpha[alpha]);
      lits.push(~linkingBeta[beta]);
      lits.push(linkingVar[sigma]);
      S.addClause_(lits);
      if (comp >= 10) {
        lits.clear();
        lits.push(linkingAlpha[alpha+1]);
        lits.push(linkingBeta[beta+1]);
        lits.push(~linkingVar[sigma+1]);
        S.addClause_(lits);
      }
    }
      }
    }
  }
  linkingAlpha.clear();
  linkingBeta.clear();
}

void Encodings::genBailleux0(vec<long long int>& weights, vec<Lit>& blockings,
          long long int max, long long int k, int comp, Solver& S,
          vec<Lit>& lits, vec<Lit>& linkingVar) {
  // koshi 20140109
  printf("c Bailleux's encoding for Cardinailty Constraints k = %lld\n", k);

  Lit one = mkLit(S.newVar());
  lits.clear();
  lits.push(one);
  S.addClause_(lits);

  genBailleux(weights,blockings,max, ~one,one, comp,S, lits, linkingVar, k);
}

/*
  Cardinaltiy Constraints:
  Robert Asin, Robert Nieuwenhuis, Albert Oliveras, Enric Rodriguez-Carbonell
  "Cardinality Networks: a theoretical and empirical study",
  Constraints (2011) 16:195-221
 */
// koshi 2013.07.01
inline void Encodings::sComparator(Lit& a, Lit& b, Lit& c1, Lit& c2,
               int comp,Solver& S, vec<Lit>& lits) {
  lits.clear();
  lits.push(~a); lits.push(~b); lits.push(c2);
  S.addClause_(lits);
  lits.clear();
  lits.push(~a); lits.push(c1);
  S.addClause_(lits);
  lits.clear();
  lits.push(~b); lits.push(c1);
  S.addClause_(lits);
  if (comp >= 10) {
    lits.clear();
    lits.push(a); lits.push(b); lits.push(~c1);
    S.addClause_(lits);
    lits.clear();
    lits.push(a); lits.push(~c2);
    S.addClause_(lits);
    lits.clear();
    lits.push(b); lits.push(~c2);
    S.addClause_(lits);
  }
}

// koshi 2013.07.01
void Encodings::genSMerge(vec<Lit>& linkA, vec<Lit>& linkB,
          Lit zero, Lit one, int comp,Solver& S,
          vec<Lit>& lits, vec<Lit>& linkingVar, long long int UB) {

  /* koshi 2013.12.10
  assert(UB > 0); is violated when k <= 1
  */

  bool lessthan = (linkA.size() <= linkB.size());
  vec<Lit> &tan = lessthan ? linkA : linkB;
  vec<Lit> &tyou = lessthan ? linkB : linkA;
  assert(tan.size() <= tyou.size());

  linkingVar.clear();
  bool dvar = (comp == 11) ? false : true;

  if (tan.size() == 0)
    for(long long int i = 0; i < tyou.size(); i++) linkingVar.push(tyou[i]);
  else if (tan.size() == 1 && tyou.size() == 1) {
    Lit c1 = mkLit(S.newVar(true,dvar));
    Lit c2 = mkLit(S.newVar(true,dvar));
    linkingVar.push(c1); linkingVar.push(c2);
    sComparator(tan[0],tyou[0], c1,c2, comp,S, lits);
  } else {
    vec<Lit> oddA,oddB, evenA,evenB;
    oddA.clear(); oddB.clear(); evenA.clear(); evenB.clear();

    long long int i;
    for(i = 0; i < tan.size(); i++) {
      if (i%2 == 0) {
    evenA.push(tan[i]); evenB.push(tyou[i]);
      } else {
    oddA.push(tan[i]); oddB.push(tyou[i]);
      }
    }
    for(; i < tyou.size(); i++) {
      if (i%2 == 0) {
    evenA.push(zero); evenB.push(tyou[i]);
      } else {
    oddA.push(zero); oddB.push(tyou[i]);
      }
    }

    // koshi 2013.07.04
    long long int UBceil = UB/2 + UB%2;
    long long int UBfloor = UB/2;
    assert(UBfloor <= UBceil);
    vec<Lit> d, e;
    genSMerge(evenA,evenB, zero,one, comp,S, lits, d, UBceil);
    genSMerge(oddA,oddB, zero,one, comp,S, lits, e, UBfloor);
    oddA.clear(); oddB.clear(); evenA.clear(); evenB.clear();

    linkingVar.push(d[0]);

    assert(d.size() >= e.size());

    while (d.size() > e.size()) e.push(zero);
    for(i = 0; i < e.size()-1; i++) {
      Lit c2i = mkLit(S.newVar(true,dvar));
      Lit c2ip1 = mkLit(S.newVar(true,dvar));
      linkingVar.push(c2i); linkingVar.push(c2ip1);
      sComparator(d[i+1],e[i], c2i,c2ip1, comp,S, lits);
    }

    linkingVar.push(e[i]);

    for (long long int i = UB+1; i < linkingVar.size(); i++) {
      lits.clear();
      lits.push(~linkingVar[i]);
      S.addClause_(lits);
    }
    long long int ssize = linkingVar.size() - UB - 1;
    if (ssize > 0) linkingVar.shrink(ssize);

    d.clear(); e.clear();
  }
}

// koshi 2013.07.01
void Encodings::genKCard(vec<long long int>& weights, vec<Lit>& blockings,
          long long int total,
          Lit zero, Lit one, int comp,Solver& S,
          vec<Lit>& lits, vec<Lit>& linkingVar, long long int UB) {

  linkingVar.clear();

  if (blockings.size() == 1) {
    long long int weight = weights[0];
    assert(weight <= UB);
    // koshi 20140121
    for(int i = 0; i<weight; i++) linkingVar.push(blockings[0]);
  } else if (blockings.size() > 1) {
    vec<Lit> linkingAlpha;
    vec<Lit> linkingBeta;

    long long int weightL = 0; long long int weightR = 0;
    vec<long long int> weightsL, weightsR;
    vec<Lit> blockingsL, blockingsR;
    const long long int half = total/2;
    (void)half; // suppress warning about unused variable
    wbsplit(half,weightL,weightR, weights,blockings,
        weightsL,blockingsL, weightsR,blockingsR);

    genKCard(weightsL,blockingsL,weightL,
        zero,one, comp,S, lits, linkingAlpha, UB);
    genKCard(weightsR,blockingsR,weightR,
        zero,one, comp,S, lits, linkingBeta, UB);

    genSMerge(linkingAlpha,linkingBeta, zero,one, comp,S, lits, linkingVar, UB);

    linkingAlpha.clear();
    linkingBeta.clear();
  }
}

// koshi 2013.07.01
void Encodings::genAsin(vec<long long int>& weights, vec<Lit>& blockings,
          long long int max, long long int k, int comp, Solver& S,
          vec<Lit>& lits, vec<Lit>& linkingVar) {
  // koshi 20140109
  printf("c Asin's encoding for Cardinailty Constraints\n");

  Lit one = mkLit(S.newVar());
  lits.clear();
  lits.push(one);
  S.addClause_(lits);

  genKCard(weights,blockings,max, ~one,one, comp,S, lits, linkingVar, k);
}


/*
  Cardinaltiy Constraints:
  Toru Ogawa, YangYang Liu, Ryuzo Hasegawa, Miyuki Koshimura, Hiroshi Fujita,
  "Modulo Based CNF Encoding of Cardinality Constraints and Its Application to
   MaxSAT Solvers",
  ICTAI 2013.
 */
// koshi 2013.10.03
void Encodings::genOgawa(long long int weightX, vec<Lit>& linkingX,
          long long int weightY, vec<Lit>& linkingY,
          long long int& total, long long int divisor,
          Lit /* zero */, Lit one, int /* comp */, Solver& S,
          vec<Lit>& lits, vec<Lit>& linkingVar, long long int /* UB */) {

  total = weightX+weightY;
  if (weightX == 0)
    for(int i = 0; i < linkingY.size(); i++) linkingVar.push(linkingY[i]);
  else if (weightY == 0)
    for(int i = 0; i < linkingX.size(); i++) linkingVar.push(linkingX[i]);
  else {
    long long int upper= total/divisor;
    long long int divisor1=divisor-1;
    /*
    printf("weightX = %lld, linkingX.size() = %d ", weightX,linkingX.size());
    printf("weightY = %lld, linkingY.size() = %d\n", weightY,linkingY.size());
    printf("upper = %lld, divisor1 = %lld\n", upper,divisor1);
    */

    linkingVar.push(one);
    for (int i = 0; i < divisor1; i++) linkingVar.push(mkLit(S.newVar()));
    linkingVar.push(one);
    for (int i = 0; i < upper; i++) linkingVar.push(mkLit(S.newVar()));
    Lit carry = mkLit(S.newVar());

    // lower part
    for (int i = 0; i < divisor; i++)
      for (int j = 0; j < divisor; j++) {
    int ij = i+j;
    lits.clear();
    lits.push(~linkingX[i]);
    lits.push(~linkingY[j]);
    if (ij < divisor) {
      lits.push(linkingVar[ij]);
      lits.push(carry);
    } else if (ij == divisor) lits.push(carry);
    else if (ij > divisor) lits.push(linkingVar[ij%divisor]);
    S.addClause_(lits);
      }

    // upper part
    for (int i = divisor; i < linkingX.size(); i++)
      for (int j = divisor; j < linkingY.size(); j++) {
    int ij = i+j-divisor;
    lits.clear();
    lits.push(~linkingX[i]);
    lits.push(~linkingY[j]);
    if (ij < linkingVar.size()) lits.push(linkingVar[ij]);
    S.addClause_(lits);
    //	printf("ij = %lld, linkingVar.size() = %lld\n",ij,linkingVar.size());
    lits.clear();
    lits.push(~carry);
    lits.push(~linkingX[i]);
    lits.push(~linkingY[j]);
    if (ij+1 < linkingVar.size()) lits.push(linkingVar[ij+1]);
    S.addClause_(lits);
      }
  }
  linkingX.clear(); linkingY.clear();
}

void Encodings::genOgawa(vec<long long int>& weights, vec<Lit>& blockings,
          long long int& total, long long int divisor,
          Lit zero, Lit one, int comp,Solver& S,
          vec<Lit>& lits, vec<Lit>& linkingVar, long long int UB) {

  linkingVar.clear();

  vec<Lit> linkingAlpha;
  vec<Lit> linkingBeta;

  if (total < divisor) {
    vec<Lit> linking;
    genBailleux(weights,blockings,total,
        zero,one, comp,S, lits, linking, UB);
    total = linking.size()-2;
    for(int i = 0; i < divisor; i++)
      if (i < linking.size()) linkingVar.push(linking[i]);
      else linkingVar.push(zero);
    linkingVar.push(one);
    linking.clear();
    //    printf("total = %lld, linkngVar.size() = %d\n", total,linkingVar.size());
  } else if (blockings.size() == 1) {
    long long int weight = weights[0];
    if (weight < UB) {
      long long int upper = weight/divisor;
      long long int lower = weight%divisor;
      long long int pad = divisor-lower-1;
      linkingVar.push(one);
      for (int i = 0; i < lower; i++) linkingVar.push(blockings[0]);
      for (int i = 0; i < pad; i++) linkingVar.push(zero);
      linkingVar.push(one);
      for (int i = 0; i < upper; i++) linkingVar.push(blockings[0]);
      total = weight;
    } else {
      lits.clear();
      lits.push(~blockings[0]);
      S.addClause_(lits);
      total = 0;
    }
  } else if (blockings.size() > 1) {
    long long int weightL = 0; long long int weightR = 0;
    vec<long long int> weightsL, weightsR;
    vec<Lit> blockingsL, blockingsR;
    const long long int half = total/2;
    (void)half; // suppress warning about unused variable
    wbsplit(half, weightL,weightR, weights,blockings,
        weightsL,blockingsL, weightsR,blockingsR);

    genOgawa(weightsL,blockingsL,weightL,divisor,
         zero,one, comp,S, lits, linkingAlpha, UB);
    genOgawa(weightsR,blockingsR,weightR,divisor,
         zero,one, comp,S, lits, linkingBeta, UB);

    weightsL.clear();blockingsL.clear();
    weightsR.clear();blockingsR.clear();

    genOgawa(weightL,linkingAlpha, weightR,linkingBeta, total,divisor,
          zero,one, comp,S, lits, linkingVar, UB);
  }
  // koshi 2013.11.12
  long long int upper = (UB-1)/divisor;
  for (long long int i = divisor+upper+1; i < linkingVar.size(); i++) {
    lits.clear();
    lits.push(~linkingVar[i]);
    S.addClause_(lits);
  }
  while (divisor+upper+2 < linkingVar.size()) linkingVar.shrink(1);
}

void Encodings::genOgawa0(int& card, // koshi 2013.12.24
           vec<long long int>& weights, vec<Lit>& blockings,
           long long int max, long long int k,
           long long int& divisor, int comp, Solver& S,
           vec<Lit>& lits, vec<Lit>& linkingVar) {
  //  koshi 20140327 assert(max >= TOTALIZER);

  /* koshi 2013.11.11
  long long int max0 = max;
  */
  long long int k0 = k;
  long long int odd = 1;
  divisor = 0;
  /* koshi 2013.11.11
  while (max0 > 0) {
    divisor++;
    max0 -= odd;
    odd += 2;
  }
  */
  while (k0 > 0) {
    divisor++;
    k0 -= odd;
    odd += 2;
  }
  printf("c max = %lld, divisor = %lld\n", max,divisor);

  // koshi 2013.12.24
  if (divisor <= 2) {
    printf("c divisor is less than or equal to 2 ");
    printf("so we use Warner's encoding, i.e. -card=warn\n");
    card = 0;
    genWarners0(weights,blockings, max,k, comp, S, lits,linkingVar);
  } else {
    // koshi 20140109
    printf("c Ogawa's encoding for Cardinality Constraints\n");

    Lit one = mkLit(S.newVar());
    lits.clear();
    lits.push(one);
    S.addClause_(lits);
    genOgawa(weights,blockings, max,divisor,
         ~one,one, comp,S, lits, linkingVar, k);
  }
}


//TODO BailW2 K-WTO
void Encodings::genBailleuxW2(vec<long long int>& weights, vec<Lit>& blockings,long long int total,Lit zero, Lit one,
        int comp,Solver& S,vec<Lit>& lits, vec<Lit>& linkingVar,vec<long long int>& linkingW , long long int UB) {

    assert(weights.size() == blockings.size());

    linkingVar.clear();
    linkingW.clear();
    bool dvar = (comp == 11) ? false : true;

    vec<Lit> linkingAlpha;
    vec<Lit> linkingBeta;

    vec<long long int> linkingWA;
    vec<long long int> linkingWB;

    if (blockings.size() == 1) {// koshi 20140121

        //1個のとき

        long long int weight = weights[0];

        if(weight >= UB){
            printf("weight(%lld) is over %lld\n" , weight , UB);
            exit(1);
        }
        //assert(weight < UB);

        linkingVar.push(one);
        linkingW.push(0);

        linkingVar.push(blockings[0]);
        linkingW.push(weights[0]);

    } else if (blockings.size() > 1) {

        //2個以上のとき

        long long int weightL = 0; long long int weightR = 0;
        vec<long long int> weightsL, weightsR;
        vec<Lit> blockingsL, blockingsR;
        const long long int half = total/2;
        (void) half; // suppress warning about unused variable 'half'
        //weightsとblockingsを半分に分ける
        wbsplit(half , weightL , weightR , weights , blockings , weightsL , blockingsL , weightsR , blockingsR);

        //LEFT
        genBailleuxW2(weightsL,blockingsL,weightL,zero,one, comp,S, lits, linkingAlpha,linkingWA , UB);

        //RIGHT
        genBailleuxW2(weightsR,blockingsR,weightR,zero,one, comp,S, lits, linkingBeta,linkingWB, UB);

        weightsL.clear();
        blockingsL.clear();
        weightsR.clear();
        blockingsR.clear();

        long long int top = ((UB < total) ? UB : total+1);
        int *table = new int[top];

        table[0] = 1;
        for (int i = 1 ; i < top ; i++){

            table[i] = 0;

        }

        int a_size = linkingWA.size();
        int b_size = linkingWB.size();

        linkingW.clear();
        linkingVar.clear();

        linkingVar.push(one);
        linkingW.push(0);

        for(int b = 1 ; b < b_size ; ++b){

            //2015 02 07
            if(linkingWB[b] < top){

                linkingVar.push(mkLit(S.newVar(true,dvar)));	//変数生成
                linkingW.push(linkingWB[b]);

                //新しく節を生成して追加
                lits.clear();
                lits.push(~linkingBeta[b]);
                lits.push(linkingVar[linkingVar.size()-1]);
                S.addClause_(lits);

                //printf("[ %d ]" , var(linkingVar[linkingVar.size()-1]));

                table[linkingWB[b]] = linkingVar.size();//1になっていたのをlinkingVar.size()に修正　2015 01 24

            }else{

                lits.clear();
                lits.push(~linkingBeta[b]);
                S.addClause_(lits);

            }

        }


        for(int a = 1 ; a < a_size ; ++a){

            long long int wa = linkingWA[a];

            if(wa >= top){
                lits.clear();
                lits.push(~linkingAlpha[a]);
                S.addClause_(lits);
                continue;

            }

            for(long long int b = 0 ; b < b_size ; ++b){

                long long int wb = linkingWB[b];

                if(wa + wb < top){

                    if(table[wa + wb] == 0){//新しい重みの和
                        linkingVar.push(mkLit(S.newVar(true,dvar)));	////変数生成
                        linkingW.push(wa + wb);
                        table[wa+wb] = linkingVar.size();	//重み(wa+wb)がlinkingVarの何番目に対応するかを記録
                        //printf("[ %d ]" , var(linkingVar[linkingVar.size()-1]));
                    }

                    //新しく節を生成して追加
                    lits.clear();
                    lits.push(~linkingAlpha[a]);
                    lits.push(~linkingBeta[b]);
                    lits.push(linkingVar[table[wa+wb]-1]);
                    S.addClause_(lits);

                }else{
                    lits.clear();
                    lits.push(~linkingAlpha[a]);
                    lits.push(~linkingBeta[b]);
                    S.addClause_(lits);
                }

            }

        }

        delete []table;

    }


    linkingAlpha.clear();
    linkingBeta.clear();
    linkingWA.clear();
    linkingWB.clear();

}

void Encodings::genBailleuxW20(vec<long long int>& weights, vec<Lit>& blockings,
          long long int max, long long int k, int comp, Solver& S,
          vec<Lit>& lits, vec<Lit>& linkingVar , vec<long long int>& linkingWeight) {
  // hayata 2014/12/17
    //printf("\nTOを構築 =====================================================\n")

    //printf("\n[bailW]\n");

    printf("c WTO encoding for Cardinailty Constraints\n");


    Lit one = mkLit(S.newVar());
    lits.clear();
    lits.push(one);
    S.addClause_(lits);

    //printf("one = %d\n" , var(one)+1);

    genBailleuxW2(weights,blockings,max, ~one,one, comp,S, lits, linkingVar,linkingWeight, k);

}




void Encodings::genCCl(Lit a, Solver& S,vec<Lit>& lits,Var varZero) { //ogawa 2013/04/02 uemura 20161129
  lits.clear();  // lits and varZero defined as global vars
  if (var(a)==varZero) {if (sign(a)==0) return;} else lits.push(a);
  S.addClause_(lits);
}

void Encodings::genCCl(Lit a, Lit b, Solver& S,vec<Lit>& lits,Var varZero) {//ogawa 2013/04/02 uemura 20161129
  lits.clear();  // lits and varZero defined as global vars
  if (var(a)==varZero) {if (sign(a)==0) return;} else lits.push(a);
  if (var(b)==varZero) {if (sign(b)==0) return;} else lits.push(b);
  S.addClause_(lits);
}

void Encodings::genCCl(Lit a, Lit b, Lit c, Solver& S,vec<Lit>& lits,Var varZero) {//ogawa 2013/04/02 uemura 20161129
  lits.clear();  // lits and varZero defined as global vars
  if (var(a)==varZero) {if (sign(a)==0) return;} else lits.push(a);
  if (var(b)==varZero) {if (sign(b)==0) return;} else lits.push(b);
  if (var(c)==varZero) {if (sign(c)==0) return;} else lits.push(c);
  S.addClause_(lits);
}

void Encodings::genCCl1(Lit a, Lit b, Lit c, Solver& S,vec<Lit>& lits,Var varZero) {//ogawa 2013/04/02 uemura 20161129
  lits.clear();  // lits and varZero defined as global vars
  printf("fe");
  if (var(a)==varZero) {if (sign(a)==0) return;} else lits.push(a);
  if (var(b)==varZero) {if (sign(b)==0) return;} else lits.push(b);
  if (var(c)==varZero) {if (sign(c)==0) return;} else lits.push(c);
  S.addClause_(lits);
}

void Encodings::genCCl(Lit a, Lit b, Lit c, Lit d, Solver& S,vec<Lit>& lits,Var varZero) {//ogawa 2013/04/02 uemura 20161129
  lits.clear();  // lits and varZero defined as global vars
  if (var(a)==varZero) {if (sign(a)==0) return;} else lits.push(a);
  if (var(b)==varZero) {if (sign(b)==0) return;} else lits.push(b);
  if (var(c)==varZero) {if (sign(c)==0) return;} else lits.push(c);
  if (var(d)==varZero) {if (sign(d)==0) return;} else lits.push(d);

  S.addClause(lits);
}

void Encodings::genCCl(Lit a, Lit b, Lit c, Lit d, Lit e, Solver& S,vec<Lit>& lits,Var varZero) {//ogawa 2013/04/02 uemura 20161129
  lits.clear();  // lits and varZero defined as global vars
  if (var(a)==varZero) {if (sign(a)==0) return;} else lits.push(a);
  if (var(b)==varZero) {if (sign(b)==0) return;} else lits.push(b);
  if (var(c)==varZero) {if (sign(c)==0) return;} else lits.push(c);
  if (var(d)==varZero) {if (sign(d)==0) return;} else lits.push(d);
  if (var(e)==varZero) {if (sign(e)==0) return;} else lits.push(e);
  S.addClause(lits);
}

//uemura 20161129
void Encodings::genKWMTO( vec<long long int>& weights ,vec<Lit>& blockings ,vec<long long int>& weightsTable,
        int from, int to, int div,Lit zero,
        vec<Lit>& lower,vec<long long int>& lowerW,vec<Lit>& upper,vec<long long int>& upperW,
        Solver& S, long long int ub,vec<Lit>& lits,Var varZero) {

    int inputsize = to-from+1;
    lower.clear();
    lowerW.clear();
    upper.clear();
    upperW.clear();

    if(inputsize == 1){
        //1個のとき

        long long int weight = weights[from];

        int low = weight % div;
        int up = weight / div;

        lower.push(zero);
        lowerW.push(0);

        if(low > 0){
            lower.push(blockings[from]);
            lowerW.push(low);
        }

        upper.push(zero);
        upperW.push(0);

        if(up > 0){
            upper.push(blockings[from]);
            upperW.push(up);
        }



    }else{

        int middle = inputsize/2;
        vec<Lit> alphaLow;
        vec<Lit> betaLow;
        vec<long long int> WalphaLow;
        vec<long long int> WbetaLow;

        vec<Lit> alphaUp;
        vec<Lit> betaUp;
        vec<long long int> WalphaUp;
        vec<long long int> WbetaUp;

        genKWMTO(weights ,blockings,weightsTable,from, from+middle-1,div,zero, alphaLow,WalphaLow,alphaUp,WalphaUp, S, ub,lits,varZero);

        genKWMTO(weights ,blockings,weightsTable,from+middle, to, div,zero,betaLow,WbetaLow,betaUp,WbetaUp, S, ub,lits,varZero);


        long long int total = weightsTable[to] - weightsTable[from] + weights[from];

        //LOWERの処理=====================================================================================


        int *tableLOW = new int[div];

        tableLOW[0] = 1;
        for (int i = 1 ; i < div ; i++){

            tableLOW[i] = 0;

        }

        int a_size = WalphaLow.size();
        int b_size = WbetaLow.size();

        lowerW.clear();
        lower.clear();

        lower.push(zero);
        lowerW.push(0);


        Lit C = mkLit(S.newVar());

        for(int a = 0 ; a < a_size ; ++a){

            long long int wa = WalphaLow[a];
            //printf("wa = %d\n",wa);

            for(long long int b = 0 ; b < b_size ; ++b){

                long long int wb = WbetaLow[b];
                //printf("wb = %d\n",wb);

                long long int wab = (wa + wb)%div;


                if(wa + wb < div){

                    if(tableLOW[wab] == 0){//新しい重みの和
                        lower.push(mkLit(S.newVar()));
                        lowerW.push(wab);
                        tableLOW[wab] = lower.size();	//重み(wa+wb)がlinkingVarの何番目に対応するかを記録
                        //printf("lower.size = %d\n",lower.size());

                    }

                    genCCl(~alphaLow[a] , ~betaLow[b] , lower[tableLOW[wab]-1] , C , S,lits,varZero);
                    //printf("ClauseLOW[-%d(a=%d) -%d(b=%d) %d c]\n" ,alphaLow[a],a,~betaLow[b], b,lower[tableLOW[wab]-1]);//arimura

                      /*for(int i = 0 ; i < Lits.size() ; ++i){
                          printf("%s%d %s " , sign(Lits[i]) == 1 ? "-" : "" , var(Lits[i]) , i == Lits.size()-1 ? "\n" : "v");
                      }*/
                }else if(wab == 0){
                    if(a!=0||b!=0){
                    genCCl(~alphaLow[a] , ~betaLow[b] , C , S,lits,varZero);
                    //printf("LOwerwab==0\n");
                    }
                      /*for(int i = 0 ; i < Lits.size() ; ++i){
                          printf("%s%d %s " , sign(Lits[i]) == 1 ? "-" : "" , var(Lits[i]) , i == Lits.size()-1 ? "\n" : "v");
                      }*/
                }else{// wa + wb > div

                    if(tableLOW[wab] == 0){//新しい重みの和

                        lower.push(mkLit(S.newVar()));
                        lowerW.push(wab);
                        tableLOW[wab] = lower.size();	//重み(wa+wb)%divがlinkingVarの何番目に対応するかを記録
                        //printf("lower.size = %d\n",lower.size());

                    }

                    genCCl(~alphaLow[a] , ~betaLow[b] , lower[tableLOW[wab]-1] , S,lits,varZero);
                    genCCl(~alphaLow[a] , ~betaLow[b] , C , S,lits,varZero);
                    //printf("ClauseLOW[-%d(a=%d) -%d(b=%d) %d c]\n" ,alphaLow[a],a,~betaLow[b], b,lower[tableLOW[wab]-1]);//arimura
                    //printf("ClauseLOW[-%d(a=%d) -%d(b=%d) c]\n" ,alphaLow[a],a,~betaLow[b], b);//arimura
                      /*for(int i = 0 ; i < Lits.size() ; ++i){
                          printf("%s%d %s " , sign(Lits[i]) == 1 ? "-" : "" , var(Lits[i]) , i == Lits.size()-1 ? "\n" : "v");
                      }*/
                }

            }

        }

        delete []tableLOW;


        WalphaLow.clear();
        WbetaLow.clear();
        alphaLow.clear();
        betaLow.clear();

        //UPPERの処理=====================================================================================

        //long long int UBU = _min(ub , total)/div;//upperの上限値
        //long long int UBU = _min(total , total)/div + 1;//upperの上限値
        long long int UBU = total/div + 1;//upperの上限値 uemura 20161129

        int *tableUP = new int[UBU+1];

        tableUP[0] = 1;
        for (int i = 1 ; i <= UBU ; i++){

            tableUP[i] = 0;

        }

        a_size = WalphaUp.size();
        b_size = WbetaUp.size();

        upperW.clear();
        upper.clear();

        upper.push(zero);
        upperW.push(0);

        for(int a = 0 ; a < a_size ; ++a){

            long long int wa = WalphaUp[a];

            for(long long int b = 0 ; b < b_size ; ++b){

                long long int wb = WbetaUp[b];

                long long int wab = wa + wb;//キャリーなしの場合

                if(UBU < wab){//超えてる場合

                    //新しく節を生成して追加
                    genCCl(~alphaUp[a] , ~betaUp[b] , S,lits,varZero);
                     /* for(int i = 0 ; i < Lits.size() ; ++i){
                          printf("%s%d %s " , sign(Lits[i]) == 1 ? "-" : "" , var(Lits[i]) , i == Lits.size()-1 ? "\n" : "v");
                      }*/
                }else{

                    if(wab > 0){

                        if(tableUP[wab] == 0){//新しい重みの和
                            upper.push(mkLit(S.newVar()));
                            upperW.push(wab);
                            tableUP[wab] = upper.size();	//重み(wa+wb)がlinkingVarの何番目に対応するかを記録
                            //printf("[ %d ]" , var(linkingVar[linkingVar.size()-1]));
                        }

                        //新しく節を生成して追加
                        genCCl(~alphaUp[a] , ~betaUp[b] , upper[tableUP[wab]-1] , S,lits,varZero);
                          /*for(int i = 0 ; i < Lits.size() ; ++i){
                              printf("%s%d %s " , sign(Lits[i]) == 1 ? "-" : "" , var(Lits[i]) , i == Lits.size()-1 ? "\n" : "v");
                          }*/
                    }

                }

                wab = wa + wb + 1;//キャリーつきの場合

                if(UBU < wab){//超えてる場合

                        genCCl(~alphaUp[a] , ~betaUp[b] , ~C , S,lits,varZero);
                          /*for(int i = 0 ; i < Lits.size() ; ++i){
                              printf("%s%d %s " , sign(Lits[i]) == 1 ? "-" : "" , var(Lits[i]) , i == Lits.size()-1 ? "\n" : "v");
                          }*/
                }else{

                    if(tableUP[wab] == 0){//新しい重みの和
                        upper.push(mkLit(S.newVar()));
                        upperW.push(wab);
                        tableUP[wab] = upper.size();	//重み(wa+wb)がlinkingVarの何番目に対応するかを記録
                    }


                    genCCl(~alphaUp[a] , ~betaUp[b] , ~C ,upper[tableUP[wab]-1] , S,lits,varZero);
                    if(wab == UBU){//test 2015 03 19
                        genCCl(~upper[tableUP[wab]-1] , ~C , S,lits,varZero);
                        //printf("OVER CARRY\n");
                    }
                      /*for(int i = 0 ; i < Lits.size() ; ++i){
                          printf("%s%d %s " , sign(Lits[i]) == 1 ? "-" : "" , var(Lits[i]) , i == Lits.size()-1 ? "\n" : "v");
                      }*/
                }

            }

        }

        /*fprintf(stderr , "LOW(%d) ",lowerW.size());
        for(int i = 1 ; i < lowerW.size(); ++i){
            fprintf(stderr , "%lld " , lowerW[i]);

        }
        fprintf(stderr , "\nUP(%d) " , upperW.size());
        for(int i = 1 ; i < upperW.size(); ++i){
            fprintf(stderr , "%lld " , upperW[i]);
        }
        fprintf(stderr , "\n");*/
        //if(carry){

        //}
        //printf("linkingVarUP.size = %d\n" , linkingVarUP.size());

        delete []tableUP;

        WalphaUp.clear();
        WbetaUp.clear();
        alphaUp.clear();
        betaUp.clear();

        //printf("C = %d " , var(C));

    }

    /*printf("\nU = { ");
    for(int i = 0 ; i < upper.size() ; ++i)
        printf("%lld " , upperW[i]);

    printf("} L = { ");
    for(int i = 0 ; i < lower.size() ; ++i)
        printf("%lld " , lowerW[i]);

    printf("}");*/

    /*printf("U = { ");
    for(int i = 0 ; i < upper.size() ; ++i)
        printf("%d(%lld) " , var(upper[i]) , upperW[i]);

    printf("} L = { ");
    for(int i = 0 ; i < lower.size() ; ++i)
        printf("%d(%lld) " , var(lower[i]) , lowerW[i]);

    printf("}\n");*/


}


void Encodings::genKWMTO0(int& /* card */, vec<long long int>& weights, vec<Lit>& blockings ,
        long long int max, long long int k,
        vec<long long int>& divisors,
        Solver& S,vec<Lit>& lits,
        vec<vec<Lit> >& linkingVars, vec<vec<long long int> >& linkingWeights){
    printf("c WMTO encoding for Cardinailty Constraints\n");

    for(int i = 0;i<2;i++){
        linkingVars.push();
        linkingWeights.push();
    }


    divisors.push(pow(max,1.0/2.0));
    printf("c p = %lld\n",divisors[0]);

    vec<long long int> weightsTable;
    long long int tmp = 0;
    int size = blockings.size();

    for(int i = 0 ; i < size ; ++i){
        tmp += weights[i];
        weightsTable.push(tmp);
    }

    Lit zero = mkLit(S.newVar());
    lits.clear();
    lits.push(zero);
    S.addClause_(lits);

    Var varZero = var(zero);

    genKWMTO( weights ,blockings , weightsTable , 0, size-1 ,divisors[0],zero, linkingVars[0] , linkingWeights[0] ,linkingVars[1] , linkingWeights[1] , S , k,lits,varZero);

}




//MRWTO UEMURA 20161112
void Encodings::genMRWTO(vec<long long int>& weights ,vec<Lit>& blockings ,vec<long long int>& weightsTable,
        int from, int to, vec<long long int>& divisors,Lit zero,
        vec<vec<Lit> >& linkingVars,vec<vec<long long int> >& linkingWeights,
        Solver& S, long long int ub,vec<Lit>& lits,Var varZero) {

    int ndigit = linkingVars.size();

    int inputsize = to-from+1;

    for(int i = 0;i<ndigit;i++){
        linkingVars[i].clear();
        linkingWeights[i].clear();
    }
    if(inputsize == 1){
        //1個のとき

        long long int tmpw = weights[from];
        int *digit = new int[ndigit];

        for(int i = 0;i<ndigit-1;i++){
            digit[i] = tmpw % divisors[i];
            tmpw = tmpw / divisors[i];
        }
        digit[ndigit-1] = tmpw;

        for(int i=0;i<ndigit;i++){
            linkingVars[i].push(zero);
            linkingWeights[i].push(0);
            if(digit[i] > 0){
                linkingVars[i].push(blockings[from]);
                linkingWeights[i].push(digit[i]);
            }
        }

        delete []digit;
    }else{

        int middle = inputsize/2;

        vec<vec<Lit> > alphalinkingVars;
        vec<vec<Lit> > betalinkingVars;
        vec<vec<long long int> > alphalinkingWeights;
        vec<vec<long long int> > betalinkingWeights;
        for(int i = 0;i<ndigit;i++){
            alphalinkingVars.push();
            betalinkingVars.push();
            alphalinkingWeights.push();
            betalinkingWeights.push();
        }

        genMRWTO(weights ,blockings,weightsTable,from, from+middle-1,divisors,zero, alphalinkingVars,alphalinkingWeights, S, ub,lits,varZero);
        genMRWTO(weights ,blockings,weightsTable,from+middle, to, divisors,zero,betalinkingVars,betalinkingWeights, S, ub,lits,varZero);

        long long int total = weightsTable[to] - weightsTable[from] + weights[from];
        //Lit *C = new Lit[ndigit];
        vec<Lit> C;

        for(int cdigit = 0;cdigit<ndigit;cdigit++){

            //最下位桁の処理=====================================================================================
            if(cdigit == 0){
                C.push();
                int div = divisors[cdigit];
                int *tableLOW = new int[div];
                tableLOW[0] = 1;
                for (int i = 1 ; i < div ; i++){
                    tableLOW[i] = 0;
                }
                int a_size = alphalinkingWeights[cdigit].size();
                int b_size = betalinkingWeights[cdigit].size();
                linkingVars[cdigit].clear();
                linkingWeights[cdigit].clear();
                linkingVars[cdigit].push(zero);
                linkingWeights[cdigit].push(0);
                C[cdigit] = mkLit(S.newVar());

                for(int a = 0 ; a < a_size ; ++a){
                    long long int wa = alphalinkingWeights[cdigit][a];

                    for(long long int b = 0 ; b < b_size ; ++b){

                        long long int wb = betalinkingWeights[cdigit][b];
                        long long int wab = (wa + wb)%div;
                        if(wa + wb < div){
                            if(tableLOW[wab] == 0){//新しい重みの和
                                linkingVars[cdigit].push(mkLit(S.newVar()));
                                linkingWeights[cdigit].push(wab);
                                tableLOW[wab] = linkingVars[cdigit].size();	//重み(wa+wb)がlinkingVars[cdigit]の何番目に対応するかを記録
                            }
                            genCCl(~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , linkingVars[cdigit][tableLOW[wab]-1] , C[cdigit] , S,lits,varZero);

                        }else if(wab == 0){
                            if(a!=0||b!=0){
                                genCCl(~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , C[cdigit] , S,lits,varZero);
                            }

                        }else{// wa + wb > div
                            if(tableLOW[wab] == 0){//新しい重みの和
                                linkingVars[cdigit].push(mkLit(S.newVar()));
                                linkingWeights[cdigit].push(wab);
                                tableLOW[wab] = linkingVars[cdigit].size();	//重み(wa+wb)%divがlinkingVarの何番目に対応するかを記録
                            }
                            genCCl(~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , linkingVars[cdigit][tableLOW[wab]-1] , S,lits,varZero);
                            genCCl(~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , C[cdigit] , S,lits,varZero);
                        }
                    }
                }
                delete []tableLOW;
                //alphalinkingWeights[cdigit].clear();
                //betalinkingWeights[cdigit].clear();
                //alphalinkingVars[cdigit].clear();
                //betalinkingVars[cdigit].clear();
            }

            //中間桁の処理====================================================================================================
            else if(cdigit > 0 && cdigit < ndigit-1){
                C.push();
                int div = divisors[cdigit];
                int *tableMIDDLE = new int[div];
                tableMIDDLE[0] = 1;
                for (int i = 1 ; i < div ; i++){
                    tableMIDDLE[i] = 0;
                }

                int a_size = alphalinkingWeights[cdigit].size();
                int b_size = betalinkingWeights[cdigit].size();
                linkingVars[cdigit].clear();
                linkingWeights[cdigit].clear();
                linkingVars[cdigit].push(zero);
                linkingWeights[cdigit].push(0);
                C[cdigit] = mkLit(S.newVar());
                for(int a = 0 ; a < a_size ; ++a){

                    long long int wa = alphalinkingWeights[cdigit][a];

                    for(long long int b = 0 ; b < b_size ; ++b){

                        long long int wb = betalinkingWeights[cdigit][b];
                        long long int wab = (wa + wb)%div;

                        if(wa + wb < div){

                            if(tableMIDDLE[wab] == 0){//新しい重みの和
                                linkingVars[cdigit].push(mkLit(S.newVar()));
                                linkingWeights[cdigit].push(wab);
                                tableMIDDLE[wab] = linkingVars[cdigit].size();	//重み(wa+wb)がlinkingVarの何番目に対応するかを記録
                            }
                            genCCl(~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , linkingVars[cdigit][tableMIDDLE[wab]-1] , C[cdigit] , S,lits,varZero);

                        }else if(wab == 0){
                            if(a!=0||b!=0){
                                genCCl(~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , C[cdigit] , S,lits,varZero);
                            }
                        }else{// wa + wb >= div

                            if(tableMIDDLE[wab] == 0){//新しい重みの和
                                linkingVars[cdigit].push(mkLit(S.newVar()));
                                linkingWeights[cdigit].push(wab);
                                tableMIDDLE[wab] = linkingVars[cdigit].size();	//重み(wa+wb)%divがlinkingVarの何番目に対応するかを記録
                            }
                            genCCl(~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , linkingVars[cdigit][tableMIDDLE[wab]-1] , S,lits,varZero);
                            genCCl(~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , C[cdigit] , S,lits,varZero);
                        }

                        wab = (wa + wb + 1)%div; //carryがある場合
                        if(wa + wb + 1 < div){
                            if(tableMIDDLE[wab] == 0){//新しい重みの和
                                linkingVars[cdigit].push(mkLit(S.newVar()));
                                linkingWeights[cdigit].push(wab);
                                tableMIDDLE[wab] = linkingVars[cdigit].size();	//重み(wa+wb)がlinkingVarの何番目に対応するかを記録
                            }
                            genCCl(~C[cdigit-1],~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , linkingVars[cdigit][tableMIDDLE[wab]-1] , C[cdigit] , S,lits,varZero);

                        }else if(wab == 0){
                            if(a!=0||b!=0){
                                genCCl(~C[cdigit-1],~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , C[cdigit] , S,lits,varZero);
                            }

                        }else{// wa + wb + 1 > div

                            if(tableMIDDLE[wab] == 0){//新しい重みの和
                                linkingVars[cdigit].push(mkLit(S.newVar()));
                                linkingWeights[cdigit].push(wab);
                                tableMIDDLE[wab] = linkingVars[cdigit].size();	//重み(wa+wb)%divがlinkingVarの何番目に対応するかを記録
                            }
                            genCCl(~C[cdigit-1],~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , linkingVars[cdigit][tableMIDDLE[wab]-1] , S,lits,varZero);
                            genCCl(~C[cdigit-1],~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , C[cdigit] , S,lits,varZero);
                        }
                    }
                }
                delete []tableMIDDLE;
                //alphalinkingWeights[cdigit].clear();
                //betalinkingWeights[cdigit].clear();
                //alphalinkingVars[cdigit].clear();
                //betalinkingVars[cdigit].clear();
            }


            //最上位桁の処理=====================================================================================
            else if(cdigit == ndigit-1){
                //long long int UBU = _min(ub , total)/div; //linkingVars[cdigit]の上限値
                long long int UBU = total;		//linkingVars[cdigit]の上限値
                for(int i = 0;i<cdigit;i++){
                    UBU /= divisors[i];
                }
                UBU++;

                int *tableUP = new int[UBU+1];
                tableUP[0] = 1;
                for (int i = 1 ; i <= UBU ; i++){
                    tableUP[i] = 0;
                }

                int a_size = alphalinkingWeights[cdigit].size();
                int b_size = betalinkingWeights[cdigit].size();
                linkingWeights[cdigit].clear();
                linkingVars[cdigit].clear();
                linkingVars[cdigit].push(zero);
                linkingWeights[cdigit].push(0);

                for(int a = 0 ; a < a_size ; ++a){

                    long long int wa = alphalinkingWeights[cdigit][a];

                    for(long long int b = 0 ; b < b_size ; ++b){

                        long long int wb = betalinkingWeights[cdigit][b];

                        long long int wab = wa + wb;//キャリーなしの場合

                        if(UBU < wab){//超えてる場合
                    //新しく節を生成して追加
                            genCCl(~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , S,lits,varZero);
                        }else{
                            if(wab > 0){
                                if(tableUP[wab] == 0){//新しい重みの和
                                    linkingVars[cdigit].push(mkLit(S.newVar()));
                                    linkingWeights[cdigit].push(wab);
                                    tableUP[wab] = linkingVars[cdigit].size();	//重み(wa+wb)がlinkingVarの何番目に対応するかを記録
                                }

                                //新しく節を生成して追加
                                genCCl(~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , linkingVars[cdigit][tableUP[wab]-1] , S,lits,varZero);
                            }
                        }

                        wab = wa + wb + 1;//キャリーつきの場合
                        if(UBU < wab){//超えてる場合
                            genCCl(~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , ~C[cdigit-1] , S,lits,varZero);
                        }else{
                            if(tableUP[wab] == 0){//新しい重みの和
                                linkingVars[cdigit].push(mkLit(S.newVar()));
                                linkingWeights[cdigit].push(wab);
                                tableUP[wab] = linkingVars[cdigit].size();	//重み(wa+wb)がlinkingVarの何番目に対応するかを記録
                            }
                            genCCl(~alphalinkingVars[cdigit][a] , ~betalinkingVars[cdigit][b] , ~C[cdigit-1] ,linkingVars[cdigit][tableUP[wab]-1] , S,lits,varZero);
                            if(wab == UBU){//test 2015 03 19
                                genCCl(~linkingVars[cdigit][tableUP[wab]-1] , ~C[cdigit-1] , S,lits,varZero);
                            }
                        }
                    }
                }
                delete[] tableUP;
            }
        }
        C.clear();
        alphalinkingWeights.clear();
        betalinkingWeights.clear();
        alphalinkingVars.clear();
        betalinkingVars.clear();
    }
}

void Encodings::genMRWTO0(int& card, vec<long long int>& weights, vec<Lit>& blockings ,
        long long int max, long long int k,
        vec<long long int>& divisors,
        Solver& S,vec<Lit>& lits,
        vec<vec<Lit> >& linkingVars, vec<vec<long long int> >& linkingWeights){

    printf("c MRWTO encoding for Cardinailty Constraints\n");

    int ndigit = 2; //mrwtoに用いる変数 基数の数 uemura 2016.11.08
    int mrdiv = 2;
    if(card == 11) {
        for(long long int i = 512;i<max;i *=16){
            ndigit++;
        }
        printf("c number of digit = %d\n",ndigit);
    }else if(card == 12){
        ndigit = 1;
        long long int wtmp = max;
        while(wtmp>mrdiv-1){
            wtmp /=mrdiv;
            ndigit++;
        }
        printf("c number of digit = %d\n",ndigit);
    }


    //linkingVars,linkingWeightSのサイズの設定
    for(int i=0;i<ndigit;i++){
        linkingVars.push();
        linkingWeights.push();
    }

    //uemura MRWTO用のdivisorsの計算
    printf("c ");
    if (card == 11){
        for(int k = 0;k<ndigit-1;k++){
            divisors.push(static_cast<long long int>(pow(max,(1.0/ndigit))));
            printf("p%d= %lld\t" ,k, divisors[k]);
        }
    printf("\n");
    }else if (card == 12){
        printf("c ");
        for(int i = 0;i<ndigit-1;i++){
            divisors.push(mrdiv);
            printf("p%d=\t%lld\t" ,i, divisors[i]);
        }
        printf ("\n");
    }

    vec<long long int> weightsTable;
    long long int tmp = 0;
    int size = blockings.size();

    for (int i = 0 ; i < size ; ++i){
        tmp += weights[i];
        weightsTable.push(tmp);
    }

    const Lit zero = mkLit(S.newVar());
    lits.clear();
    lits.push(zero);
    S.addClause_(lits);

    const Var varZero = var(zero);

    genMRWTO(weights ,blockings , weightsTable , 0, size-1 ,divisors,
            zero, linkingVars , linkingWeights , S , k,lits,varZero);
}
