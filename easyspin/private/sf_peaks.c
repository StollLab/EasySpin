/*
sf_peaks  Bin the peaks of a coherence transfer pathway into a buffer

  sf_peaks(buffer,dt,IncScheme,idxFreeL,idxFreeR,Ea,Eb,G,D)
  sf_peaks(...,M1L,...,M(K-1)L,M1R,...,M(K-1)R)

  Input:
    buffer      complex frequency-domain buffer, modified in place: vector
                (1 dimension) or nPoints1 x nPoints2 matrix (2 dimensions)
    dt          time increment of each dimension (microseconds)
    IncScheme   incrementation scheme, K x nDimensions: interval k gets
                longer by IncScheme(k,d)*dt(d) per step along dimension d;
                one row per incremented interval (K <= 4), in sequence
                order. E.g. [1] 3p ESEEM, [1;1] 2p ESEEM, [1 0;0 1]
                HYSCORE, [1 0;0 1;0 1;1 0] 2D-CP, [1;2;1] Carr-Purcell
    idxFreeL, idxFreeR   left and right manifold (1 = alpha, 2 = beta) of
                each incremented interval
    Ea, Eb      energies of the alpha and beta manifolds (MHz)
    G           state before the first incremented interval
    D           detection operator
    MkL, MkR    left and right transfer matrices between intervals k and
                k+1 (G, D, MkL, MkR: nStates x nStates)

  Physics: during interval k, the density matrix element rho(a,b), with
  left state a in manifold idxFreeL(k) and right state b in manifold
  idxFreeR(k), oscillates with the frequency EL(a)-ER(b). Between
  intervals, rho -> MkL*rho*MkR, and the signal is trace(D*rho). Expanded,
  the signal is a sum of peaks, one for each choice of states
  a1,b1,...,aK,bK, with the amplitude

    D(bK,aK) * M(K-1)L(aK,a(K-1)) * ... * M1L(a2,a1) * G(a1,b1)
             * M1R(b1,b2) * ... * M(K-1)R(b(K-1),bK)

  and the frequency nu_d = sum_k IncScheme(k,d)*(EL_k(ak)-ER_k(bk)) along
  dimension d. There are nStates^(2K) peaks; for K > 4, saffron computes
  the signal point by point instead (Opt.Binning).

  Binning: a peak goes into bin fmod(-nu*dt,1)*nPoints, rounded to the
  nearest integer with ties away from zero, negative bins wrapped around.
  The - sign: quantum evolution goes with exp(-i*omega*t), Fourier theory
  with exp(+i*omega*t). The rounding is odd in nu, so conjugate peak pairs
  (nu,a), (-nu,conj(a)) give exactly conjugate-symmetric buffers.

  Algorithm: nested loops over the left and right state of each interval.
  Each matrix element of the chain is multiplied in within the inner one
  of the two loops over its states. Intervals that change dimension 2 are
  looped outside, so that the innermost loops add to one buffer column.
  Shortcuts: peaks with zero amplitude are skipped early (pulses that keep
  the electron in its manifold give diagonal transfer matrices). For
  K >= 2, where the outer loops reuse them, the two matrix elements with
  the innermost right state are tabulated as one product, and the
  dimension-1 bins are tabulated if only the innermost interval, or (for
  K >= 3) the two innermost intervals, change dimension 1. For K = 1, each
  table entry would be used once, so nothing is tabulated.
 */

#include "mex.h"
#include <math.h>

#define MAX_INTERVALS 4  /* nStates^(2K) peaks: binning limited to K <= 4 */

/* Each interval has a left state a and a right state b. The current
   states are kept in one array, state[]: state[LEFT(k)] is the left state
   of interval k, state[RIGHT(k)] its right state. */
#define LEFT(k)  (2*(k))
#define RIGHT(k) (2*(k)+1)

typedef mxComplexDouble Complex;

static Complex cmul(Complex x, Complex y)
{
  Complex z;
  z.real = x.real*y.real - x.imag*y.imag;
  z.imag = x.real*y.imag + x.imag*y.real;
  return z;
}

static int isZero(Complex x)
{
  return x.real==0 && x.imag==0;
}

/* adds the amplitude of a peak to its bin */
static void addPeak(Complex *bin, Complex peak)
{
  bin->real += peak.real;
  bin->imag += peak.imag;
}

/* One matrix of the amplitude chain, with the positions in state[] of the
   states that select its row and its column, e.g. G(a1,b1): row LEFT(0),
   col RIGHT(0) */
typedef struct {
  const Complex *M;
  int row, col;
} ChainMatrix;

typedef struct {
  int nStates, nIntervals, nDims;
  int nPoints[2];                /* buffer size along each dimension */
  double dt[2];                  /* time increment of each dimension (µs) */
  Complex *buffer;

  /* Incremented intervals, numbered in sequence order */
  const double *EL[MAX_INTERVALS];   /* energies of the left manifold (MHz) */
  const double *ER[MAX_INTERVALS];   /* energies of the right manifold (MHz) */
  double incScheme[MAX_INTERVALS][2];  /* incrementation scheme, IncScheme(k,d) */
  int loopOrder[MAX_INTERVALS];      /* intervals from outer to inner loop */

  /* Current states, and the chain matrices multiplied in within the loop
     over each state (at most two per state) */
  ChainMatrix chain[2*MAX_INTERVALS][2];
  int nChain[2*MAX_INTERVALS];
  int state[2*MAX_INTERVALS];

  /* Innermost interval: its right state b appears in two matrix elements
     of the chain. For K >= 2, their product is tabulated over b and the
     other states u, v of the two elements (u = state[rightOther[0]], v =
     state[rightOther[1]]) at rightProduct[(u*nStates+v)*nStates+b].
     Example HYSCORE: the innermost interval is t1, b1 appears in G(a1,b1)
     and M1R(b1,b2), so u = a1 and v = b2. NULL for K = 1: the two
     elements G(a1,b1), D(b1,a1) are multiplied in per peak. */
  Complex *rightProduct;
  int rightOther[2];
  const Complex *G, *D;          /* state before the first interval, detection operator */

  /* Tabulated dimension-1 bins of the innermost one or two intervals, or
     NULL (always for K = 1); with two intervals, one n x n block per state
     pair of the second-innermost interval */
  int *bins1;
  int bins1Intervals;

  int bin2PerPeak;               /* 1 if the innermost interval changes dimension 2 */
} Pathway;

/*-------------------------------------------------------------------*/
/* Bin of a peak with phase x = -nu*dt (in cycles): fmod(x,1)*nPoints,
   rounded to the nearest integer with ties away from zero, negative bins
   wrapped around; for finite x and nPoints>0, 0 <= bin < nPoints. The
   fractional part x-trunc(x) is exact and equal to fmod(x,1), but faster;
   doubles with |x| >= 2^52 are integers (fractional part 0), and the cast
   to long long needs |x| < 2^63. y-b is exact, so the rounding is exact
   and odd in x. */
static int binIndex(double x, int nPoints)
{
  double frac = (fabs(x)<4503599627370496.0) ? x - (double)(long long)x : 0.0;
  double y = frac*nPoints;
  int b = (int)y;
  if (y-b>=0.5) b++;
  else if (y-b<=-0.5) b--;
  if (b<0) b += nPoints;
  if (b>=nPoints) b -= nPoints;
  return b;
}

/* Multiplies amp by the chain matrices assigned to the loop over state[s]
   (all their states are set) */
static Complex multiplyChain(const Pathway *P, int s, Complex amp)
{
  int i, n = P->nStates;
  for (i=0;i<P->nChain[s];i++) {
    const ChainMatrix *m = &P->chain[s][i];
    amp = cmul(amp,m->M[P->state[m->row] + P->state[m->col]*n]);
  }
  return amp;
}

/*-------------------------------------------------------------------*/
/* Loops over the states of the innermost interval and adds the peaks to
   the buffer. amp is the amplitude of the outer part of the chain, nu1
   and nu2 the frequency sums of the outer intervals. */
static void innermostInterval(Pathway *P, Complex amp, double nu1, double nu2)
{
  /* Loop constants as local variables: the compiler cannot keep them in
     registers if they are read via P, since the buffer writes might change
     them */
  const int n = P->nStates, k = P->loopOrder[P->nIntervals-1];
  const double *EL = P->EL[k], *ER = P->ER[k];
  const double scheme1 = P->incScheme[k][0], scheme2 = P->incScheme[k][1];
  const double dt1 = P->dt[0], dt2 = P->dt[1];
  const int nPoints1 = P->nPoints[0], nPoints2 = P->nPoints[1];
  const int bin2PerPeak = P->bin2PerPeak;
  const int *bins1Block = P->bins1;  /* tabulated bins of this interval, n x n */
  const Complex *G = P->G, *D = P->D;
  Complex *column = P->buffer;
  int a, b;

  /* the dimension-2 frequency is complete unless this interval changes it */
  if (P->nDims==2 && !bin2PerPeak)
    column += (size_t)binIndex(-nu2*dt2,nPoints2)*nPoints1;

  /* bins tabulated over two intervals: the block for the current states
     of the second-innermost interval */
  if (P->bins1Intervals==2) {
    int outer = P->loopOrder[P->nIntervals-2];
    bins1Block += ((size_t)P->state[LEFT(outer)]*n + P->state[RIGHT(outer)])*n*n;
  }

  for (a=0;a<n;a++) {
    const Complex *products = NULL;
    const int *bins1 = NULL;
    Complex ampA;
    P->state[LEFT(k)] = a;
    ampA = multiplyChain(P,LEFT(k),amp);
    if (isZero(ampA)) continue;
    if (P->rightProduct) products = P->rightProduct +
      ((size_t)P->state[P->rightOther[0]]*n + P->state[P->rightOther[1]])*n;
    if (bins1Block) bins1 = bins1Block + (size_t)a*n;

    if (!products) {
      /* K = 1, nothing tabulated: the peak amplitude is G(a1,b1)*D(b1,a1),
         the bins are computed from the peak frequencies */
      for (b=0;b<n;b++) {
        double nu = EL[a]-ER[b];
        int bin1 = binIndex(-(nu1 + scheme1*nu)*dt1,nPoints1);
        int bin2 = bin2PerPeak ? binIndex(-(nu2 + scheme2*nu)*dt2,nPoints2) : 0;
        addPeak(&column[bin1 + (size_t)bin2*nPoints1],cmul(G[a+b*n],D[b+a*n]));
      }
    }
    else if (bins1 && !bin2PerPeak) {
      /* common case: tabulated dimension-1 bins, one buffer column */
      for (b=0;b<n;b++)
        addPeak(&column[bins1[b]],cmul(ampA,products[b]));
    }
    else {
      /* bins computed from the peak frequencies */
      for (b=0;b<n;b++) {
        double nu = EL[a]-ER[b];
        int bin1 = bins1 ? bins1[b] : binIndex(-(nu1 + scheme1*nu)*dt1,nPoints1);
        int bin2 = bin2PerPeak ? binIndex(-(nu2 + scheme2*nu)*dt2,nPoints2) : 0;
        addPeak(&column[bin1 + (size_t)bin2*nPoints1],cmul(ampA,products[b]));
      }
    }
  }
}

/* Loops over the states of the interval at loop position pos, then over
   the intervals inside it. amp is the amplitude of the chain so far, nu1
   and nu2 the frequency sums of the intervals so far. */
static void loopInterval(Pathway *P, int pos, Complex amp, double nu1, double nu2)
{
  int n = P->nStates, k, a, b;

  if (pos==P->nIntervals-1) {
    innermostInterval(P,amp,nu1,nu2);
    return;
  }

  k = P->loopOrder[pos];
  for (a=0;a<n;a++) {
    Complex ampA;
    P->state[LEFT(k)] = a;
    ampA = multiplyChain(P,LEFT(k),amp);
    if (isZero(ampA)) continue;
    for (b=0;b<n;b++) {
      Complex ampB;
      double nu = P->EL[k][a] - P->ER[k][b];
      P->state[RIGHT(k)] = b;
      ampB = multiplyChain(P,RIGHT(k),ampA);
      if (isZero(ampB)) continue;
      loopInterval(P,pos+1,ampB,nu1 + P->incScheme[k][0]*nu,nu2 + P->incScheme[k][1]*nu);
    }
  }
}

/*-------------------------------------------------------------------*/
/* Assigns chain matrix M(row,col) to the inner one of the loops over its
   two states, the first point where both states are set */
static void addChainMatrix(Pathway *P, const int *loopDepth, const Complex *M, int row, int col)
{
  int s = (loopDepth[row]>loopDepth[col]) ? row : col;
  if (P->nChain[s]>=2) mexErrMsgTxt("Internal error: more than two chain matrices per state.");
  P->chain[s][P->nChain[s]].M = M;
  P->chain[s][P->nChain[s]].row = row;
  P->chain[s][P->nChain[s]].col = col;
  P->nChain[s]++;
}

/* Element of chain matrix m for state b at position s in state[] and
   state u at its other position */
static Complex chainElement(const ChainMatrix *m, int s, int b, int u, int n)
{
  return (m->row==s) ? m->M[b + u*n] : m->M[u + b*n];
}

/* Complex n x n matrix from an input argument; real input is copied with
   zero imaginary parts (mxCalloc, freed when the MEX function returns) */
static const Complex *getMatrix(const mxArray *A, int n, const char *name)
{
  size_t q, N = mxGetNumberOfElements(A);
  const double *re;
  Complex *X;
  if (!mxIsDouble(A) || N!=(size_t)n*n)
    mexErrMsgIdAndTxt("EasySpin:sf_peaks","%s must be a double array with nStates x nStates elements.",name);
  if (mxIsComplex(A)) return mxGetComplexDoubles(A);
  re = mxGetDoubles(A);
  X = mxCalloc(N,sizeof(Complex));
  for (q=0;q<N;q++) X[q].real = re[q];
  return X;
}

/*===================================================================*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  Pathway P;
  int i, k, pos, d, K, n, loopDepth[2*MAX_INTERVALS];
  const double *dt, *incScheme, *freeL, *freeR, *Ea, *Eb;
  const Complex *G, *D, *ML[MAX_INTERVALS], *MR[MAX_INTERVALS];
  Complex unitAmplitude = {1.0, 0.0};

  /* Inputs
     --------------------------------------------------------------- */
  if (nrhs<9)
    mexErrMsgTxt("Insufficient number of input arguments.");
  if (!mxIsDouble(prhs[0]) || !mxIsComplex(prhs[0]))
    mexErrMsgTxt("buffer must be a complex double array.");
  for (i=1;i<=6;i++)
    if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]))
      mexErrMsgTxt("dt, IncScheme, idxFreeL, idxFreeR, Ea and Eb must be real double arrays.");

  P.buffer = mxGetComplexDoubles(prhs[0]);
  dt = mxGetDoubles(prhs[1]);
  P.nDims = (int)mxGetNumberOfElements(prhs[1]);
  if (P.nDims<1 || P.nDims>2)
    mexErrMsgTxt("dt must have 1 or 2 elements.");
  incScheme = mxGetDoubles(prhs[2]);
  K = (int)mxGetM(prhs[2]);
  if ((int)mxGetN(prhs[2])!=P.nDims)
    mexErrMsgTxt("IncScheme must have one column per element of dt.");
  if (K<1 || K>MAX_INTERVALS)
    mexErrMsgTxt("IncScheme must have 1 to 4 rows (incremented intervals).");
  if (nrhs!=9+2*(K-1))
    mexErrMsgTxt("Wrong number of transfer matrices.");
  P.nIntervals = K;

  if (P.nDims==1) {
    P.nPoints[0] = (int)mxGetNumberOfElements(prhs[0]);
    P.nPoints[1] = 1;
  }
  else {
    P.nPoints[0] = (int)mxGetM(prhs[0]);
    P.nPoints[1] = (int)mxGetN(prhs[0]);
  }
  if (P.nPoints[0]<1 || P.nPoints[1]<1)
    mexErrMsgTxt("buffer must not be empty.");
  for (d=0;d<2;d++)
    P.dt[d] = (d<P.nDims) ? dt[d] : 0;

  freeL = mxGetDoubles(prhs[3]);
  freeR = mxGetDoubles(prhs[4]);
  if ((int)mxGetNumberOfElements(prhs[3])!=K || (int)mxGetNumberOfElements(prhs[4])!=K)
    mexErrMsgTxt("idxFreeL and idxFreeR must have one element per row of IncScheme.");
  for (k=0;k<K;k++)
    if ((freeL[k]!=1 && freeL[k]!=2) || (freeR[k]!=1 && freeR[k]!=2))
      mexErrMsgTxt("idxFreeL and idxFreeR must be 1 (alpha) or 2 (beta).");
  Ea = mxGetDoubles(prhs[5]);
  Eb = mxGetDoubles(prhs[6]);
  n = (int)mxGetNumberOfElements(prhs[5]);
  if ((int)mxGetNumberOfElements(prhs[6])!=n)
    mexErrMsgTxt("Ea and Eb must have the same number of elements.");
  P.nStates = n;

  /* with finite energies and time increments, binIndex keeps all bins
     within the buffer, so the binning loops need no range checks */
  for (i=0;i<n;i++)
    if (!mxIsFinite(Ea[i]) || !mxIsFinite(Eb[i]))
      mexErrMsgTxt("Ea and Eb must be finite.");
  for (d=0;d<P.nDims;d++)
    if (!mxIsFinite(dt[d]))
      mexErrMsgTxt("dt must be finite.");
  for (i=0;i<K*P.nDims;i++)
    if (!mxIsFinite(incScheme[i]))
      mexErrMsgTxt("IncScheme must be finite.");

  G = getMatrix(prhs[7],n,"G");
  D = getMatrix(prhs[8],n,"D");
  P.G = G;
  P.D = D;
  for (k=0;k<K-1;k++) {
    ML[k] = getMatrix(prhs[9+k],n,"MkL");
    MR[k] = getMatrix(prhs[9+(K-1)+k],n,"MkR");
  }

  for (k=0;k<K;k++) {
    P.EL[k] = (freeL[k]==1) ? Ea : Eb;
    P.ER[k] = (freeR[k]==1) ? Ea : Eb;
    for (d=0;d<2;d++)
      P.incScheme[k][d] = (d<P.nDims) ? incScheme[k+d*K] : 0;
  }

  /* Loop order: intervals that change dimension 2 outside, the others
     inside, each in sequence order
     --------------------------------------------------------------- */
  pos = 0;
  for (k=0;k<K;k++) if (P.incScheme[k][1]!=0) P.loopOrder[pos++] = k;
  for (k=0;k<K;k++) if (P.incScheme[k][1]==0) P.loopOrder[pos++] = k;
  for (pos=0;pos<K;pos++) {
    loopDepth[LEFT(P.loopOrder[pos])] = 2*pos;      /* loop over a */
    loopDepth[RIGHT(P.loopOrder[pos])] = 2*pos+1;   /* loop over b, inside */
  }
  P.bin2PerPeak = P.nDims==2 && P.incScheme[P.loopOrder[K-1]][1]!=0;

  /* Chain matrices, each assigned to the inner one of its two loops
     --------------------------------------------------------------- */
  for (i=0;i<2*K;i++) P.nChain[i] = 0;
  addChainMatrix(&P,loopDepth,G,LEFT(0),RIGHT(0));               /* G(a1,b1) */
  for (k=0;k<K-1;k++) {
    addChainMatrix(&P,loopDepth,ML[k],LEFT(k+1),LEFT(k));        /* MkL(a(k+1),ak) */
    addChainMatrix(&P,loopDepth,MR[k],RIGHT(k),RIGHT(k+1));      /* MkR(bk,b(k+1)) */
  }
  addChainMatrix(&P,loopDepth,D,RIGHT(K-1),LEFT(K-1));           /* D(bK,aK) */

  /* Products of the two chain matrices that contain the right state b of
     the innermost interval (K >= 2)
     --------------------------------------------------------------- */
  P.rightProduct = NULL;
  if (K>=2) {
    int s = RIGHT(P.loopOrder[K-1]), u, v, b;
    const ChainMatrix *m1 = &P.chain[s][0], *m2 = &P.chain[s][1];
    if (P.nChain[s]!=2)
      mexErrMsgTxt("Internal error: the innermost right state must be in two chain matrices.");
    P.rightOther[0] = (m1->row==s) ? m1->col : m1->row;
    P.rightOther[1] = (m2->row==s) ? m2->col : m2->row;
    P.rightProduct = mxMalloc((size_t)n*n*n*sizeof(Complex));
    for (u=0;u<n;u++)
      for (v=0;v<n;v++)
        for (b=0;b<n;b++)
          P.rightProduct[((size_t)u*n+v)*n+b] = cmul(chainElement(m1,s,b,u,n),chainElement(m2,s,b,v,n));
    P.nChain[s] = 0;  /* multiplied in via rightProduct */
  }

  /* Dimension-1 bins of the innermost one or two intervals (K >= 2),
     tabulated with the same expression as in innermostInterval, so the
     bins are identical
     --------------------------------------------------------------- */
  P.bins1 = NULL;
  P.bins1Intervals = 0;
  if (K>=2) {
    int first1 = K-1;  /* loop position of the first interval that changes dimension 1 */
    int inner = P.loopOrder[K-1], a, b, a2, b2;
    for (pos=0;pos<K;pos++)
      if (P.incScheme[P.loopOrder[pos]][0]!=0) { first1 = pos; break; }
    if (first1==K-1) {
      /* only the innermost interval: n x n table */
      P.bins1Intervals = 1;
      P.bins1 = mxMalloc((size_t)n*n*sizeof(int));
      for (a=0;a<n;a++)
        for (b=0;b<n;b++)
          P.bins1[a*n+b] = binIndex(-P.incScheme[inner][0]*(P.EL[inner][a]-P.ER[inner][b])*P.dt[0],P.nPoints[0]);
    }
    else if (first1==K-2 && K>=3) {
      /* the two innermost intervals, reused for the states of the outer
         ones: n^4 table, one n x n block per state pair (a2,b2) of the
         second-innermost interval */
      int outer = P.loopOrder[K-2];
      P.bins1Intervals = 2;
      P.bins1 = mxMalloc((size_t)n*n*n*n*sizeof(int));
      for (a2=0;a2<n;a2++)
        for (b2=0;b2<n;b2++) {
          double nu1 = P.incScheme[outer][0]*(P.EL[outer][a2]-P.ER[outer][b2]);
          int *block = P.bins1 + (size_t)(a2*n+b2)*n*n;
          for (a=0;a<n;a++)
            for (b=0;b<n;b++)
              block[a*n+b] = binIndex(-(nu1 + P.incScheme[inner][0]*(P.EL[inner][a]-P.ER[inner][b]))*P.dt[0],P.nPoints[0]);
        }
    }
  }

  /* Peaks
     --------------------------------------------------------------- */
  loopInterval(&P,0,unitAmplitude,0.0,0.0);

  if (P.bins1) mxFree(P.bins1);
  if (P.rightProduct) mxFree(P.rightProduct);
}
