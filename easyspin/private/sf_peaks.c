/*
sf_peaks  Bin all peaks of a coherence transfer pathway into a spectral buffer

  sf_peaks(buffer,dt,IncScheme,idxL,idxR,Ea,Eb,G,D)
  sf_peaks(...,M1L,...,M(K-1)L,M1R,...,M(K-1)R)

  Input:
    buffer      complex frequency-domain buffer, modified in place: vector
                (1 dimension) or nPoints1 x nPoints2 matrix (2 dimensions)
    dt          time increment of each dimension (microseconds)
    IncScheme   incrementation scheme, K x nDimensions: interval k gets
                longer by IncScheme(k,d)*dt(d) per step along dimension d;
                one row per incremented interval (K <= 4), in sequence
                order. Examples:
                  3p ESEEM     [1]
                  2p ESEEM     [1; 1]
                  HYSCORE      [1 0; 0 1]
                  2D-CP        [1 0; 0 1; 0 1; 1 0]
                  Carr-Purcell [1; 2; 1]
                An interval can change along more than one dimension.
    idxL, idxR  left and right manifold index (1 = alpha, 2 = beta) of
                each incremented interval
    Ea, Eb      energies of the alpha and beta manifolds (MHz)
    G           density matrix before the first incremented interval
    D           detection operator matrix
    MkL, MkR    left and right transfer matrices between intervals k and
                k+1 (G, D, MkL, MkR: nStates x nStates)

  During interval k, the density matrix element rho(l,r), with left
  state l in manifold idxL(k) and right state r in manifold idxR(k),
  oscillates with the frequency EL_k(l)-ER_k(r), where EL_k and ER_k are
  the energies (Ea or Eb) of these two manifolds. Transfers between
  intervals are rho -> MkL*rho*MkR, and the signal is trace(D*rho).
  Expanded, the signal is a sum of peaks, one for each choice of states
  l1,r1,...,lK,rK, with the amplitude

    D(rK,lK) * M(K-1)L(lK,l(K-1)) * ... * M1L(l2,l1) * G(l1,r1)
             * M1R(r1,r2) * ... * M(K-1)R(r(K-1),rK)

  and the frequency nu_d = sum_k IncScheme(k,d)*(EL_k(lk)-ER_k(rk)) along
  dimension d. There are nStates^(2K) peaks. For more than 4 incremented
  intervals, direct point-by-point evaluation in the time domain is more
  efficient (see Opt.Binning in saffron).

  A peak goes into bin fmod(-nu*dt,1)*nPoints, rounded to the nearest integer
  with ties away from zero, negative bins wrapped around.
  The - sign accounts for quantum evolution going with exp(-i*omega*t),
  whereas Fourier theory uses exp(+i*omega*t). The rounding is odd in nu,
  so conjugate peak pairs (nu,A), (-nu,conj(A)) give exactly a
  conjugate-symmetric spectrum.

  Algorithm: Nested loops over the left and right state of each interval.
  Each matrix element of the amplitude is multiplied in within the inner
  one of the two loops over its states. Intervals that change dimension 2
  are looped outside, so that the innermost loops add to one buffer column.

  Shortcuts: peaks with zero amplitude are skipped early (pulses that keep
  the electron in its manifold give diagonal transfer matrices). For
  K >= 2, where the outer loops reuse them, the two matrix elements with
  the innermost right state are tabulated as one product, and the
  dimension-1 bins are tabulated if only the innermost interval, or (for
  K >= 3) the two innermost intervals, change dimension 1. For K = 1, each
  table entry would be used once, so nothing is tabulated.

  Tabulated bins vs bins computed per peak: the bin along a dimension
  depends on the states of all intervals that change that dimension. If
  these are only the innermost one or two intervals, the bins form a table
  over their states (nStates^2 or nStates^4 entries). The table is the same
  for all states of the outer intervals, so it is computed once per call,
  and each peak looks up its bin. Otherwise, the bin also depends on the
  outer states, a table would have as many entries as there are peaks, and
  binIndex is evaluated for each peak. Examples are 2p ESEEM [1; 1],
  Carr-Purcell [1; 2; 1], and an innermost interval that changes dimension
  2. Along dimension 2, the bin is otherwise computed once per buffer
  column. Tables and per-peak evaluation use the same expression, so the
  bins are identical.
 */

#include "mex.h"
#include <math.h>

/* For K incremented intervals, there are nStates^(2K) peaks. Binning is
   limited to K <= 4; for larger K, time-domain evaluation is more efficient */
#define MAX_INTERVALS 4

/* Each interval has a left state l and a right state r. The current
   states are kept in one array, state[]: state[LEFT(k)] is the left state
   of interval k, state[RIGHT(k)] its right state. */
#define LEFT(k)  (2*(k))
#define RIGHT(k) (2*(k)+1)

typedef mxComplexDouble Complex;

static Complex cmul(Complex x, Complex y) {
  Complex z;
  z.real = x.real*y.real - x.imag*y.imag;
  z.imag = x.real*y.imag + x.imag*y.real;
  return z;
}

static int isZero(Complex x) {
  return x.real==0 && x.imag==0;
}

/* adds the amplitude of a peak to its bin */
static void addPeak(Complex *bin, Complex peak) {
  bin->real += peak.real;
  bin->imag += peak.imag;
}

/* One matrix of the pathway (G, D, or a transfer matrix MkL, MkR), with
   the positions in state[] of the states that select its row and its
   column, e.g. G(l1,r1): row LEFT(0), col RIGHT(0) */
typedef struct {
  const Complex *M;
  int row, col;
} PathwayMatrix;

typedef struct {
  int nStates, nIntervals;
  /* In 1D, dimension 2 has nPoints 1, dt 0 and IncScheme 0, so its bin is 0 */
  int nPoints[2];                /* buffer size along each dimension */
  double dt[2];                  /* time increment of each dimension (µs) */
  Complex *buffer;

  /* Incremented intervals, numbered in sequence order */
  const double *EL[MAX_INTERVALS];   /* energies of the left manifold (MHz) */
  const double *ER[MAX_INTERVALS];   /* energies of the right manifold (MHz) */
  double incScheme[MAX_INTERVALS][2];  /* incrementation scheme, IncScheme(k,d) */
  int loopOrder[MAX_INTERVALS];      /* intervals from outer to inner loop */

  /* Current states, and the matrices multiplied in within the loop over
     each state (at most two per state) */
  PathwayMatrix matrices[2*MAX_INTERVALS][2];
  int nMatrices[2*MAX_INTERVALS];
  int state[2*MAX_INTERVALS];

  /* Innermost interval: its right state r appears in two matrix elements
     of the amplitude. For K >= 2, their product is tabulated over r and
     the other states u, v of the two elements (u = state[rightOther[0]],
     v = state[rightOther[1]]) at rightProduct[(u*nStates+v)*nStates+r].
     Example HYSCORE: the innermost interval is t1, r1 appears in G(l1,r1)
     and M1R(r1,r2), so u = l1 and v = r2. NULL for K = 1: the two
     elements G(l1,r1), D(r1,l1) are multiplied in per peak. */
  Complex *rightProduct;
  int rightOther[2];
  const Complex *G, *D;          /* density matrix before the first interval, detection operator */

  /* Tabulated dimension-1 bins of the innermost one or two intervals, or
     NULL (always for K = 1); with two intervals, one n x n block per state
     pair of the second-innermost interval */
  int *bins1;
  int bins1Intervals;

  int bin2PerPeak;               /* 1 if the innermost interval changes dimension 2 */
} Pathway;

/* Bin of a peak with phase x = -nu*dt (in cycles): fmod(x,1)*nPoints,
   rounded to the nearest integer with ties away from zero, negative bins
   wrapped around; for finite x and nPoints>0, 0 <= bin < nPoints. The
   fractional part x-trunc(x) is exact and equal to fmod(x,1), but faster;
   doubles with |x| >= 2^52 are integers (fractional part 0), and the cast
   to long long needs |x| < 2^63. y-b is exact, so the rounding is exact
   and odd in x. */
static int binIndex(double x, int nPoints) {
  double frac = (fabs(x)<4503599627370496.0) ? x - (double)(long long)x : 0.0;  /* 2^52 */
  double y = frac*nPoints;
  int b = (int)y;
  if (y-b>=0.5) b++;
  else if (y-b<=-0.5) b--;
  if (b<0) b += nPoints;
  if (b>=nPoints) b -= nPoints;
  return b;
}

/* Multiplies amp by the matrix elements assigned to the loop over state[s]
   (all their states are set) */
static Complex multiplyMatrices(const Pathway *P, int s, Complex amp) {
  int i, n = P->nStates;
  for (i=0;i<P->nMatrices[s];i++) {
    const PathwayMatrix *m = &P->matrices[s][i];
    amp = cmul(amp,m->M[P->state[m->row] + P->state[m->col]*n]);
  }
  return amp;
}

/* Loops over the states of the innermost interval and adds the peaks to
   the buffer. amp is the amplitude of the outer intervals, nu1 and nu2
   the frequency sums of the outer intervals. */
static void innermostInterval(Pathway *P, Complex amp, double nu1, double nu2) {
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
  int l, r;

  /* The dimension-2 frequency is complete unless this interval changes it */
  if (!bin2PerPeak)
    column += (size_t)binIndex(-nu2*dt2,nPoints2)*nPoints1;

  /* bins tabulated over two intervals: the block for the current states
     of the second-innermost interval */
  if (P->bins1Intervals==2) {
    int outer = P->loopOrder[P->nIntervals-2];
    bins1Block += ((size_t)P->state[LEFT(outer)]*n + P->state[RIGHT(outer)])*n*n;
  }

  for (l=0;l<n;l++) {
    const Complex *products = NULL;
    const int *bins1 = NULL;
    Complex ampL;
    P->state[LEFT(k)] = l;
    ampL = multiplyMatrices(P,LEFT(k),amp);
    if (isZero(ampL)) continue;
    /* rows of the tables for the current states (NULL if not tabulated) */
    if (P->rightProduct) {
      int u = P->state[P->rightOther[0]], v = P->state[P->rightOther[1]];
      products = P->rightProduct + ((size_t)u*n + v)*n;
    }
    if (bins1Block) bins1 = bins1Block + (size_t)l*n;

    if (P->nIntervals==1) {
      /* K = 1, nothing tabulated: the peak amplitude is G(l1,r1)*D(r1,l1),
         the bins are computed from the peak frequencies */
      for (r=0;r<n;r++) {
        double nu = EL[l]-ER[r];
        int bin1 = binIndex(-(nu1 + scheme1*nu)*dt1,nPoints1);
        int bin2 = bin2PerPeak ? binIndex(-(nu2 + scheme2*nu)*dt2,nPoints2) : 0;
        addPeak(&column[bin1 + (size_t)bin2*nPoints1],cmul(G[l+r*n],D[r+l*n]));
      }
    }
    else if (bins1 && !bin2PerPeak) {
      /* tabulated dimension-1 bins, one buffer column (e.g. HYSCORE, 2D-CP) */
      for (r=0;r<n;r++)
        addPeak(&column[bins1[r]],cmul(ampL,products[r]));
    }
    else {
      /* bins computed from the peak frequencies (e.g. 2p ESEEM,
         Carr-Purcell), except tabulated dimension-1 bins */
      for (r=0;r<n;r++) {
        double nu = EL[l]-ER[r];
        int bin1 = bins1 ? bins1[r] : binIndex(-(nu1 + scheme1*nu)*dt1,nPoints1);
        int bin2 = bin2PerPeak ? binIndex(-(nu2 + scheme2*nu)*dt2,nPoints2) : 0;
        addPeak(&column[bin1 + (size_t)bin2*nPoints1],cmul(ampL,products[r]));
      }
    }
  }
}

/* Loops (recursively) over the states of the interval at loop position pos,
   then over the intervals inside it. amp is the amplitude so far, nu1 and nu2
   are the frequency sums of the intervals so far. */
static void loopInterval(Pathway *P, int pos, Complex amp, double nu1, double nu2) {
  int n = P->nStates, k, l, r;

  if (pos==P->nIntervals-1) {
    innermostInterval(P,amp,nu1,nu2);
    return;
  }

  k = P->loopOrder[pos];
  for (l=0;l<n;l++) {
    Complex ampL;
    P->state[LEFT(k)] = l;
    ampL = multiplyMatrices(P,LEFT(k),amp);
    if (isZero(ampL)) continue;
    for (r=0;r<n;r++) {
      Complex ampR;
      double nu = P->EL[k][l] - P->ER[k][r];
      P->state[RIGHT(k)] = r;
      ampR = multiplyMatrices(P,RIGHT(k),ampL);
      if (isZero(ampR)) continue;
      loopInterval(P,pos+1,ampR,nu1 + P->incScheme[k][0]*nu,nu2 + P->incScheme[k][1]*nu);
    }
  }
}

/* Assigns matrix M(row,col) to the inner one of the loops over its two
   states, the first point where both states are set */
static void addPathwayMatrix(Pathway *P, const int *loopDepth, const Complex *M, int row, int col) {
  int s = (loopDepth[row]>loopDepth[col]) ? row : col;
  PathwayMatrix *m;
  if (P->nMatrices[s]>=2) mexErrMsgTxt("Internal error: more than two pathway matrices per state.");
  m = &P->matrices[s][P->nMatrices[s]];
  m->M = M;
  m->row = row;
  m->col = col;
  P->nMatrices[s]++;
}

/* Element of matrix m for state r at position s in state[] and state u
   at its other position */
static Complex matrixElement(const PathwayMatrix *m, int s, int r, int u, int n) {
  return (m->row==s) ? m->M[r + u*n] : m->M[u + r*n];
}

/* Complex n x n matrix from an input argument; real input is copied with
   zero imaginary parts (mxCalloc, freed when the MEX function returns) */
static const Complex *getMatrix(const mxArray *A, int n, const char *name) {
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
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  Pathway P;
  int i, k, pos, d, K, n, nDims, loopDepth[2*MAX_INTERVALS];
  const double *dt, *incScheme, *idxL, *idxR, *Ea, *Eb;
  const Complex *ML[MAX_INTERVALS], *MR[MAX_INTERVALS];
  Complex unitAmplitude = {1.0, 0.0};

  /* Check input arguments */
  if (nrhs<9)
    mexErrMsgTxt("At least 9 input arguments are required.");
  if (!mxIsDouble(prhs[0]) || !mxIsComplex(prhs[0]))
    mexErrMsgTxt("buffer (1st input) must be a complex double array.");
  for (i=1;i<=6;i++)
    if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]))
      mexErrMsgTxt("dt, IncScheme, idxL, idxR, Ea and Eb (inputs 2-7) must be real double arrays.");

  /* Read and check buffer, dt, incrementation scheme, number of points */
  P.buffer = mxGetComplexDoubles(prhs[0]);
  dt = mxGetDoubles(prhs[1]);
  nDims = (int)mxGetNumberOfElements(prhs[1]);
  if (nDims<1 || nDims>2)
    mexErrMsgTxt("dt (2nd input) must have 1 or 2 elements.");
  incScheme = mxGetDoubles(prhs[2]);
  K = (int)mxGetM(prhs[2]);
  if ((int)mxGetN(prhs[2])!=nDims)
    mexErrMsgTxt("IncScheme must have one column per element of dt.");
  if (K<1 || K>MAX_INTERVALS)
    mexErrMsgTxt("IncScheme must have 1 to 4 rows (incremented intervals).");
  if (nrhs!=9+2*(K-1))
    mexErrMsgTxt("Wrong number of transfer matrices for the given incrementation scheme.");
  P.nIntervals = K;
  if (nDims==1) {
    P.nPoints[0] = (int)mxGetNumberOfElements(prhs[0]);
    P.nPoints[1] = 1;
  }
  else {
    P.nPoints[0] = (int)mxGetM(prhs[0]);
    P.nPoints[1] = (int)mxGetN(prhs[0]);
  }
  if (P.nPoints[0]<1 || P.nPoints[1]<1)
    mexErrMsgTxt("buffer must not be empty.");
  /* In 1D, dimension 2 has one point and time increment zero */
  P.dt[0] = dt[0];
  P.dt[1] = (nDims==2) ? dt[1] : 0;

  /* Read and check left and right manifold of each incremented interval */
  idxL = mxGetDoubles(prhs[3]);
  idxR = mxGetDoubles(prhs[4]);
  if ((int)mxGetNumberOfElements(prhs[3])!=K || (int)mxGetNumberOfElements(prhs[4])!=K)
    mexErrMsgTxt("idxL and idxR must have one element per incremented interval.");
  for (k=0;k<K;k++)
    if ((idxL[k]!=1 && idxL[k]!=2) || (idxR[k]!=1 && idxR[k]!=2))
      mexErrMsgTxt("idxL and idxR must be 1 (alpha) or 2 (beta).");

  /* Get energies of the alpha and beta manifolds */
  Ea = mxGetDoubles(prhs[5]);
  Eb = mxGetDoubles(prhs[6]);
  n = (int)mxGetNumberOfElements(prhs[5]);
  if ((int)mxGetNumberOfElements(prhs[6])!=n)
    mexErrMsgTxt("Ea and Eb must have the same number of elements.");
  P.nStates = n;

  /* Assert energies and time increments are finite */
  for (i=0;i<n;i++)
    if (!mxIsFinite(Ea[i]) || !mxIsFinite(Eb[i]))
      mexErrMsgTxt("Ea and Eb must be finite.");
  for (d=0;d<nDims;d++)
    if (!mxIsFinite(dt[d]))
      mexErrMsgTxt("dt must be finite.");
  for (i=0;i<K*nDims;i++)
    if (!mxIsFinite(incScheme[i]))
      mexErrMsgTxt("IncScheme must be finite.");

  /* Get initial density matrix, detection operator, transfer matrices */
  P.G = getMatrix(prhs[7],n,"G");
  P.D = getMatrix(prhs[8],n,"D");
  for (k=0;k<K-1;k++) {
    ML[k] = getMatrix(prhs[9+k],n,"MkL");
    MR[k] = getMatrix(prhs[9+(K-1)+k],n,"MkR");
  }

  /* Get energies of the left and right manifold for each interval */
  for (k=0;k<K;k++) {
    P.EL[k] = (idxL[k]==1) ? Ea : Eb;
    P.ER[k] = (idxR[k]==1) ? Ea : Eb;
  }

  /* Copy the incrementation scheme (K x nDims, column-major) to the
     K x 2 array P.incScheme, with IncScheme(k,d) at P.incScheme[k][d].
     In 1D, dimension 2 is set to zero, so the code can always use both
     dimensions. */
  for (k=0;k<K;k++) {
    P.incScheme[k][0] = incScheme[k];
    P.incScheme[k][1] = (nDims==2) ? incScheme[k+K] : 0;
  }

  /* Determine canonical loop order (from outermost to innermost): first all
     intervals that change along dimension 2 (outer loops), then the others
     (inner loops), all in sequence order */
  pos = 0;
  for (k=0;k<K;k++)
    if (P.incScheme[k][1]!=0) {  /* dimension 2 changes interval k */
      P.loopOrder[pos] = k;
      pos++;
    }
  for (k=0;k<K;k++)
    if (P.incScheme[k][1]==0) {  /* dimension 2 does not change interval k */
      P.loopOrder[pos] = k;
      pos++;
    }

  /* Nesting depth of the loop over each state, 0 = outermost */
  for (pos=0;pos<K;pos++) {
    loopDepth[LEFT(P.loopOrder[pos])] = 2*pos;      /* loop over l */
    loopDepth[RIGHT(P.loopOrder[pos])] = 2*pos+1;   /* loop over r, inside */
  }

  /* Dimension-2 bin per peak or per buffer column: if the innermost
     interval does not change dimension 2, the dimension-2 bin is fixed
     within the innermost loops, and all their peaks go into one buffer
     column. If it changes dimension 2, the dimension-2 bin differs between
     these peaks and is computed for each peak. With the loop order above,
     this happens only if all intervals change dimension 2, e.g. [1 1] or
     [0 1; 1 1]. */
  P.bin2PerPeak = P.incScheme[P.loopOrder[K-1]][1]!=0;

  /* G, D and transfer matrices, each assigned to the inner one of its two
     loops */
  for (i=0;i<2*K;i++) P.nMatrices[i] = 0;
  addPathwayMatrix(&P,loopDepth,P.G,LEFT(0),RIGHT(0));               /* G(l1,r1) */
  for (k=0;k<K-1;k++) {
    addPathwayMatrix(&P,loopDepth,ML[k],LEFT(k+1),LEFT(k));        /* MkL(l(k+1),lk) */
    addPathwayMatrix(&P,loopDepth,MR[k],RIGHT(k),RIGHT(k+1));      /* MkR(rk,r(k+1)) */
  }
  addPathwayMatrix(&P,loopDepth,P.D,RIGHT(K-1),LEFT(K-1));           /* D(rK,lK) */

  /* Products of the two matrix elements that contain the right state r
     of the innermost interval (K >= 2) */
  P.rightProduct = NULL;
  if (K>=2) {
    int s = RIGHT(P.loopOrder[K-1]), u, v, r;
    const PathwayMatrix *m1 = &P.matrices[s][0], *m2 = &P.matrices[s][1];
    if (P.nMatrices[s]!=2)
      mexErrMsgTxt("Internal error: the innermost right state must be in two pathway matrices.");
    P.rightOther[0] = (m1->row==s) ? m1->col : m1->row;
    P.rightOther[1] = (m2->row==s) ? m2->col : m2->row;
    P.rightProduct = mxMalloc((size_t)n*n*n*sizeof(Complex));
    for (u=0;u<n;u++)
      for (v=0;v<n;v++)
        for (r=0;r<n;r++)
          P.rightProduct[((size_t)u*n+v)*n+r] = cmul(matrixElement(m1,s,r,u,n),matrixElement(m2,s,r,v,n));
    P.nMatrices[s] = 0;  /* multiplied in via rightProduct */
  }

  /* Dimension-1 bins (K >= 2). They depend on the states of the intervals
     from the outermost one that changes dimension 1 inwards. If these are
     the innermost one or two intervals, the bins are tabulated (see the
     header). The expressions are the same as in innermostInterval, so the
     bins are identical. */
  P.bins1 = NULL;
  P.bins1Intervals = 0;
  if (K>=2) {
    const int inner = P.loopOrder[K-1];
    const double *EL = P.EL[inner], *ER = P.ER[inner];
    const double scheme1 = P.incScheme[inner][0], dt1 = P.dt[0];
    const int nPoints1 = P.nPoints[0];
    int nIntervals1 = 1;  /* number of innermost intervals the bins depend on */
    int l, r, l2, r2;
    for (pos=0;pos<K-1;pos++)
      if (P.incScheme[P.loopOrder[pos]][0]!=0) {
        nIntervals1 = K-pos;
        break;
      }

    if (nIntervals1==1) {
      /* only the innermost interval: n x n table */
      P.bins1Intervals = 1;
      P.bins1 = mxMalloc((size_t)n*n*sizeof(int));
      for (l=0;l<n;l++)
        for (r=0;r<n;r++)
          P.bins1[l*n+r] = binIndex(-scheme1*(EL[l]-ER[r])*dt1,nPoints1);
    }
    else if (nIntervals1==2 && K>=3) {
      /* The two innermost intervals, reused for the states of the outer
         ones: n^4 table, one n x n block per state pair (l2,r2) of the
         second-innermost interval. For K = 2, each entry would be used
         once, so nothing is tabulated. */
      const int outer = P.loopOrder[K-2];
      P.bins1Intervals = 2;
      P.bins1 = mxMalloc((size_t)n*n*n*n*sizeof(int));
      for (l2=0;l2<n;l2++)
        for (r2=0;r2<n;r2++) {
          double nu1 = P.incScheme[outer][0]*(P.EL[outer][l2]-P.ER[outer][r2]);
          int *block = P.bins1 + (size_t)(l2*n+r2)*n*n;
          for (l=0;l<n;l++)
            for (r=0;r<n;r++)
              block[l*n+r] = binIndex(-(nu1 + scheme1*(EL[l]-ER[r]))*dt1,nPoints1);
        }
    }
  }

  /* Peaks */
  loopInterval(&P,0,unitAmplitude,0.0,0.0);

  if (P.bins1) mxFree(P.bins1);
  if (P.rightProduct) mxFree(P.rightProduct);
}
