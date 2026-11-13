/*
sf_peaks  Bin the peaks of a coherence transfer pathway into a buffer

  sf_peaks(IncSchemeID,buffer,dt,idxFreeL,idxFreeR,Ea,Eb,G,D)
  sf_peaks(...,M1L,M1R)
  sf_peaks(...,M1L,M2L,M1R,M2R)
  sf_peaks(...,M1L,M2L,M3L,M1R,M2R,M3R)

  Computes the amplitude and frequency of each peak of a coherence transfer
  pathway with incremented delays, and adds the amplitudes to the
  frequency-domain buffer at the bins of the peak frequencies.

  Input:
    IncSchemeID  incrementation scheme, the dimension of each incremented
                 delay:
                    1   [1]        3   [1 -1]      14   [1 1 2]
                    2   [1 1]     11   [1 2]       15   [1 2 2 1]
                   12   [1 2 1]                    17   [1 1 2 2]
    buffer       complex buffer: 1 dimension: vector; 2 dimensions:
                 nPoints1 x nPoints2 matrix; must be stored as complex,
                 e.g. complex(zeros(...))
    dt           time increment of each dimension (microseconds)
    idxFreeL     left-side manifold (1 = alpha, 2 = beta) of each incremented
                 delay
    idxFreeR     right-side manifold (1 = alpha, 2 = beta) of each
                 incremented delay
    Ea, Eb       energies of the alpha and beta manifolds (MHz), real
                 vectors with nStates elements
    G            state before the first incremented delay, nStates x nStates
    D            detection operator, nStates x nStates
    M1L, M2L,... left-side mixing matrices between consecutive incremented
                 delays, nStates x nStates; one less than incremented delays
    M1R, M2R,... right-side mixing matrices, same as for the left side

  Output:
    none; buffer is modified in place

  Binning: a peak with frequency nu goes into bin fmod(-nu*dt,1)*nPoints,
  rounded to the nearest integer with ties away from zero, with negative
  bins wrapped around. This rounding is odd in nu, so peaks at -nu and +nu
  land in mirrored bins and conjugate peak pairs (nu,a), (-nu,conj(a)) give
  exactly conjugate-symmetric buffers.
 */

#include "mex.h"
#include <math.h>

/*-------------------------------------------------------------------*/
/* bin of a peak with phase x = -nu*dt (in cycles): fmod(x,1)*nPoints,
   rounded to the nearest integer with ties away from zero, negative bins
   wrapped around; for finite x and nPoints>0, 0 <= bin < nPoints. The
   fractional part x-trunc(x) is exact and equal to fmod(x,1), but faster;
   doubles with |x| >= 2^52 are integers (fractional part 0), and the cast
   to long long needs |x| < 2^63. y-b is exact, so the rounding is exact and
   odd in x. */
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

/*-------------------------------------------------------------------*/
/* bins of nu = (El1[a]-Er1[b]) + (El2[c]-Er2[d]) for all states a,b,c,d,
   in the order bins[((a*nStates+b)*nStates+c)*nStates+d], allocated with
   mxMalloc */
static int *pairSumBins(const double *El1, const double *Er1,
  const double *El2, const double *Er2, int nStates, double dt, int nPoints)
{
  int a, b, c, d, q = 0;
  int *bins = mxMalloc((size_t)nStates*nStates*nStates*nStates*sizeof(int));
  for (a=0;a<nStates;a++)
    for (b=0;b<nStates;b++)
      for (c=0;c<nStates;c++)
        for (d=0;d<nStates;d++)
          bins[q++] = binIndex(-((El1[a]-Er1[b]) + (El2[c]-Er2[d]))*dt,nPoints);
  return bins;
}

/*-------------------------------------------------------------------*/
/* complex nStates x nStates matrix from an input argument; real input
   is copied with zero imaginary parts (mxCalloc, freed when the MEX
   function returns) */
static mxComplexDouble *getMatrix(const mxArray *A, int nStates, const char *name)
{
  size_t q, N = mxGetNumberOfElements(A);
  const double *re;
  mxComplexDouble *X;
  if (!mxIsDouble(A) || N!=(size_t)nStates*nStates)
    mexErrMsgIdAndTxt("EasySpin:sf_peaks","%s must be a double array with nStates x nStates elements.",name);
  if (mxIsComplex(A)) return mxGetComplexDoubles(A);
  re = mxGetDoubles(A);
  X = mxCalloc(N,sizeof(mxComplexDouble));
  for (q=0;q<N;q++) X[q].real = re[q];
  return X;
}

/*===================================================================*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int nStates, nPoints1, nPoints2, nFreeEvolutions, nDimensions, nMix;
  mxComplexDouble *G, *D, *T1l = NULL, *T1r = NULL, *T2l = NULL, *T2r = NULL, *T3l = NULL, *T3r = NULL;
  mxComplexDouble *spec, *col;
  double *Ea, *Eb, *E1left, *E1right, *E2left, *E2right, *E3left, *E3right, *E4left, *E4right;
  double rAmp, iAmp, rAmp0,iAmp0, rAmp1, iAmp1, rAmp2, iAmp2, rAmp3, iAmp3, rAmp4, iAmp4, rAmp5, iAmp5;
  double *freeL, *freeR;
  int id2, idx, q, *nu1, *nu2, *bins1;

  int IncSchemeID;
  int i, j, k, l, ij, jk, kl, li, kj, il;
  int m, n, o, p, lm, mn, no, op, pi, ni;
  double *dt, dt1, dt2;
  
  /*-----------------------------------------------------------------
     Check number of input arguments
    ----------------------------------------------------------------- */
  if (nrhs<9)
    mexErrMsgTxt("Insufficient number of input arguments.");
  
  /*-----------------------------------------------------------------
     Read input arguments
    ----------------------------------------------------------------- */

  /* IncSchemeID ... identifies the incrementation scheme */
  IncSchemeID = mxGetScalar(prhs[0]);
  if (mxGetNumberOfElements(prhs[0])!=1)
    mexErrMsgTxt("IncSchemeID must be a scalar!");
  switch (IncSchemeID) {
    case  1: nFreeEvolutions = 1; nDimensions = 1; break;
    case  2: nFreeEvolutions = 2; nDimensions = 1; break;
    case  3: nFreeEvolutions = 2; nDimensions = 1; break;
    case 11: nFreeEvolutions = 2; nDimensions = 2; break;
    case 12: nFreeEvolutions = 3; nDimensions = 2; break;
    case 14: nFreeEvolutions = 3; nDimensions = 2; break;
    case 15: nFreeEvolutions = 4; nDimensions = 2; break;
    case 17: nFreeEvolutions = 4; nDimensions = 2; break;
    default: mexErrMsgTxt("Unrecognized incrementation scheme.");
  }
  nMix = nFreeEvolutions-1;
  if (nrhs!=9+2*nMix)
    mexErrMsgTxt("Wrong number of input arguments.");

  /* buffer ... complex frequency-domain buffer, modified in place; with
     the interleaved complex API, MATLAB passes complex arrays without
     copying them */
  if (!mxIsDouble(prhs[1]) || !mxIsComplex(prhs[1]))
    mexErrMsgTxt("buffer must be a complex double array.");
  spec = mxGetComplexDoubles(prhs[1]);
  if (nDimensions==1) {
    nPoints1 = mxGetM(prhs[1]);
    if (nPoints1==1) nPoints1 = mxGetN(prhs[1]);
    nPoints2 = 1;
  }
  else {
    nPoints1 = mxGetM(prhs[1]);
    nPoints2 = mxGetN(prhs[1]);
  }
  if ((nPoints1<1)||(nPoints2<1))
    mexErrMsgTxt("buffer must not be empty.");

  /* dt, freeL, freeR, Ea, Eb are real double arrays; mxGetDoubles returns
     NULL for anything else */
  for (i=2;i<=6;i++)
    if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]))
      mexErrMsgTxt("dt, idxFreeL, idxFreeR, Ea and Eb must be real double arrays.");

  /* dt ... time increment */
  dt = mxGetDoubles(prhs[2]);
  if (mxGetNumberOfElements(prhs[2])!=nDimensions)
    mexErrMsgTxt("dt has wrong number of elements!");
  dt1 = dt[0];
  if (nDimensions==2) dt2 = dt[1];

  /* freeL ... index for left-side propagators (1=alpha,2=beta) */
  freeL = mxGetDoubles(prhs[3]);
  if (mxGetNumberOfElements(prhs[3])!=nFreeEvolutions)
    mexErrMsgTxt("Wrong number of left evolution intervals.");

  /* freeR ... index for right-side propagators (1=alpha,2=beta) */
  freeR = mxGetDoubles(prhs[4]);
  if (mxGetNumberOfElements(prhs[4])!=nFreeEvolutions)
    mexErrMsgTxt("Wrong number of right evolution intervals.");

  /* Ea, Eb ... energies of the alpha and beta manifolds */
  Ea = mxGetDoubles(prhs[5]);
  Eb = mxGetDoubles(prhs[6]);
  nStates = mxGetNumberOfElements(prhs[5]);
  if (nStates!=mxGetNumberOfElements(prhs[6]))
    mexErrMsgTxt("Ea and Eb must have the same number of elements.");

  /* with finite energies and time increments, binIndex keeps all bins
     within the buffer, so the binning loops need no range checks */
  for (i=0;i<nStates;i++)
    if (!mxIsFinite(Ea[i]) || !mxIsFinite(Eb[i]))
      mexErrMsgTxt("Ea and Eb must be finite.");
  for (i=0;i<nDimensions;i++)
    if (!mxIsFinite(dt[i]))
      mexErrMsgTxt("dt must be finite.");

  /* G ... state before the first incremented delay, D ... detection
     operator, then the mixing matrices M1L, M2L, ..., M1R, M2R, ... */
  G = getMatrix(prhs[7],nStates,"G");
  D = getMatrix(prhs[8],nStates,"D");
  if (nMix>=1) {
    T1l = getMatrix(prhs[9],nStates,"M1L");
    T1r = getMatrix(prhs[9+nMix],nStates,"M1R");
  }
  if (nMix>=2) {
    T2l = getMatrix(prhs[10],nStates,"M2L");
    T2r = getMatrix(prhs[10+nMix],nStates,"M2R");
  }
  if (nMix>=3) {
    T3l = getMatrix(prhs[11],nStates,"M3L");
    T3r = getMatrix(prhs[11+nMix],nStates,"M3R");
  }

  switch (IncSchemeID) {

    case 1: /* [1], e.g. three-pulse ESEEM */
      /*-------------------------------------------------------------*/
      /* Compute amplitudes and bin peaks */
      /* Amp = T[kj]*S[jk]; */
      E1left  = (freeL[0]==1) ? Ea : Eb;
      E1right = (freeR[0]==1) ? Ea : Eb;
      for (j=0;j<nStates;j++) {
        for (k=0;k<nStates;k++) {
          jk = j+k*nStates;
          kj = k+j*nStates;
          
          /* the - signs here account for the fact that the quantum
           * evolution goes with (-omega_ij)t, but Fourier theory
           * uses (+omega_ij)t */
          idx = binIndex(-(E1left[j]-E1right[k])*dt1,nPoints1);
      
          /* compute amplitude T[kj]*S[jk] */
          rAmp = D[kj].real*G[jk].real - D[kj].imag*G[jk].imag;
          iAmp = D[kj].real*G[jk].imag + D[kj].imag*G[jk].real;
          
          spec[idx].real += rAmp;
          spec[idx].imag += iAmp;
        }
      }
      break;
    
    case 2: /* [1 1], e.g. two-pulse ESEEM, 1D-CP */
      /*-------------------------------------------------------------*/
      /* compute amplitudes and frequencies and bin peaks */
      /* Amp = M1L[ij]*S[jk]*M1R[kl]*T[li]; */
      E1left  = (freeL[0]==1) ? Ea : Eb;
      E1right = (freeR[0]==1) ? Ea : Eb;
      E2left  = (freeL[1]==1) ? Ea : Eb;
      E2right = (freeR[1]==1) ? Ea : Eb;
      for (j=0;j<nStates;j++) {
        for (k=0;k<nStates;k++) {
          jk = j+k*nStates;
          for (i=0;i<nStates;i++) {
            ij = i+j*nStates;
            rAmp = T1l[ij].real*G[jk].real - T1l[ij].imag*G[jk].imag;
            iAmp = T1l[ij].real*G[jk].imag + T1l[ij].imag*G[jk].real;
            for (l=0;l<nStates;l++) {
              kl = k+l*nStates;
              li = l+i*nStates;
              il = i+l*nStates;

              idx = binIndex(-((E2left[i]-E2right[l])+(E1left[j]-E1right[k]))*dt1,nPoints1);

              rAmp1 = rAmp*T1r[kl].real - iAmp*T1r[kl].imag;
              iAmp1 = rAmp*T1r[kl].imag + iAmp*T1r[kl].real;
              rAmp2 = rAmp1*D[li].real - iAmp1*D[li].imag;
              iAmp2 = rAmp1*D[li].imag + iAmp1*D[li].real;
              
              spec[idx].real += rAmp2;
              spec[idx].imag += iAmp2;
            }
          }
        }
      }
      break;
      
    case 11: /* [1 2], e.g. HYSCORE*/
      /*-------------------------------------------------------------*/
      E1left  = (freeL[0]==1) ? Ea : Eb;
      E1right = (freeR[0]==1) ? Ea : Eb;
      E2left  = (freeL[1]==1) ? Ea : Eb;
      E2right = (freeR[1]==1) ? Ea : Eb;
      
      nu1 = mxCalloc(nStates*nStates,sizeof(int));
      nu2 = mxCalloc(nStates*nStates,sizeof(int));
      if ((!nu1) || (!nu2))
        mexErrMsgTxt("Could not allocate memory for binning.");

      /* pre-compute nuclear frequencies */
      idx = 0;
      for (j=0;j<nStates;j++) {
        for (i=0;i<nStates;i++) {

          /* the - signs here account for the fact that the quantum
           * evolution goes with (-omega_ij)t, but Fourier theory
           * uses (+omega_ij)t */
          nu1[idx] = binIndex(-(E1left[i]-E1right[j])*dt1,nPoints1);
          nu2[idx] = binIndex(-(E2left[i]-E2right[j])*dt2,nPoints2);

          idx++;
        }
      }

      /* compute amplitudes and bin peaks */
      /* Amp = Mp[ij]*S[jk]*Mm[kl]*T[li]; */
      /* the outer loops (i,l) fix the buffer column (dimension 2), the
         inner loops (j,k) stay within it */
      for (i=0;i<nStates;i++) {
        for (l=0;l<nStates;l++) {
          li = l+i*nStates;
          il = i+l*nStates;
          col = spec + nu2[il]*nPoints1;
          for (j=0;j<nStates;j++) {
            ij = i+j*nStates;
            rAmp = D[li].real*T1l[ij].real - D[li].imag*T1l[ij].imag;
            iAmp = D[li].real*T1l[ij].imag + D[li].imag*T1l[ij].real;
            for (k=0;k<nStates;k++) {
              jk = j+k*nStates;
              kl = k+l*nStates;
              rAmp1 = rAmp*G[jk].real - iAmp*G[jk].imag;
              iAmp1 = rAmp*G[jk].imag + iAmp*G[jk].real;
              rAmp2 = rAmp1*T1r[kl].real - iAmp1*T1r[kl].imag;
              iAmp2 = rAmp1*T1r[kl].imag + iAmp1*T1r[kl].real;
              col[nu1[jk]].real += rAmp2;
              col[nu1[jk]].imag += iAmp2;
            }
          }
        }
      }

      mxFree(nu1);
      mxFree(nu2);
      break;
    
    case 3: /* [1, -1] */
      E1left  = (freeL[0]==1) ? Ea : Eb;
      E1right = (freeR[0]==1) ? Ea : Eb;
      E2left  = (freeL[1]==1) ? Ea : Eb;
      E2right = (freeR[1]==1) ? Ea : Eb;
      for (j=0;j<nStates;j++) {
        for (k=0;k<nStates;k++) {
          jk = j+k*nStates;
          for (i=0;i<nStates;i++) {
            ij = i+j*nStates;
            rAmp = T1l[ij].real*G[jk].real - T1l[ij].imag*G[jk].imag;
            iAmp = T1l[ij].real*G[jk].imag + T1l[ij].imag*G[jk].real;
            for (l=0;l<nStates;l++) {
              kl = k+l*nStates;
              li = l+i*nStates;
              il = i+l*nStates;

              idx = binIndex(-(-E2left[i]+E1left[j]-E1right[k]+E2right[l])*dt1,nPoints1);

              rAmp1 = rAmp*T1r[kl].real - iAmp*T1r[kl].imag;
              iAmp1 = rAmp*T1r[kl].imag + iAmp*T1r[kl].real;
              rAmp2 = rAmp1*D[li].real - iAmp1*D[li].imag;
              iAmp2 = rAmp1*D[li].imag + iAmp1*D[li].real;
              
              spec[idx].real += rAmp2;
              spec[idx].imag += iAmp2;
            }
          }
        }
      }
      break;
      
    case 12: /* [1 2 1], e.g. 2D three-pulse ESEEM */
      
      /* Amp = D[ni]*T2l[ij]*T1l[jk]*G[kl]*T1r[lm]*T2r[mn] */
      E1left  = (freeL[0]==1) ? Ea : Eb;
      E2left  = (freeL[1]==1) ? Ea : Eb;
      E3left  = (freeL[2]==1) ? Ea : Eb;
      E1right = (freeR[0]==1) ? Ea : Eb;
      E2right = (freeR[1]==1) ? Ea : Eb;
      E3right = (freeR[2]==1) ? Ea : Eb;
      
      /* bins of nu1 = E1l[k]-E1r[l] + E3l[i]-E3r[n], in the order (k,l,i,n) */
      bins1 = pairSumBins(E1left,E1right,E3left,E3right,nStates,dt1,nPoints1);

      for (k=0;k<nStates;k++) {
        for (l=0;l<nStates;l++) {
          kl = k + l*nStates;
          q = (k*nStates + l)*nStates*nStates;
          for (j=0;j<nStates;j++) {
            jk = j + k*nStates;
            /* Amp0 = T1l[jk]*G[kl] */
            rAmp0 = T1l[jk].real*G[kl].real - T1l[jk].imag*G[kl].imag;
            iAmp0 = T1l[jk].real*G[kl].imag + T1l[jk].imag*G[kl].real;
            for (m=0;m<nStates;m++) {
              lm = l + m*nStates;
              /* Amp1 = Amp0*T1r[lm] */
              rAmp1 = rAmp0*T1r[lm].real - iAmp0*T1r[lm].imag;
              iAmp1 = rAmp0*T1r[lm].imag + iAmp0*T1r[lm].real;
              
              /* nu2 = (E2l[j]-E2r[m]) */
              id2 = binIndex(-(E2left[j]-E2right[m])*dt2,nPoints2);
              
              for (i=0;i<nStates;i++) {
                ij = i + j*nStates;
                /* Amp2 = T2l[ij]*Amp1 */
                rAmp2 = T2l[ij].real*rAmp1 - T2l[ij].imag*iAmp1;
                iAmp2 = T2l[ij].real*iAmp1 + T2l[ij].imag*rAmp1;
                for (n=0;n<nStates;n++) {
                  mn = m + n*nStates;
                  ni = n + i*nStates;
                  /* Amp3 = Amp2*T2r[mn] */
                  rAmp3 = rAmp2*T2r[mn].real - iAmp2*T2r[mn].imag;
                  iAmp3 = rAmp2*T2r[mn].imag + iAmp2*T2r[mn].real;
                  rAmp = D[ni].real*rAmp3 - D[ni].imag*iAmp3;
                  iAmp = D[ni].real*iAmp3 + D[ni].imag*rAmp3;
                  
                  /* nu1 = E1l[k]-E1r[l] + E3l[i]-E3r[n] */
                  idx = bins1[q + i*nStates + n] + id2*nPoints1;
                  
                  spec[idx].real += rAmp;
                  spec[idx].imag += iAmp;
                }
              }
              
            }
          }
        }
      }
      mxFree(bins1);
      break;
    
    case 14: /* [1 1 2], e.g. Hubrich's CF-NF(2) */
      
      /* Amp = D[ni]*T2l[ij]*T1l[jk]*G[kl]*T1r[lm]*T2r[mn] */
      E1left  = (freeL[0]==1) ? Ea : Eb;
      E2left  = (freeL[1]==1) ? Ea : Eb;
      E3left  = (freeL[2]==1) ? Ea : Eb;
      E1right = (freeR[0]==1) ? Ea : Eb;
      E2right = (freeR[1]==1) ? Ea : Eb;
      E3right = (freeR[2]==1) ? Ea : Eb;
      
      /* bins of nu1 = E1l[k]-E1r[l]+E2l[j]-E2r[m], in the order (j,m,k,l)
         of the inner loops below */
      bins1 = pairSumBins(E2left,E2right,E1left,E1right,nStates,dt1,nPoints1);

      /* the outer loops (i,n) fix the buffer column (dimension 2), the
         inner loops (j,m,k,l) stay within it */
      for (i=0;i<nStates;i++) {
        for (n=0;n<nStates;n++) {
          ni = n + i*nStates;
          /* nu2 = (E3l[i]-E3r[n]) */
          id2 = binIndex(-(E3left[i]-E3right[n])*dt2,nPoints2);
          q = 0;
          for (j=0;j<nStates;j++) {
            ij = i + j*nStates;
            /* Amp0 = D[ni]*T2l[ij] */
            rAmp0 = D[ni].real*T2l[ij].real - D[ni].imag*T2l[ij].imag;
            iAmp0 = D[ni].real*T2l[ij].imag + D[ni].imag*T2l[ij].real;
            for (m=0;m<nStates;m++) {
              mn = m + n*nStates;
              /* Amp1 = Amp0*T2r[mn] */
              rAmp1 = rAmp0*T2r[mn].real - iAmp0*T2r[mn].imag;
              iAmp1 = rAmp0*T2r[mn].imag + iAmp0*T2r[mn].real;
              for (k=0;k<nStates;k++) {
                jk = j + k*nStates;
                /* Amp2 = Amp1*T1l[jk] */
                rAmp2 = rAmp1*T1l[jk].real - iAmp1*T1l[jk].imag;
                iAmp2 = rAmp1*T1l[jk].imag + iAmp1*T1l[jk].real;
                for (l=0;l<nStates;l++) {
                  kl = k + l*nStates;
                  lm = l + m*nStates;
                  /* Amp = Amp2*G[kl]*T1r[lm], via Amp5 = Amp2*G[kl] */
                  rAmp5 = rAmp2*G[kl].real - iAmp2*G[kl].imag;
                  iAmp5 = rAmp2*G[kl].imag + iAmp2*G[kl].real;
                  rAmp = rAmp5*T1r[lm].real - iAmp5*T1r[lm].imag;
                  iAmp = rAmp5*T1r[lm].imag + iAmp5*T1r[lm].real;

                  idx = bins1[q++] + id2*nPoints1;
                  
                  spec[idx].real += rAmp;
                  spec[idx].imag += iAmp;
                }
              }
            }
          }
        }
      }
      mxFree(bins1);
      break;

    case 15: /* [1 2 2 1], e.g. 2D-CP */

      /* Amp = D[pi]*T3l[ij]*T2l[jk]*T1l[kl]*G[lm]*T1r[mn]*T2r[no]*T3r[op] */
      E1left  = (freeL[0]==1) ? Ea : Eb;
      E2left  = (freeL[1]==1) ? Ea : Eb;
      E3left  = (freeL[2]==1) ? Ea : Eb;
      E4left  = (freeL[3]==1) ? Ea : Eb;
      
      E1right = (freeR[0]==1) ? Ea : Eb;
      E2right = (freeR[1]==1) ? Ea : Eb;
      E3right = (freeR[2]==1) ? Ea : Eb;
      E4right = (freeR[3]==1) ? Ea : Eb;
      
      /* bins of nu1 = E1l[l]-E1r[m]+E4l[i]-E4r[p], in the order (l,m,i,p)
         of the inner loops below */
      bins1 = pairSumBins(E1left,E1right,E4left,E4right,nStates,dt1,nPoints1);
              
      /* the outer loops (k,n,j,o) fix the buffer column (dimension 2), the
         inner loops (l,m,i,p) stay within it */
      for (k=0;k<nStates;k++) {
        for (n=0;n<nStates;n++) {
          for (j=0;j<nStates;j++) {
            jk = j + k*nStates;
            for (o=0;o<nStates;o++) {
              no = n + o*nStates;
              /* nu2 = E2l[k]-E2r[n] + E3l[j]-E3r[o] */
              id2 = binIndex(-((E2left[k]-E2right[n]) + (E3left[j]-E3right[o]))*dt2,nPoints2);
              /* Amp0 = T2l[jk]*T2r[no] */
              rAmp0 = T2l[jk].real*T2r[no].real - T2l[jk].imag*T2r[no].imag;
              iAmp0 = T2l[jk].real*T2r[no].imag + T2l[jk].imag*T2r[no].real;
              q = 0;
              for (l=0;l<nStates;l++) {
                kl = k + l*nStates;
                /* Amp1 = Amp0*T1l[kl] */
                rAmp1 = rAmp0*T1l[kl].real - iAmp0*T1l[kl].imag;
                iAmp1 = rAmp0*T1l[kl].imag + iAmp0*T1l[kl].real;
                for (m=0;m<nStates;m++) {
                  lm = l + m*nStates;
                  mn = m + n*nStates;
                  /* Amp2 = Amp1*G[lm]*T1r[mn], via Amp5 = Amp1*G[lm] */
                  rAmp5 = rAmp1*G[lm].real - iAmp1*G[lm].imag;
                  iAmp5 = rAmp1*G[lm].imag + iAmp1*G[lm].real;
                  rAmp2 = rAmp5*T1r[mn].real - iAmp5*T1r[mn].imag;
                  iAmp2 = rAmp5*T1r[mn].imag + iAmp5*T1r[mn].real;
                  for (i=0;i<nStates;i++) {
                    ij = i + j*nStates;
                    /* Amp3 = Amp2*T3l[ij] */
                    rAmp3 = rAmp2*T3l[ij].real - iAmp2*T3l[ij].imag;
                    iAmp3 = rAmp2*T3l[ij].imag + iAmp2*T3l[ij].real;
                    for (p=0;p<nStates;p++) {
                      pi = p + i*nStates;
                      op = o + p*nStates;
                      
                      /* Amp4 = Amp3*D[pi], Amp = Amp4*T3r[op] */
                      rAmp4 = rAmp3*D[pi].real - iAmp3*D[pi].imag;
                      iAmp4 = rAmp3*D[pi].imag + iAmp3*D[pi].real;
                      rAmp = rAmp4*T3r[op].real - iAmp4*T3r[op].imag;
                      iAmp = rAmp4*T3r[op].imag + iAmp4*T3r[op].real;
                      
                      idx = bins1[q++] + id2*nPoints1;
                      
                      spec[idx].real += rAmp;
                      spec[idx].imag += iAmp;
                    }
                  }
                }
              }
            }
          }
        }
      }
      mxFree(bins1);
      break;
      
    /*-------------------------------------------------------------------*/
    case 17: /* [1 1 2 2] e.g. 2D refocused primary ESEEM */
      
      /* Amp = D[pi]*T3l[ij]*T2l[jk]*T1l[kl]*G[lm]*T1r[mn]*T2r[no]*T3r[op] */
      E1left  = (freeL[0]==1) ? Ea : Eb;
      E2left  = (freeL[1]==1) ? Ea : Eb;
      E3left  = (freeL[2]==1) ? Ea : Eb;
      E4left  = (freeL[3]==1) ? Ea : Eb;
      E1right = (freeR[0]==1) ? Ea : Eb;
      E2right = (freeR[1]==1) ? Ea : Eb;
      E3right = (freeR[2]==1) ? Ea : Eb;
      E4right = (freeR[3]==1) ? Ea : Eb;
      
      /* bins of nu1 = E1l[l]-E1r[m]+E2l[k]-E2r[n], in the order (k,n,l,m)
         of the inner loops below */
      bins1 = pairSumBins(E2left,E2right,E1left,E1right,nStates,dt1,nPoints1);
              
      /* the outer loops (j,o,i,p) fix the buffer column (dimension 2), the
         inner loops (k,n,l,m) stay within it */
      for (j=0;j<nStates;j++) {
        for (o=0;o<nStates;o++) {
          for (i=0;i<nStates;i++) {
            ij = i + j*nStates;
            for (p=0;p<nStates;p++) {
              op = o + p*nStates;
              pi = p + i*nStates;
              /* nu2 = (E3l[j]-E3r[o]) + (E4l[i]-E4r[p]) */
              id2 = binIndex(-((E3left[j]-E3right[o]) + (E4left[i]-E4right[p]))*dt2,nPoints2);
              /* Amp0 = D[pi]*T3l[ij]*T3r[op], via Amp5 = D[pi]*T3l[ij] */
              rAmp5 = D[pi].real*T3l[ij].real - D[pi].imag*T3l[ij].imag;
              iAmp5 = D[pi].real*T3l[ij].imag + D[pi].imag*T3l[ij].real;
              rAmp0 = rAmp5*T3r[op].real - iAmp5*T3r[op].imag;
              iAmp0 = rAmp5*T3r[op].imag + iAmp5*T3r[op].real;
              q = 0;
              for (k=0;k<nStates;k++) {
                jk = j + k*nStates;
                /* Amp1 = Amp0*T2l[jk] */
                rAmp1 = rAmp0*T2l[jk].real - iAmp0*T2l[jk].imag;
                iAmp1 = rAmp0*T2l[jk].imag + iAmp0*T2l[jk].real;
                for (n=0;n<nStates;n++) {
                  no = n + o*nStates;
                  /* Amp2 = Amp1*T2r[no] */
                  rAmp2 = rAmp1*T2r[no].real - iAmp1*T2r[no].imag;
                  iAmp2 = rAmp1*T2r[no].imag + iAmp1*T2r[no].real;
                  for (l=0;l<nStates;l++) {
                    kl = k + l*nStates;
                    /* Amp3 = Amp2*T1l[kl] */
                    rAmp3 = rAmp2*T1l[kl].real - iAmp2*T1l[kl].imag;
                    iAmp3 = rAmp2*T1l[kl].imag + iAmp2*T1l[kl].real;
                    for (m=0;m<nStates;m++) {
                      lm = l + m*nStates;
                      mn = m + n*nStates;

                      /* Amp = Amp3*G[lm]*T1r[mn], via Amp5 = Amp3*G[lm] */
                      rAmp5 = rAmp3*G[lm].real - iAmp3*G[lm].imag;
                      iAmp5 = rAmp3*G[lm].imag + iAmp3*G[lm].real;
                      rAmp = rAmp5*T1r[mn].real - iAmp5*T1r[mn].imag;
                      iAmp = rAmp5*T1r[mn].imag + iAmp5*T1r[mn].real;

                      idx = bins1[q++] + id2*nPoints1;
                      
                      spec[idx].real += rAmp;
                      spec[idx].imag += iAmp;
                    }
                  }
                }
              }
            }
          }
        }
      }
      mxFree(bins1);
      break;

      default:
      mexErrMsgTxt("Incrementation scheme is currently not supported.");
      break;

  }

}
