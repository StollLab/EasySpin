/*
sf_peaks  Bin the peaks of a coherence transfer pathway into a buffer

  sf_peaks(IncSchemeID,bufferRe,bufferIm,dt,idxFreeL,idxFreeR,Ea,Eb,G,D)
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
    bufferRe     buffer, real part: 1 dimension: vector; 2 dimensions:
                 nPoints1 x nPoints2 matrix
    bufferIm     buffer, imaginary part, same size as bufferRe
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
    none; bufferRe and bufferIm are modified in place

  Binning: a peak with frequency nu goes into bin fmod(-nu*dt,1)*nPoints,
  rounded toward zero, with negative bins wrapped around. This rounding is
  odd in nu, so peaks at -nu and +nu land in mirrored bins and conjugate
  peak pairs (nu,a), (-nu,conj(a)) give exactly conjugate-symmetric buffers.
 */

#include "mex.h"
#include <math.h>

/*-------------------------------------------------------------------*/
/* bin of a peak with phase x = -nu*dt (in cycles): fmod(x,1)*nPoints,
   rounded toward zero, negative bins wrapped around; x-trunc(x) is exact
   and equal to fmod(x,1) */
static int binIndex(double x, int nPoints)
{
  int b = (int)((x - (double)(long long)x)*nPoints);
  if (b<0) b += nPoints;
  if (b>=nPoints) b -= nPoints;  /* fractional part rounded up to 1 */
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

/*===================================================================*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int nPeaks, nStates, nPoints1, nPoints2, nFreeEvolutions, nDimensions;
  double *rG, *iG, *rD, *iD, *rT1l, *iT1l, *rT1r, *iT1r;
  double *rT2l, *iT2l, *rT2r, *iT2r, *rT3l, *iT3l, *rT3r, *iT3r;
  double *Ea, *Eb, *E1left, *E1right, *E2left, *E2right, *E3left, *E3right, *E4left, *E4right;
  double rAmp, iAmp, rAmp0,iAmp0, rAmp1, iAmp1, rAmp2, iAmp2, rAmp3, iAmp3, rAmp4, iAmp4, rAmp5, iAmp5;
  double *rSpec, *iSpec, *rCol, *iCol;
  double *freeL, *freeR;
  int id1, id2, idx, q, *nu1, *nu2, *bins1;

  int a;
  mxArray *imagZeros;

  int IncSchemeID;
  int i, j, k, l, ij, jk, kl, li, kj, il;
  int m, n, o, p, lm, mn, no, op, pi, ni;
  double *dt, dt1, dt2;
  
  /*-----------------------------------------------------------------
     Check number of input and output arguments
    ----------------------------------------------------------------- */
  if (nrhs<10)
    mexErrMsgTxt("Insufficient number of input arguments.");
  
  /*-----------------------------------------------------------------
     Read input arguments
    ----------------------------------------------------------------- */

  a = 0;
  /* IncSchemeID ... identifies the incrementation scheme */
  /*    1    [1]         */
  /*    2    [1 1]       */
  /*   11    [1 2]       */
  IncSchemeID = mxGetScalar(prhs[a]);
  if (mxGetNumberOfElements(prhs[a])!=1)
    mexErrMsgTxt("IncSchemeID must be a scalar!");
  switch (IncSchemeID) {
    case  1: nFreeEvolutions = 1; nDimensions = 1; break;
    case  2: nFreeEvolutions = 2; nDimensions = 1; break;
    case  3: nFreeEvolutions = 2; nDimensions = 1; break;
    case 11: nFreeEvolutions = 2; nDimensions = 2; break;
    case 12: nFreeEvolutions = 3; nDimensions = 2; break;
    case 13: nFreeEvolutions = 3; nDimensions = 2; break;
    case 14: nFreeEvolutions = 3; nDimensions = 2; break;
    case 15: nFreeEvolutions = 4; nDimensions = 2; break;
    case 17: nFreeEvolutions = 4; nDimensions = 2; break;
    default: mexErrMsgTxt("Unrecognized incrementation scheme.");
  }
  
  a++;
  /* Spec: spectral storage array */
  rSpec = mxGetPr(prhs[a]);
  if (nDimensions==1) {
    nPoints1 = mxGetM(prhs[a]);
    if (nPoints1==1) nPoints1 = mxGetN(prhs[a]);
    nPoints2 = 1;
  }
  else {
    nPoints1 = mxGetM(prhs[a]);
    nPoints2 = mxGetN(prhs[a]);
  }

  a++;
  /* Spec: spectral storage array */
  iSpec = mxGetPr(prhs[a]);
  if (mxGetNumberOfElements(prhs[a])!=nPoints1*nPoints2)
      mexErrMsgTxt("buffIm has wrong number of elements!");

  
  a++;
  /* dt ... time increment */
  dt = mxGetPr(prhs[a]);
  if (mxGetNumberOfElements(prhs[a])!=nDimensions)
      mexErrMsgTxt("dt has wrong number of elements!");
  dt1 = dt[0];
  if (nDimensions==2) dt2 = dt[1];
    

  a++;
  /* freeL ... index for left-side propagators (1=alpha,2=beta) */
  freeL = mxGetPr(prhs[a]);
  if (mxGetNumberOfElements(prhs[a])!=nFreeEvolutions)
    mexErrMsgTxt("Wrong number of left evolution intervals.");
  
  a++;
  /* freeR ... index for right-side propagators (1=alpha,2=beta) */
  freeR = mxGetPr(prhs[a]);
  if (mxGetNumberOfElements(prhs[a])!=nFreeEvolutions)
    mexErrMsgTxt("Wrong number of right evolution intervals.");

  a++;
  /* Ea ... energies for alpha manifold */
  Ea = mxGetPr(prhs[a]);
  if (mxIsComplex(prhs[a])) mexErrMsgTxt("Ea must be real!");
  nStates = mxGetNumberOfElements(prhs[a]);
  
  a++;
  /* Eb ... energies for beta manifold */
  Eb = mxGetPr(prhs[a]);
  if (mxIsComplex(prhs[a])) mexErrMsgTxt("Eb must be real!");
  if (nStates!=mxGetNumberOfElements(prhs[a]))
    mexErrMsgTxt("Ea and Eb must have the same number of elements.");

  imagZeros = mxCreateDoubleMatrix(nStates,nStates,mxREAL);
  
  a++;
  /* G ... Start coherence matrix */
  if (mxGetNumberOfElements(prhs[a])!=nStates*nStates)
    mexErrMsgTxt("G has wrong size!");
  rG = mxGetPr(prhs[a]);
  iG = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
  
  a++;
  /* D ... detection matrix */
  if (mxGetNumberOfElements(prhs[a])!=nStates*nStates)
    mexErrMsgTxt("D has wrong size!");
  rD = mxGetPr(prhs[a]);
  iD = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
  
  /* get mixing matrices if necessary */
  if (nFreeEvolutions==2) {
    
    if (nrhs!=12) mexErrMsgTxt("Wrong number of input arguments.");
    
    /* T1l, T1r ... mixing matrices */
    a++;
    if (mxGetNumberOfElements(prhs[a])!=nStates*nStates)
      mexErrMsgTxt("M1l has wrong size!");
    rT1l = mxGetPr(prhs[a]);
    iT1l = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
    a++;
    if (mxGetNumberOfElements(prhs[a])!=nStates*nStates)
      mexErrMsgTxt("M1r has wrong size!");
    rT1r = mxGetPr(prhs[a]);
    iT1r = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
  }
  else if (nFreeEvolutions==1) {
    if (nrhs!=10) mexErrMsgTxt("Wrong number of input arguments.");
  }
  else if (nFreeEvolutions==4) {
    if (nrhs!=16) mexErrMsgTxt("Wrong number of input arguments.");

    a++;
    if (mxGetNumberOfElements(prhs[a])!=nStates*nStates)
      mexErrMsgTxt("M1l has wrong size!");
    rT1l = mxGetPr(prhs[a]);
    iT1l = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
    
    a++;
    if (mxGetNumberOfElements(prhs[a])!=nStates*nStates)
      mexErrMsgTxt("M2l has wrong size!");
    rT2l = mxGetPr(prhs[a]);
    iT2l = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
    
    a++;
    if (mxGetNumberOfElements(prhs[a])!=nStates*nStates)
      mexErrMsgTxt("M3l has wrong size!");
    rT3l = mxGetPr(prhs[a]);
    iT3l = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
    
    a++;
    if (mxGetNumberOfElements(prhs[a])!=nStates*nStates)
      mexErrMsgTxt("M1r has wrong size!");
    rT1r = mxGetPr(prhs[a]);
    iT1r = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
    
    a++;
    if (mxGetNumberOfElements(prhs[a])!=nStates*nStates)
      mexErrMsgTxt("M2r has wrong size!");
    rT2r = mxGetPr(prhs[a]);
    iT2r = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
    
    a++;
    if (mxGetNumberOfElements(prhs[a])!=nStates*nStates)
      mexErrMsgTxt("M3r has wrong size!");
    rT3r = mxGetPr(prhs[a]);
    iT3r = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
  }
  else if (nFreeEvolutions==3) {
    if (nrhs!=14) mexErrMsgTxt("Wrong number of input arguments.");

    a++;
    if (mxGetNumberOfElements(prhs[a])!=nStates*nStates)
      mexErrMsgTxt("M1l has wrong size!");
    rT1l = mxGetPr(prhs[a]);
    iT1l = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
    
    a++;
    if (mxGetNumberOfElements(prhs[a])!=nStates*nStates)
      mexErrMsgTxt("M2l has wrong size!");
    rT2l = mxGetPr(prhs[a]);
    iT2l = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
    
    a++;
    if (mxGetNumberOfElements(prhs[a])!=nStates*nStates)
      mexErrMsgTxt("M1r has wrong size!");
    rT1r = mxGetPr(prhs[a]);
    iT1r = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
    
    a++;
    if (mxGetNumberOfElements(prhs[a])!=nStates*nStates)
      mexErrMsgTxt("M2r has wrong size!");
    rT2r = mxGetPr(prhs[a]);
    iT2r = mxIsComplex(prhs[a]) ? mxGetPi(prhs[a]) : mxGetPr(imagZeros);
  }
  else
    mexErrMsgTxt("Unsupported number of free evolution periods.");
  
  nPeaks = 0;
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
          
          /* the - signs here account for the fact the the quantum
           * evolution goes with (-omega_ij)t, but Fourier theory
           * uses (+omega_ij)t */
          idx = binIndex(-(E1left[j]-E1right[k])*dt1,nPoints1);
      
          /* compute amplitude T[kj]*S[jk] */
          rAmp = rD[kj]*rG[jk] - iD[kj]*iG[jk];
          iAmp = rD[kj]*iG[jk] + iD[kj]*rG[jk];
          
          rSpec[idx] += rAmp;
          iSpec[idx] += iAmp;
          nPeaks++;
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
            rAmp = rT1l[ij]*rG[jk] - iT1l[ij]*iG[jk];
            iAmp = rT1l[ij]*iG[jk] + iT1l[ij]*rG[jk];
            for (l=0;l<nStates;l++) {
              kl = k+l*nStates;
              li = l+i*nStates;
              il = i+l*nStates;

              idx = binIndex(-((E2left[i]-E2right[l])+(E1left[j]-E1right[k]))*dt1,nPoints1);

              rAmp1 = rAmp*rT1r[kl] - iAmp*iT1r[kl];
              iAmp1 = rAmp*iT1r[kl] + iAmp*rT1r[kl];
              rAmp2 = rAmp1*rD[li] - iAmp1*iD[li];
              iAmp2 = rAmp1*iD[li] + iAmp1*rD[li];
              
              rSpec[idx] += rAmp2;
              iSpec[idx] += iAmp2;
              nPeaks++;
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

          /* the - signs here account for the fact the the quantum
           * evolution goes with (-omega_ij)t, but Fourier theory
           * uses (+omega_ij)t */
          id1 = binIndex(-(E1left[i]-E1right[j])*dt1,nPoints1);
          id2 = binIndex(-(E2left[i]-E2right[j])*dt2,nPoints2);

          nu1[idx] = id1;
          nu2[idx] = id2;
          if ((id1<0)||(id1>nPoints1-1))
            mexErrMsgTxt("id1 is out of range.");
          if ((id2<0)||(id2>nPoints2-1))
            mexErrMsgTxt("id2 is out of range.");

          idx++;
        }
      }

      /* compute amplitudes and bin peaks */
      /* Amp = Mp[ij]*S[jk]*Mm[kl]*T[li]; */
      /* the outer loops (i,l) fix the buffer column (dimension 2), the
         inner loops (j,k) stay within it */
      nPeaks = 0;
      for (i=0;i<nStates;i++) {
        for (l=0;l<nStates;l++) {
          li = l+i*nStates;
          il = i+l*nStates;
          rCol = rSpec + nu2[il]*nPoints1;
          iCol = iSpec + nu2[il]*nPoints1;
          for (j=0;j<nStates;j++) {
            ij = i+j*nStates;
            rAmp = rD[li]*rT1l[ij] - iD[li]*iT1l[ij];
            iAmp = rD[li]*iT1l[ij] + iD[li]*rT1l[ij];
            for (k=0;k<nStates;k++) {
              jk = j+k*nStates;
              kl = k+l*nStates;
              rAmp1 = rAmp*rG[jk] - iAmp*iG[jk];
              iAmp1 = rAmp*iG[jk] + iAmp*rG[jk];
              rAmp2 = rAmp1*rT1r[kl] - iAmp1*iT1r[kl];
              iAmp2 = rAmp1*iT1r[kl] + iAmp1*rT1r[kl];
              rCol[nu1[jk]] += rAmp2;
              iCol[nu1[jk]] += iAmp2;
              nPeaks++;
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
            rAmp = rT1l[ij]*rG[jk] - iT1l[ij]*iG[jk];
            iAmp = rT1l[ij]*iG[jk] + iT1l[ij]*rG[jk];
            for (l=0;l<nStates;l++) {
              kl = k+l*nStates;
              li = l+i*nStates;
              il = i+l*nStates;

              idx = binIndex(-(-E2left[i]+E1left[j]-E1right[k]+E2right[l])*dt1,nPoints1);

              rAmp1 = rAmp*rT1r[kl] - iAmp*iT1r[kl];
              iAmp1 = rAmp*iT1r[kl] + iAmp*rT1r[kl];
              rAmp2 = rAmp1*rD[li] - iAmp1*iD[li];
              iAmp2 = rAmp1*iD[li] + iAmp1*rD[li];
              
              rSpec[idx] += rAmp2;
              iSpec[idx] += iAmp2;
              nPeaks++;
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
            rAmp0 = rT1l[jk]*rG[kl] - iT1l[jk]*iG[kl];
            iAmp0 = rT1l[jk]*iG[kl] + iT1l[jk]*rG[kl];
            for (m=0;m<nStates;m++) {
              lm = l + m*nStates;
              /* Amp1 = Amp0*T1r[lm] */
              rAmp1 = rAmp0*rT1r[lm] - iAmp0*iT1r[lm];
              iAmp1 = rAmp0*iT1r[lm] + iAmp0*rT1r[lm];
              
              /* nu2 = (E2l[j]-E2r[m]) */
              id2 = binIndex(-(E2left[j]-E2right[m])*dt2,nPoints2);
              
              for (i=0;i<nStates;i++) {
                ij = i + j*nStates;
                /* Amp2 = T2l[ij]*Amp1 */
                rAmp2 = rT2l[ij]*rAmp1 - iT2l[ij]*iAmp1;
                iAmp2 = rT2l[ij]*iAmp1 + iT2l[ij]*rAmp1;
                for (n=0;n<nStates;n++) {
                  mn = m + n*nStates;
                  ni = n + i*nStates;
                  /* Amp3 = Amp2*T2r[mn] */
                  rAmp3 = rAmp2*rT2r[mn] - iAmp2*iT2r[mn];
                  iAmp3 = rAmp2*iT2r[mn] + iAmp2*rT2r[mn];
                  rAmp = rD[ni]*rAmp3 - iD[ni]*iAmp3;
                  iAmp = rD[ni]*iAmp3 + iD[ni]*rAmp3;
                  
                  /* nu1 = E1l[k]-E1r[l] + E3l[i]-E3r[n] */
                  id1 = bins1[q + i*nStates + n];
                  
                  idx = id1 + id2*nPoints1;
                  if ((idx<0)||(idx>nPoints1*nPoints2))
                    mexErrMsgTxt("idx out of range.");
                  
                  rSpec[idx] += rAmp;
                  iSpec[idx] += iAmp;
                  nPeaks++;
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
            rAmp0 = rD[ni]*rT2l[ij] - iD[ni]*iT2l[ij];
            iAmp0 = rD[ni]*iT2l[ij] + iD[ni]*rT2l[ij];
            for (m=0;m<nStates;m++) {
              mn = m + n*nStates;
              /* Amp1 = Amp0*T2r[mn] */
              rAmp1 = rAmp0*rT2r[mn] - iAmp0*iT2r[mn];
              iAmp1 = rAmp0*iT2r[mn] + iAmp0*rT2r[mn];
              for (k=0;k<nStates;k++) {
                jk = j + k*nStates;
                /* Amp2 = Amp1*T1l[jk] */
                rAmp2 = rAmp1*rT1l[jk] - iAmp1*iT1l[jk];
                iAmp2 = rAmp1*iT1l[jk] + iAmp1*rT1l[jk];
                for (l=0;l<nStates;l++) {
                  kl = k + l*nStates;
                  lm = l + m*nStates;
                  /* Amp = Amp2*G[kl]*T1r[lm], via Amp5 = Amp2*G[kl] */
                  rAmp5 = rAmp2*rG[kl] - iAmp2*iG[kl];
                  iAmp5 = rAmp2*iG[kl] + iAmp2*rG[kl];
                  rAmp = rAmp5*rT1r[lm] - iAmp5*iT1r[lm];
                  iAmp = rAmp5*iT1r[lm] + iAmp5*rT1r[lm];

                  idx = bins1[q++] + id2*nPoints1;
                  if ((idx<0)||(idx>nPoints1*nPoints2))
                    mexErrMsgTxt("idx out of range.");
                  
                  rSpec[idx] += rAmp;
                  iSpec[idx] += iAmp;
                  nPeaks++;
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
              rAmp0 = rT2l[jk]*rT2r[no] - iT2l[jk]*iT2r[no];
              iAmp0 = rT2l[jk]*iT2r[no] + iT2l[jk]*rT2r[no];
              q = 0;
              for (l=0;l<nStates;l++) {
                kl = k + l*nStates;
                /* Amp1 = Amp0*T1l[kl] */
                rAmp1 = rAmp0*rT1l[kl] - iAmp0*iT1l[kl];
                iAmp1 = rAmp0*iT1l[kl] + iAmp0*rT1l[kl];
                for (m=0;m<nStates;m++) {
                  lm = l + m*nStates;
                  mn = m + n*nStates;
                  /* Amp2 = Amp1*G[lm]*T1r[mn], via Amp5 = Amp1*G[lm] */
                  rAmp5 = rAmp1*rG[lm] - iAmp1*iG[lm];
                  iAmp5 = rAmp1*iG[lm] + iAmp1*rG[lm];
                  rAmp2 = rAmp5*rT1r[mn] - iAmp5*iT1r[mn];
                  iAmp2 = rAmp5*iT1r[mn] + iAmp5*rT1r[mn];
                  for (i=0;i<nStates;i++) {
                    ij = i + j*nStates;
                    /* Amp3 = Amp2*T3l[ij] */
                    rAmp3 = rAmp2*rT3l[ij] - iAmp2*iT3l[ij];
                    iAmp3 = rAmp2*iT3l[ij] + iAmp2*rT3l[ij];
                    for (p=0;p<nStates;p++) {
                      pi = p + i*nStates;
                      op = o + p*nStates;
                      
                      /* Amp4 = Amp3*D[pi], Amp = Amp4*T3r[op] */
                      rAmp4 = rAmp3*rD[pi] - iAmp3*iD[pi];
                      iAmp4 = rAmp3*iD[pi] + iAmp3*rD[pi];
                      rAmp = rAmp4*rT3r[op] - iAmp4*iT3r[op];
                      iAmp = rAmp4*iT2r[op] + iAmp4*rT3r[op];
                      
                      id1 = bins1[q++];
                      if ((id1<0)||(id1>nPoints1-1))
                        mexErrMsgTxt("id1 is out of range.");
                      if ((id2<0)||(id2>nPoints2-1))
                        mexErrMsgTxt("id2 is out of range.");
                      idx = id1 + id2*nPoints1;
                      
                      rSpec[idx] += rAmp;
                      iSpec[idx] += iAmp;
                      nPeaks++;
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
              rAmp5 = rD[pi]*rT3l[ij] - iD[pi]*iT3l[ij];
              iAmp5 = rD[pi]*iT3l[ij] + iD[pi]*rT3l[ij];
              rAmp0 = rAmp5*rT3r[op] - iAmp5*iT3r[op];
              iAmp0 = rAmp5*iT3r[op] + iAmp5*rT3r[op];
              q = 0;
              for (k=0;k<nStates;k++) {
                jk = j + k*nStates;
                /* Amp1 = Amp0*T2l[jk] */
                rAmp1 = rAmp0*rT2l[jk] - iAmp0*iT2l[jk];
                iAmp1 = rAmp0*iT2l[jk] + iAmp0*rT2l[jk];
                for (n=0;n<nStates;n++) {
                  no = n + o*nStates;
                  /* Amp2 = Amp1*T2r[no] */
                  rAmp2 = rAmp1*rT2r[no] - iAmp1*iT2r[no];
                  iAmp2 = rAmp1*iT2r[no] + iAmp1*rT2r[no];
                  for (l=0;l<nStates;l++) {
                    kl = k + l*nStates;
                    /* Amp3 = Amp2*T1l[kl] */
                    rAmp3 = rAmp2*rT1l[kl] - iAmp2*iT1l[kl];
                    iAmp3 = rAmp2*iT1l[kl] + iAmp2*rT1l[kl];
                    for (m=0;m<nStates;m++) {
                      lm = l + m*nStates;
                      mn = m + n*nStates;

                      /* Amp = Amp3*G[lm]*T1r[mn], via Amp5 = Amp3*G[lm] */
                      rAmp5 = rAmp3*rG[lm] - iAmp3*iG[lm];
                      iAmp5 = rAmp3*iG[lm] + iAmp3*rG[lm];
                      rAmp = rAmp5*rT1r[mn] - iAmp5*iT1r[mn];
                      iAmp = rAmp5*iT1r[mn] + iAmp5*rT1r[mn];

                      idx = bins1[q++] + id2*nPoints1;
                      if ((idx<0)||(idx>nPoints1*nPoints2))
                        mexErrMsgTxt("idx out of range.");
                      
                      rSpec[idx] += rAmp;
                      iSpec[idx] += iAmp;
                      nPeaks++;
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

  
  mxDestroyArray(imagZeros);
  
}
