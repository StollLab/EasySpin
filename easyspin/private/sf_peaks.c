/*
sf_peaks(IncSchemeID,bufferRe,bufferIm,dt,idxFreeL,idxFreeR,Ea,Eb,G,D,M1p,M1m)
 */

#include "mex.h"
#include <math.h>

/*===================================================================*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int nPeaks, nStates, nPoints1, nPoints2, nFreeEvolutions, nDimensions;
  double *rG, *iG, *rD, *iD, *rT1l, *iT1l, *rT1r, *iT1r;
  double *rT2l, *iT2l, *rT2r, *iT2r, *rT3l, *iT3l, *rT3r, *iT3r;
  double *Ea, *Eb, *E1left, *E1right, *E2left, *E2right, *E3left, *E3right, *E4left, *E4right;
  double rAmp, iAmp, rAmp0,iAmp0, rAmp1, iAmp1, rAmp2, iAmp2, rAmp3, iAmp3, rAmp4, iAmp4, rAmp5, iAmp5;
  double *rSpec, *iSpec;
  double *freeL, *freeR;
  double id1, id2, nu, nu1_, nu2_;
  int idx, *nu1, *nu2;

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
          nu = E1left[j]-E1right[k];
          id1 = fmod(-nu*dt1,1.0)*nPoints1;
          if (id1>=0) id1 = floor(id1);
          else { id1 = ceil(id1); if (id1<0) id1 += nPoints1; }
          idx = id1;
      
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

              nu = (E2left[i]-E2right[l])+(E1left[j]-E1right[k]);

              id1 = fmod(-nu*dt1,1.0)*nPoints1;
              if (id1>=0) id1 = floor(id1);
              else { id1 = ceil(id1); if (id1<0) id1 += nPoints1; }
              idx = id1;

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
          id1 = fmod(-(E1left[i]-E1right[j])*dt1,1.0)*nPoints1;
          id2 = fmod(-(E2left[i]-E2right[j])*dt2,1.0)*nPoints2;

          if (id1>=0) id1 = floor(id1);
          else { id1 = ceil(id1); if (id1<0) id1 += nPoints1; }
          if (id2>=0) id2 = floor(id2);
          else { id2 = ceil(id2); if (id2<0) id2 += nPoints2; }

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
      nPeaks = 0;
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
              rAmp1 = rAmp*rT1r[kl] - iAmp*iT1r[kl];
              iAmp1 = rAmp*iT1r[kl] + iAmp*rT1r[kl];
              rAmp2 = rAmp1*rD[li] - iAmp1*iD[li];
              iAmp2 = rAmp1*iD[li] + iAmp1*rD[li];
              idx = nu1[jk] + nu2[il]*nPoints1;
              rSpec[idx] += rAmp2;
              iSpec[idx] += iAmp2;
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

              nu = -E2left[i]+E1left[j]-E1right[k]+E2right[l];

              id1 = fmod(-nu*dt1,1.0)*nPoints1;
              if (id1>=0) id1 = floor(id1);
              else { id1 = ceil(id1); if (id1<0) id1 += nPoints1; }
              idx = id1;

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
      
      for (k=0;k<nStates;k++) {
        for (l=0;l<nStates;l++) {
          kl = k + l*nStates;
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
              nu2_ = (E2left[j]-E2right[m]);
              id2 = fmod(-nu2_*dt2,1.0)*nPoints2;
              if (id2>=0) id2 = floor(id2);
              else { id2 = ceil(id2); if (id2<0) id2 += nPoints2; }
              
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
                  nu1_ = (E1left[k]-E1right[l]) + (E3left[i]-E3right[n]);
                  id1 = fmod(-nu1_*dt1,1.0)*nPoints1;
                  if (id1>=0) id1 = floor(id1);
                  else { id1 = ceil(id1); if (id1<0) id1 += nPoints1; }
                  
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
      break;
    
    case 14: /* [1 1 2], e.g. Hubrich's CF-NF(2) */
      
      /* Amp = D[ni]*T2l[ij]*T1l[jk]*G[kl]*T1r[lm]*T2r[mn] */
      E1left  = (freeL[0]==1) ? Ea : Eb;
      E2left  = (freeL[1]==1) ? Ea : Eb;
      E3left  = (freeL[2]==1) ? Ea : Eb;
      E1right = (freeR[0]==1) ? Ea : Eb;
      E2right = (freeR[1]==1) ? Ea : Eb;      
      E3right = (freeR[2]==1) ? Ea : Eb;      
      
      for (k=0;k<nStates;k++) {
        for (l=0;l<nStates;l++) {
          kl = k + l*nStates;
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
              /* nu1 = E1l[k]-E1r[l]+E2l[j]-E2r[m] */
              nu1_ = (E1left[k]-E1right[l]) + (E2left[j]-E2right[m]);
              id1 = fmod(-nu1_*dt1,1.0)*nPoints1;
              if (id1>=0) id1 = floor(id1);
              else { id1 = ceil(id1); if (id1<0) id1 += nPoints1; }
              
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
                  
                  /* nu2 = (E3l[i]-E3r[n]) */
                  nu2_ = (E3left[i]-E3right[n]);
                  id2 = fmod(-nu2_*dt2,1.0)*nPoints2;
                  if (id2>=0) id2 = floor(id2);
                  else { id2 = ceil(id2); if (id2<0) id2 += nPoints2; }
                  
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
      
      for (l=0;l<nStates;l++) {
        for (m=0;m<nStates;m++) {
          lm = l + m*nStates;
          for (i=0;i<nStates;i++) {
            for (p=0;p<nStates;p++) {
              pi = p + i*nStates;
              /* nu1 = E1l[l]-E1r[m]+E4l[i]-E4r[p] */
              nu1_ = (E1left[l]-E1right[m]) + (E4left[i]-E4right[p]);
              id1 = fmod(-nu1_*dt1,1.0)*nPoints1;
              if (id1>=0) id1 = floor(id1);
              else { id1 = ceil(id1); if (id1<0) id1 += nPoints1; }
              
              for (k=0;k<nStates;k++) {
                kl = k + l*nStates;
                /* Amp0 = T1l[kl]*G[lm] */
                rAmp0 = rT1l[kl]*rG[lm] - iT1l[kl]*iG[lm];
                iAmp0 = rT1l[kl]*iG[lm] + iT1l[kl]*rG[lm];     
                for (n=0;n<nStates;n++) {
                  mn = m + n*nStates;
                  /* Amp1 = Amp0*T1r[mn] */
                  rAmp1 = rAmp0*rT1r[mn] - iAmp0*iT1r[mn];
                  iAmp1 = rAmp0*iT1r[mn] + iAmp0*rT1r[mn];
                  for (j=0;j<nStates;j++) {
                    ij = i + j*nStates;
                    jk = j + k*nStates;
                    /* Amp2 = D[pi]*T3l[ij]*T2l[jk]*Amp1 */
                    rAmp2 = rT2l[jk]*rAmp1 - iT2l[jk]*iAmp1;
                    iAmp2 = rT2l[jk]*iAmp1 + iT2l[jk]*rAmp1;
                    rAmp3 = rT3l[ij]*rAmp2 - iT3l[ij]*iAmp2;
                    iAmp3 = rT3l[ij]*iAmp2 + iT3l[ij]*rAmp2;
                    rAmp2 = rD[pi]*rAmp3 - iD[pi]*iAmp3;
                    iAmp2 = rD[pi]*iAmp3 + iD[pi]*rAmp3;
                    for (o=0;o<nStates;o++) {
                      no = n + o*nStates;
                      op = o + p*nStates;
                      
                      /* Amp = Amp2*T2r[no]*T3r[op] */
                      rAmp3 = rAmp2*rT2r[no] - iAmp2*iT2r[no];
                      iAmp3 = rAmp2*iT2r[no] + iAmp2*rT2r[no];
                      rAmp = rAmp3*rT3r[op] - iAmp3*iT3r[op];
                      iAmp = rAmp3*iT2r[op] + iAmp3*rT3r[op];
                      
                      /* nu2 = E2l[k]-E2r[n] + E3l[j]-E3r[o] */
                      nu2_ = (E2left[k]-E2right[n]) + (E3left[j]-E3right[o]);
                      id2 = fmod(-nu2_*dt2,1.0)*nPoints2;
                      if (id2>=0) id2 = floor(id2);
                      else { id2 = ceil(id2); if (id2<0) id2 += nPoints2; }
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
      
      for (l=0;l<nStates;l++) {
        for (m=0;m<nStates;m++) {
          lm = l + m*nStates;
          for (k=0;k<nStates;k++) {
            kl = k + l*nStates;
            /* Amp0 = T1l[kl]*G[lm] */
            rAmp0 = rT1l[kl]*rG[lm] - iT1l[kl]*iG[lm];
            iAmp0 = rT1l[kl]*iG[lm] + iT1l[kl]*rG[lm];
            for (n=0;n<nStates;n++) {
              mn = m + n*nStates;
              /* nu1 = E1l[l]-E1r[m]+E2l[k]-E2r[n] */
              nu1_ = (E1left[l]-E1right[m]) + (E2left[k]-E2right[n]);
              id1 = fmod(-nu1_*dt1,1.0)*nPoints1;
              if (id1>=0) id1 = floor(id1);
              else { id1 = ceil(id1); if (id1<0) id1 += nPoints1; }
              /* Amp1 = Amp0*T1r[mn] */
              rAmp1 = rAmp0*rT1r[mn] - iAmp0*iT1r[mn];
              iAmp1 = rAmp0*iT1r[mn] + iAmp0*rT1r[mn];
              
              for (j=0;j<nStates;j++) {
                jk = j + k*nStates;
                /* Amp2 = T2l[jk]*Amp1 */
                rAmp2 = rT2l[jk]*rAmp1 - iT2l[jk]*iAmp1;
                iAmp2 = rT2l[jk]*iAmp1 + iT2l[jk]*rAmp1;
                for (o=0;o<nStates;o++) {
                  no = n + o*nStates;
                  /* Amp3 = Amp2*T2r[no] */
                  rAmp3 = rAmp2*rT2r[no] - iAmp2*iT2r[no];
                  iAmp3 = rAmp2*iT2r[no] + iAmp2*rT2r[no];
                  for (i=0;i<nStates;i++) {
                    ij = i + j*nStates;
                    /* Amp4 = T3l[ij]*Amp3 */
                    rAmp4 = rT3l[ij]*rAmp3 - iT3l[ij]*iAmp3;
                    iAmp4 = rT3l[ij]*iAmp3 + iT3l[ij]*rAmp3;
                    for (p=0;p<nStates;p++) {
                      op = o + p*nStates;
                      pi = p + i*nStates;
                      
                      /* Amp = D[pi*Amp4*T3r[op] */
                      rAmp5 = rAmp4*rT3r[op] - iAmp4*iT3r[op];
                      iAmp5 = rAmp4*iT3r[op] + iAmp4*rT3r[op];
                      rAmp = rD[pi]*rAmp5 - iD[pi]*iAmp5;
                      iAmp = rD[pi]*iAmp5 + iD[pi]*rAmp5;
                      
                      /* nu2 = (E3l[k]-E3r[n]) + (E4l[j]-E4r[o]) */
                      nu2_ = (E3left[j]-E3right[o]) + (E4left[i]-E4right[p]);
                      id2 = fmod(-nu2_*dt2,1.0)*nPoints2;
                      if (id2>=0) id2 = floor(id2);
                      else { id2 = ceil(id2); if (id2<0) id2 += nPoints2; }
              
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
        }
      }
      break;

      default:
      mexErrMsgTxt("Incrementation scheme is currently not supported.");
      break;      

  }

  
  mxDestroyArray(imagZeros);
  
}
