#include <mex.h>

/*
 * multinucstick(B0,nStates,shifts,startPos,deltaPos,nPoints)
 *
 * shifts: array of size max(2*I+1) x nNuclei
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  int iNuc, nNuclei, nPeaks, i;
  double *nStates, *shifts;
  double p_, centerPos;
  int *iState, M, idx;
  double *position, *spectrum;
  double startPos, deltaPos, nPoints;
  int Debug = 0;
  
  centerPos = mxGetScalar(prhs[0]);
  
  nStates = mxGetPr(prhs[1]);
  nNuclei = mxGetNumberOfElements(prhs[1]);

  /* shifts array: size max(2*I+1) x nNuclei */
  if (mxGetN(prhs[2])!=nNuclei)
    mexErrMsgTxt("Wrong size of shifts array!");
  shifts = mxGetPr(prhs[2]);
  M = mxGetM(prhs[2]);
     
  startPos = mxGetScalar(prhs[3]);
  deltaPos = mxGetScalar(prhs[4]);
  nPoints = mxGetScalar(prhs[5]);

  if (deltaPos<0)
    mexErrMsgTxt("delta cannot be negative.");
    
  centerPos = (centerPos - startPos)/deltaPos;
  for (i=0;i<mxGetNumberOfElements(prhs[2]);i++)
     shifts[i]/=deltaPos;
  
  plhs[0] = mxCreateDoubleMatrix(1,nPoints,mxREAL);
  spectrum = mxGetPr(plhs[0]);
  
  iState = (int*)mxMalloc(nNuclei*sizeof(int));
  for (i=0;i<nNuclei;i++) iState[i]=0;
  
  position = (double*)mxMalloc(nNuclei*sizeof(double));
  for (i=0;i<nNuclei;i++) position[i]=0;
  
  iNuc = 0;
  nPeaks = 0;
  while (iNuc>=0) {
    if (iNuc==nNuclei) { /*if vertex, process and then backtrack */
      idx = position[iNuc-1];
      if ((idx>=0)&&(idx<nPoints)) {
        spectrum[idx]+=1;
        nPeaks++;
      }
      iNuc--;
    }
    if (iState[iNuc]<nStates[iNuc]) { /* if not last state, go deeper, next sibling */
      p_ = (iNuc>0) ? position[iNuc-1] : centerPos;
      position[iNuc] = p_ + shifts[iNuc*M+iState[iNuc]];
      iState[iNuc]++;
      iNuc++;
    }
    else { /* if last state, reset and then backtrack */
      iState[iNuc] = 0;
      iNuc--;
    }
  }
  
  if (Debug)
    mexPrintf("%d peaks\n",nPeaks);
  mxFree(position);
  mxFree(iState);
    
} /* void mexFunction  */
