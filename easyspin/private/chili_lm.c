#include <math.h>
#include <mex.h>

#include "jjj.h"

bool Display = false;

__inline int isodd(int k) { return (k % 2); }
__inline int parity(int k) { return (isodd(k) ? -1 : +1); }
__inline int mini(int a,int b) { return ((a<b) ? a : b); }

int nNuclei;

long blockSize, allocatedSize;

const double sqrt12 = 0.70710678118655; /* sqrt(1.0/2.0); */
const double sqrt13 = 0.57735026918963; /* sqrt(1.0/3.0); */
const double sqrt23 = 0.81649658092773; /* sqrt(2.0/3.0); */

double *ridx, *cidx;
double *MatrixRe, *MatrixIm;

int nRows, nElements, nBasis;

int Lemax, Lomax, Kmax, Mmax;
int jKmin, pSmin, pImax, pIbmax, deltaK;
int MeirovitchSymm;

struct SystemStruct {
  double EZI0, *ReEZI2, *ImEZI2;
  double *d2psi;
  double DirTilt;
  double I, NZI0, HFI0, *ReHFI2, *ImHFI2;
  double Ib, NZI0b, HFI0b, *ReHFI2b, *ImHFI2b;
};

struct DiffusionStruct {
  double* xlk;
  double Rxx, Ryy, Rzz;
  double Exchange;
  int maxL;
};

struct SystemStruct Sys;
struct DiffusionStruct Diff;

#include "chili_lm0.inc" /* functions for S=1/2 and no nuclear spins */
#include "chili_lm1.inc" /* functions for S=1/2 and one nuclear spin */
#include "chili_lm2.inc" /* functions for S=1/2 and two nuclear spins */


/*============================================================================ */
/*============================================================================ */
/*============================================================================ */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  mxArray *T, *D;
  double *R;
  double *basisopts, *allocationOptions;
  int idxS, idxD, idx;
  if (nrhs==1) {
    double *jm = mxGetPr(prhs[0]);
    /*mexPrintf("%0.10f\n",jjj(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5])); */
    return;
  }

  if (nrhs!=4) mexErrMsgTxt("4 input arguments expected.");
  if (nlhs!=4) mexErrMsgTxt("4 output arguments expected.");

  /* Parse spin system input structure */
  if (Display) mexPrintf("Parsing system structure...\n");
  idxS = 0;
  nNuclei = (int)mxGetScalar(mxGetField(prhs[idxS],0,"nNuclei"));
  if (nNuclei>2)
    mexErrMsgTxt("chili_lm0123 works only for 0, 1 or 2 nuclear spins.");  
  Sys.EZI0 = mxGetScalar(mxGetField(prhs[idxS],0,"EZ0"));
  Sys.DirTilt = mxGetScalar(mxGetField(prhs[idxS],0,"DirTilt"));
  Sys.d2psi = mxGetPr(mxGetField(prhs[idxS],0,"d2psi"));
  T = mxGetField(prhs[idxS],0,"EZ2");
  Sys.ReEZI2 = mxGetPr(T);
  Sys.ImEZI2 = mxGetPi(T);
  
  if (nNuclei>=1) {
    Sys.I = mxGetScalar(mxGetField(prhs[idxS],0,"I"));
    Sys.NZI0 = mxGetScalar(mxGetField(prhs[idxS],0,"NZ0"));
    Sys.HFI0 = mxGetScalar(mxGetField(prhs[idxS],0,"HF0"));
    T = mxGetField(prhs[idxS],0,"HF2");
    Sys.ReHFI2 = mxGetPr(T);
    Sys.ImHFI2 = mxGetPi(T);
  }
  if (nNuclei==2) {
    Sys.Ib = mxGetScalar(mxGetField(prhs[idxS],0,"Ib"));
    Sys.NZI0b = mxGetScalar(mxGetField(prhs[idxS],0,"NZ0b"));
    Sys.HFI0b = mxGetScalar(mxGetField(prhs[idxS],0,"HF0b"));
    T = mxGetField(prhs[idxS],0,"HF2b");
    Sys.ReHFI2b = mxGetPr(T);
    Sys.ImHFI2b = mxGetPi(T);
  }
  
  /* Parsing basis set input structure */
  if (Display) mexPrintf("Parsing basis structure...\n");
  basisopts = mxGetPr(prhs[1]);
  
  Lemax = (int)basisopts[0];
  Lomax = (int)basisopts[1];
  Kmax = (int)basisopts[2];
  Mmax = (int)basisopts[3];
  jKmin = (int)basisopts[4];
  pSmin = (int)basisopts[5];
  deltaK = (int)basisopts[6];
  MeirovitchSymm = (int)basisopts[7];
  pImax = (nNuclei>=1) ? (int)basisopts[8] : 0;
  pIbmax = (nNuclei>=2) ? (int)basisopts[9] : 0;
  
  if (Display)
    mexPrintf("  (%d %d %d %d) jKmin %d, pSmin %d, deltaK %d; pImax %d, pIbmax %d\n",
      Lemax,Lomax,Kmax,Mmax,jKmin,pSmin,deltaK,pImax,pIbmax);
  
  /* Parse diffusion input structure */
  if (Display) mexPrintf("Parsing diffusion structure...\n");
  idxD = 2;
  Diff.Exchange = mxGetScalar(mxGetField(prhs[idxD],0,"Exchange"));
  Diff.xlk = mxGetPr(mxGetField(prhs[idxD],0,"xlk"));
  Diff.maxL = (int)mxGetScalar(mxGetField(prhs[idxD],0,"maxL"));
  R = mxGetPr(mxGetField(prhs[idxD],0,"Diff"));
  Diff.Rxx = R[0];
  Diff.Ryy = R[1];
  Diff.Rzz = R[2];

  /* Parse allocation settings */
  allocationOptions = mxGetPr(prhs[3]);
  blockSize = (long)allocationOptions[0];
  if (Display) {
    mexPrintf("  allocation block size: %d\n",blockSize);
  }
  
  /* allocate memory for matrix element indices and values */
  allocatedSize = blockSize;
  ridx = mxCalloc(allocatedSize,sizeof(double));
  cidx = mxCalloc(allocatedSize,sizeof(double));
  MatrixRe = mxCalloc(allocatedSize,sizeof(double));
  MatrixIm = mxCalloc(allocatedSize,sizeof(double));
  if (!(cidx && ridx && MatrixRe && MatrixIm))
    mexErrMsgTxt("Could not allocate initial arrays for Liouvillian.");
  
  /* calculate matrix elements */
  if (Display) mexPrintf("  starting matrix calculation...\n");
  if (nNuclei==0)
    makematrix0();
  else if (nNuclei==1)
    makematrix1();
  else if (nNuclei==2)
    makematrix2();
  if (Display) mexPrintf("  finishing matrix calculation...\n");
  
  /* reallocate (shrink) arrays to correct size */
  allocatedSize = nElements;
  MatrixRe = mxRealloc(MatrixRe,allocatedSize*sizeof(double));
  MatrixIm = mxRealloc(MatrixIm,allocatedSize*sizeof(double));
  ridx = mxRealloc(ridx,allocatedSize*sizeof(double));
  cidx = mxRealloc(cidx,allocatedSize*sizeof(double));
  if (!(cidx && ridx && MatrixRe && MatrixIm))
    mexErrMsgTxt("Could not reallocate arrays for Liouvillian to final size.");
  
  /* adjust indices from 0-based to 1-based */
  for (idx=0;idx<nElements;idx++) {
    ridx[idx]+=1;
    cidx[idx]+=1;
  }
  
  /* allocate and assign mex function output arrays */
  plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
  mxSetPr(plhs[0],ridx);
  mxSetM(plhs[0],allocatedSize);
  mxSetN(plhs[0],1);
  
  plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL);
  mxSetPr(plhs[1],cidx);
  mxSetM(plhs[1],allocatedSize);
  mxSetN(plhs[1],1);
  
  plhs[2] = mxCreateDoubleMatrix(0,0,mxCOMPLEX);
  mxSetPr(plhs[2],MatrixRe);
  mxSetPi(plhs[2],MatrixIm);
  mxSetM(plhs[2],allocatedSize);
  mxSetN(plhs[2],1);
  
  plhs[3] = mxCreateDoubleScalar(nRows);
  
  return;
}
