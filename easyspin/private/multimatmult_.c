/*
=======================
Code adapted from mmx package on File Exchange
www.mathworks.com/matlabcentral/fileexchange/37515-mmx-multithreaded-matrix-operations-on-n-d-matrices
=======================
*/

#include "mex.h"
#include "string.h"
#include <ctype.h>
#include <math.h>

/*
================
global variables
================
*/
double *A, *B, *C;
int rA, cA, rB ,cB, rC, cC, strideA, strideB, strideC;

/*
==================================================
straightforward implementations of matrix multiply
==================================================
*/
void mulMatMat(double *A, const int rA, const int cA, 
               double *B, const int rB, const int cB,
               double *C) {
  int i, j, k;
  double *a, *c, tmp;  
  for( i=0; i<cB; i++ ){
    c = C + i*rA;
    for( k=0; k<cA; k++ ){
      tmp = B[i*cA+k];
      a = A + k*rA;
  /*        c   = C + i*rA; */
      for( j=0; j<rA; j++ ){
        c[j] += tmp * a[j];
      }
    }
  }
}

/*
=============
mexFunction()
=============
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize  Andim, Bndim, Cndim;
  mwSize *Adims, *Bdims, *Cdims;
  double *Ai, *Bi;
  int i, j, k, N;

  /*
  ==============
  process inputs
  ==============
  */

  if( nrhs != 2 ) {
    mexErrMsgTxt("Only 2 inputs can be accepted.");
  }

  /*
  type check
  */
  if ( (!mxIsDouble(prhs[0])) || (!mxIsDouble(prhs[1])) ) {
    mexErrMsgTxt("Only inputs of type 'double' are supported.");
  }

  /*
  get rA, cA, rB, cB
  */
  A     = mxGetPr(prhs[0]);
  Andim = mxGetNumberOfDimensions(prhs[0]);
  Adims = (mwSize *) mxGetDimensions(prhs[0]);
  rA    = Adims[0];
  cA    = Adims[1];    

  B     = mxGetPr(prhs[1]);
  Bndim = mxGetNumberOfDimensions(prhs[1]);
  Bdims = (mwSize *) mxGetDimensions(prhs[1]);
  rB    = Bdims[0];
  cB    = Bdims[1]; 


  /*
  ================
  dimension checks
  ================
  */
  
  if ( (Andim != Bndim) ) {  /* || (Andim != 3) ) { */
    mexErrMsgTxt("Inputs should have the same number of dimensions.");
  }
/*  
  for( i=0; i<Andim; i++ ) {
    if (Adims[i] != Bdims[i]) {
      mexErrMsgTxt("Inputs should have the same shape.");
    }
  }
*/
  /*
  ===============
  process outputs
  =============== 
  */


  Cndim    = Andim;
  Cndim    = (Cndim > 3) ? Cndim : 3;
  Cdims    = (mwSize *) mxMalloc( Cndim * sizeof(mwSize) );

  /*
  set Cdims[0,1]
  */
  rC = rA;
  cC = cB;
  
  Cdims[0] = rC;
  Cdims[1] = cC;

  for( i=2; i<Cndim; i++ ) {
    Cdims[i] = Adims[i];  /* assumes that A and B are the same size */
  }

  /*
  stride sizes
  */
  strideA    = rA*cA;
  strideB    = rB*cB;
  strideC    = rC*cC;

  /*
  N is the total number of matrix operations
  */
  N  = 1;
  for( i=2; i<Cndim; i++ ) {
    N *= Cdims[i];
  }

  /*
  allocate C
  */
  plhs[0] = mxCreateNumericArray(Cndim, Cdims, mxDOUBLE_CLASS, mxREAL);
  nlhs    = 1;
  C       = mxGetPr(plhs[0]);

  for( i=0; i<N; i++ ){
    /*
    pointers to scheduled data
    */
    Ai = A + strideA*i;
    Bi = B + strideB*i;
    /*
    excecute the task
    */
    mulMatMat(Ai, rA, cA, Bi, rB, cB, C + strideC*i);
  }


  mxFree(Cdims);
}