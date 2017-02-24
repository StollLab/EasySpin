// =======================
// Code adapted from mmx package on File Exchange
// http://www.mathworks.com/matlabcentral/fileexchange/37515-mmx-multithreaded-matrix-operations-on-n-d-matrices
// =======================

#include "mex.h"
#include "string.h"
#include <ctype.h>
#include <math.h>

// function declarations
#include "matrix_fun_simple.c"


// ================
// global variables
// ================
double *A, *B, *C;
int rA, cA, rB ,cB, rC, cC, strideA, strideB, strideC;

// =============
// mexFunction()
// =============
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize  Andim, Bndim, Cndim;
  mwSize *Adims, *Bdims, *Cdims;
  double *Ai, *Bi;
  int i, j, k, N;

  // no input: print documentation
  if( nrhs==0 ) {
    mexPrintf(  "=======================================================\n"
          "===  mmx(): fast, multithreaded, n-D multiplication ===\n"
          "Basic usage:\nThe command   C = mmx('mult',A,B);\n"
          "is equivalent to the matlab loop\n"
          " for i=1:N\n    C(:,:,i) = A(:,:,i)*B(:,:,i);\n end\n"
          "===== Type 'help mmx' for detailed information. =======\n"
          "=======================================================\n");   
  }

  // just getting help or setting the thread count, exit now
  if( nrhs < 2 ) {
    return;   
  }

  // ==============
  // process inputs
  // ==============   

  if( nrhs != 2 ) {
    mexErrMsgTxt("Only 2 inputs can be accepted.");
  }


  // type check
  if ( (!mxIsDouble(prhs[0])) || (!mxIsDouble(prhs[1])) ) {
    mexErrMsgTxt("Only inputs of type 'double' are supported.");
  }

  // get rA, cA, rB, cB
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


  // ================
  // dimension checks
  // ================  
  
  if ( (Andim != Bndim) ) { //|| (Andim != 3) ) {
    mexErrMsgTxt("Inputs should have the same number of dimensions.");
  }
  
  for( i=0; i<Andim; i++ ) {
    if (Adims[i] != Bdims[i]) {
      mexErrMsgTxt("Inputs should have the same shape.");
    }
  }

  // ===============
  // process outputs
  // =============== 


  Cndim    = Andim;
  Cndim    = (Cndim > 3) ? Cndim : 3;
  Cdims    = (mwSize *) mxMalloc( Cndim * sizeof(mwSize) );

  // set Cdims[0,1]
  rC = rA;
  cC = cB;

  for( i=0; i<Cndim; i++ ) {
    Cdims[i] = Adims[i];
  }

  // stride sizes
  strideA    = rA*cA;
  strideB    = rB*cB;
  strideC    = rC*cC;

  // N is the total number of matrix operations
  N  = 1;
  for( i=2; i<Cndim; i++ ) {
    N *= Cdims[i];
  }

  // if one of the output dimensions is 0 we're done, goodbye
  if ( Cdims[0]*Cdims[1]*N == 0 ) {
    return;
  }


  // allocate C
  plhs[0] = mxCreateNumericArray(Cndim, Cdims, mxDOUBLE_CLASS, mxREAL);
  nlhs    = 1;
  C       = mxGetPr(plhs[0]);

  for( i=0; i<N; i++ ){
    // pointers to scheduled data
    Ai = A + strideA*i;
    Bi = B + strideB*i;
    // excecute the task
    mulMatMat(Ai, rA, cA, Bi, rB, cB, C + strideC*i);
  }


  mxFree(Cdims);
}

