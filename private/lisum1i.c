/* lisum1c = LInear SUMmation, 1D, programmed in C

  y = lisum1ic(T,PosT,WidT,Pos,Amp,Wid,x)

  Accumulates line shapes into a 1D spectrum by
  interpolatively copying from a shape integral template
  into a spectrum vector. Use an integral as
  shape template, because it is differentiated
  during the binning!

    T: line shape integral template defined over x = 1:length(T)
    PosT: line center in the template
    WidT: line width in the template
    Pos: line positions in spectrum
    Amp: line amplitudes in spectrum
    Wid: line widths in spectrum
    x:  x axis of spectrum

  The x axis of T is assumed to be 1:length(T).
  PosT is the Matlab index of the center. Pos, Wid
  and Amp must contain the same number of elements,
  they are accessed as Pos(:) etc.
  y will be a row vector the same length as x.

 */

#include <math.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  /* input parameters */
  double WidT, PosT;
  double *T, *Pos, *Wid, *Amp, *x;
  /* output parameters */
  double *y;
  /* working variables */
  long nT, nx, nPdat, nIdat, nWdat;
  long idx, first, last, k, idxT;
  double alpha, beta, delta, deltaT;
  double left, right, WidS, Amplitude, pT, remT, Value, ValueOld;
  bool Verbose, InfEncountered, NaNEncountered, NegEncountered;
  /*int Harmonic; */

  Verbose = false;

  /* Argument checks. */
  if (nrhs!=7) mexErrMsgTxt("Wrong number of input parameters!");
  if (nlhs!=1) mexErrMsgTxt("Wrong number of output parameters!");

  /* Get all input arguments. */
  T = mxGetPr(prhs[0]);
  PosT = mxGetScalar(prhs[1]) - 1; /* Matlab -> C array indices */
  WidT = mxGetScalar(prhs[2]);
  Pos = mxGetPr(prhs[3]);
  Amp = mxGetPr(prhs[4]);
  Wid = mxGetPr(prhs[5]);
  x = mxGetPr(prhs[6]);
  /*Harmonic = 0;*/

  /* Get vector lengths. */
  nPdat = mxGetNumberOfElements(prhs[3]);
  nIdat = mxGetNumberOfElements(prhs[4]);
  nWdat = mxGetNumberOfElements(prhs[5]);

  if ((nPdat!=nIdat)||(nPdat!=nWdat))
    mexErrMsgTxt("Pos, Amp, Wid must have the same number of elements!");
  
  nT = mxGetNumberOfElements(prhs[0]);
  nx = mxGetNumberOfElements(prhs[6]);

  /* Allocate result array. */
  plhs[0] = mxCreateDoubleMatrix(1,nx,mxREAL);
  y = mxGetPr(plhs[0]);

  delta = x[1]-x[0];
  deltaT = 1;
  beta = delta/deltaT;

  /*if (Verbose) mexPrintf("begin of lisum1ic loop\n"); */
  /* Loop over all peaks to add. */
  InfEncountered = false;
  NaNEncountered = false;
  NegEncountered = false;
  for (k=0; k<nPdat; k++) {

    /* Skip if NaN or inf. */
    if (mxIsNaN(Pos[k])||mxIsNaN(Amp[k])||mxIsNaN(Wid[k])) {NaNEncountered=true; continue;}
    if (mxIsInf(Pos[k])||mxIsInf(Amp[k])||mxIsInf(Wid[k])) {InfEncountered=true; continue;}
    
    /* Set to zero if negative width */
    WidS = Wid[k];
    if (WidS<0) {NegEncountered=true; WidS=0;}
    
    /* If width is zero, set it to a small value */
    if (WidS<delta/100) WidS = delta/100;

    /* Scaling factor between template and spectrum abscissae */
    alpha = WidT/WidS * beta;
    
    /* Left and right border of line in spectrum. */
    left = (Pos[k]-x[0])/delta - PosT/alpha;
    left -= 0.5; /* Integration */
    right = left + (nT-1)/alpha;

    /* Skip peak if completely outside. */
    if ((left>=nx-1)||(right<=0)) continue; 

    /* Chop if necessary. */
    first = (long)((left<0) ? 0 : ceil(left));
    last = (long)((right>nx-1) ? nx-1 : ceil(right));

    /* scaled amplitude */
    Amplitude = Amp[k]/beta;
    /*if (Harmonic>0) { */
    /*  Amplitude *= WidT/WidS; */
    /*  if (Harmonic>1) Amplitude *= WidT/WidS; */
    /*} */
    /*if (Verbose) */
    /*  mexPrintf("  setup> left %g  right %g   first %d  last %d    A %g\n",left,right,first,last,Amplitude); */

    /* Loop over all points in spectrum to update. */
    pT = alpha*(first-left);
    idxT = (long)(pT-alpha);
    remT = pT - alpha - idxT;
    ValueOld = (idxT<0) ? T[0] : (T[idxT]*(1-remT) + T[idxT+1]*remT);
    /*if (Verbose) mexPrintf("    start: pT %5f    %11e ...\n",pT-alpha,ValueOld); */
    
    for (idx=first; idx<last; idx++, pT+=alpha, ValueOld=Value) {
      idxT = (long)pT;
      remT = pT - idxT;
      Value = T[idxT]*(1-remT) + T[idxT+1]*remT;
      /*if (Verbose) */
      /*  mexPrintf("    idxS %4d  pT %5f    %11e - %-11e...",idx,pT,Value,ValueOld); */
      y[idx] += Amplitude*(Value-ValueOld);
      /*if (Verbose) mexPrintf("added\n"); */
    }
    idxT = (long)pT;
    remT = pT - idxT;
    Value = idxT<nT-2 ? T[idxT]*(1-remT) + T[idxT+1]*remT : T[nT-1];
    /*if (Verbose) */
    /*    mexPrintf("  L:idxS %4d  pT %5f    %11e - %-11e...",idx,pT,Value,ValueOld); */
    y[idx] += Amplitude*(Value-ValueOld);
    /*if (Verbose) mexPrintf("added\n"); */
  } /* loop over all lines */
  /*if (Verbose) mexPrintf("end of lisum1ic loop\n"); */
  if (InfEncountered) mexPrintf("********** lisum1ic: Inf encountered!! Please report! **************\n");
  if (NegEncountered) mexPrintf("********** lisum1ic: Negative width encountered!! Please report! **************\n");
} /* void mexFunction */
