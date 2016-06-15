/*
========================================================================
Spectrum = projectzones(Position, Ampl, SegWeights, x)
========================================================================

Calculates weighted function value distribution of scalar function of
one variable. "Position" are the function values, "Ampl" the weights.
"SegWeights" are additional weights, be to multiplied with Ampl.
x specifies the range over which the distribution is to be calculated.

The projection is integrative, i.e. a bin in the output contains the
integral over the distribution between this bin and the next. This
is NOT a sampled distribution. For flat distributions wrt to the
value grid, it doesn't make much difference.

   Position      1xN or Nx1 vector of double
   Ampl          1xN or Nx1 vector of double or a 1x1 double
                 real or complex
   SegWeights    1x(N-1) or (N-1)x1 vector of double
   xData         1xM or Mx1 vector of double
========================================================================
*/

#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[],
         int nrhs, const mxArray *prhs[])
{
  bool anisoAmp, Cplx;
  double meanAmp, rmeanAmp, imeanAmp;
  double delta, left, right, dum;
  double Height, rHeight, iHeight;
  double *Position, *Ampl, *rAmpl, *iAmpl, *SegWeights, *x;
  double *Spectrum, *rSpectrum, *iSpectrum;
  long iSeg, nSeg, idx, first, last, nPoints, nPeaks, nAmpl;

  if (nrhs!=4)
    mexErrMsgTxt("Wrong number of input parameters!");
  if (nlhs!=1)
    mexErrMsgTxt("Wrong number of output parameters!");

  /* get all input parameters */
  Position = mxGetPr(prhs[0]);
  nPeaks = mxGetNumberOfElements(prhs[0]);
  nAmpl = mxGetNumberOfElements(prhs[1]);
  Cplx = mxIsComplex(prhs[1]);
  SegWeights = mxGetPr(prhs[2]);
  nSeg = mxGetNumberOfElements(prhs[2]);
  x = mxGetPr(prhs[3]);
  nPoints = mxGetNumberOfElements(prhs[3]);
  delta = x[1] - x[0];

  if (nSeg!=nPeaks-1)
    mexErrMsgTxt("Wrong number of segment weights. Should be nPeaks-1");
  if ((nAmpl!=1)&&(nAmpl!=nPeaks))
    mexErrMsgTxt("Wrong number of amplitudes! Should be 1 or nPeaks.");
  if (mxIsComplex(prhs[0]))
    mexErrMsgTxt("Peak positions are complex! Only real are allowed.");
  if (mxIsComplex(prhs[2]))
    mexErrMsgTxt("Segment weights are complex! Only real are allowed.");

  /* allocate result array */
  if (Cplx) {
    rAmpl = mxGetPr(prhs[1]);
    iAmpl = mxGetPi(prhs[1]);
    rmeanAmp = rAmpl[0]/delta;
    imeanAmp = iAmpl[0]/delta;
    plhs[0] = mxCreateDoubleMatrix(1,nPoints,mxCOMPLEX);
    rSpectrum = mxGetPr(plhs[0]);
    iSpectrum = mxGetPi(plhs[0]);
  }
  else {
    Ampl = mxGetPr(prhs[1]);
    meanAmp = Ampl[0]/delta;
    plhs[0] = mxCreateDoubleMatrix(1,nPoints,mxREAL);
    Spectrum = mxGetPr(plhs[0]);
  }

  anisoAmp = (nAmpl!=1);

  for (iSeg=0; iSeg<nSeg; iSeg++) {

    left = Position[iSeg];
    right = Position[iSeg+1];
    
    if (mxIsNaN(left) || mxIsNaN(right)) continue;
    
    if (left>right) {dum = left; left = right; right = dum;}
        
    left  = (left -x[0])/delta;
    right = (right-x[0])/delta;
    first = (long)(left);
    last =  (long)(right);

    /*mexPrintf("left %g  right %g  first %d  last %d\n",left, right,first,last);*/

    /* skip segment is completely outside */
    if ((first>=nPoints)||(last<0)) continue;

    if (Cplx) {
      /* get possible separate amplitude for this segment */
      if (anisoAmp) {
        rmeanAmp = (rAmpl[iSeg]+rAmpl[iSeg+1])/2/delta;
        imeanAmp = (iAmpl[iSeg]+iAmpl[iSeg+1])/2/delta;
      }
      if (first==last) {
        rSpectrum[first] += rmeanAmp*SegWeights[iSeg];
        iSpectrum[first] += imeanAmp*SegWeights[iSeg];
      }
      else {
        /* compute height of projection rectangle */
        rHeight = rmeanAmp*SegWeights[iSeg]/(right-left);
        iHeight = imeanAmp*SegWeights[iSeg]/(right-left);
        /* set first element or chop */
        if (first>=0) {
          rSpectrum[first] += rHeight*(first+1-left);
          iSpectrum[first] += iHeight*(first+1-left);
        }
        else first = -1;
        /* set last element or chop */
        if (last<nPoints) {
          rSpectrum[last] += rHeight*(right-last);
          iSpectrum[last] += iHeight*(right-last);
        }
        else last = nPoints;

        /* set possible intermediate elements */
        for (idx=first+1; idx<last; idx++) {
          rSpectrum[idx] += rHeight;
          iSpectrum[idx] += iHeight;
        }
      } /* if (first==last) else */
    }
    else {
      /* get possible separate amplitude for this segment */
      if (anisoAmp)
        meanAmp = (Ampl[iSeg]+Ampl[iSeg+1])/2/delta;
      if (first==last)
        Spectrum[first] += meanAmp*SegWeights[iSeg];
      else {
        /* compute height of projection rectangle */
        Height = meanAmp*SegWeights[iSeg]/(right-left);
        /*mexPrintf("   right-left  %g,   height  %g,  first+1-left %g\n",right-left,Height,first+1-left); */
        /* set first element or chop */
        if (first>=0)
          Spectrum[first] += Height*(first+1-left);
        else
          first = -1;
        /* set last element or chop */
        if (last<nPoints)
          Spectrum[last] += Height*(right-last);
        else
          last = nPoints;
        /* set possible intermediate elements */
        for (idx=first+1; idx<last; idx++)
          Spectrum[idx] += Height;
      } /* if (first==last) else */
    } /* if (Cplx) */
    
  } /* for iSeg */
  
} /* mexFunction */
