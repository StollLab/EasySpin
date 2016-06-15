/*
========================================================================
Spectrum = projecttriangles(Tri,Areas,Fun,Amplitude,x)
========================================================================

  Given data Fun with a triangulation Tri and triangle
  areas Areas, projecttriangles computes the function value
  distribution over the range given in x, including additional
  weights from Amplitude for each data point.

  Tri:       triangulation, 3xN uint32 array
  Areas:     1xN or Nx1 array of areas for the triangles
  Fun:       1xM or Mx1 array of position data
  Amplitude: 1xM or Mx1 or 1x1 array of amplitude data
  xData:     1xK array, x axis vector for spectrum

========================================================================
*/

/* guarantee portability between 32 and 64-bit systems */
#ifdef bit64
typedef int INT32;
#else
typedef int INT32;
#endif

#include <mex.h>
#include "math.h"

void triproject(double spec[], INT32 tri[], double wei[],
         double fun[], double amp[],
         double x[], INT32 nPoints,
         INT32 nTri, INT32 nLines, INT32 isoInt)
{
  double Amplitude, left, middle, right, delta, dum;
  double Area, Width1, Width2, Width, f0;
  INT32 iTri, idx1, idx2, idx3, first1, last1, first2, last2, idx;

  Amplitude = amp[0];
  delta = x[1] - x[0];
  for (iTri=0; iTri<nTri; iTri++) {

    /* retrieve indices from triangulation */
    idx1 = tri[3*iTri]-1;
    idx2 = tri[3*iTri+1]-1;
    idx3 = tri[3*iTri+2]-1;

    /* get position data */
    left = fun[idx1];
    middle = fun[idx2];
    right = fun[idx3];
    
    Area = wei[iTri];
    
    /* skip triangles with NaN positions */
    if (mxIsNaN(left) || mxIsNaN(middle) || mxIsNaN(right)) continue;
    
    /* sort and scale positions */
    if (right<left) {dum=left; left=right; right=dum;}
    if (middle<left) {dum=left; left=middle; middle=dum;}
    else
      if (middle>right) {dum=right; right=middle; middle=dum;}
    left   = (left  -x[0])/delta;
    middle = (middle-x[0])/delta;
    right  = (right -x[0])/delta;
    
    Width = right-left;
    Width1 = middle-left;
    Width2 = right-middle;

    /* compute mean amplitude if necessary */
    if (!isoInt) Amplitude = (amp[idx1]+amp[idx2]+amp[idx3])/3;

    first1 = (INT32)left;
    last1 = (INT32)middle;
    first2 = last1;
    last2 = (INT32)right;

    /* if left triangle has non-zero width*/
    if (Width1>0)
      /* if left triangle lies at least partially within range */
      if ((first1<nPoints)&&(last1>=0)) {
        /* compute prefactor (Height = 2*Amplitude*Area/Width) */
        f0 = (2*Amplitude*Area/Width)/Width1/delta;
        /* Integral over spectrum (mind x axis!) should be equal to sum of weights. */
        if (first1==last1)
          spec[first1] += f0*Width1*Width1/2;
        else {
          /* update or chop first bin */
          if (first1>=0)
            spec[first1] += f0*(first1+1-left)/2*(first1+1-left);
          else first1 = -1;
          /* update or chop last bin */
          if (last1<nPoints)
            spec[last1] += f0*((last1+middle)/2-left)*(middle-last1);
          else last1 = nPoints;
          /* update intermediate bins */
          left -= 0.5;
          for (idx=first1+1; idx<last1; idx++)
            spec[idx] += f0*(idx-left);
        }
      }

    /* if right triangle has non-zero width */
    if (Width2>0)
      /* if right triangle lies at least partially within range */
      if ((first2<nPoints)&&(last2>=0)) {
        /* compute prefactor (Height = 2*Amplitude*Area/Width) */
        f0 = (2*Amplitude*Area/Width)/Width2/delta;
        /* Integral over spectrum (mind x axis!) should be equal to sum of weights. */
        if (first2==last2)
          spec[first2] += f0*Width2*Width2/2;
        else {
          /* update or chop first bin */
          if (first2>=0)
            spec[first2] += f0*(first2+1-middle)*(right-(middle+first2+1)/2);
          else first2 = -1;
          /* update or chop last bin */
          if (last2<nPoints)
            spec[last2] += f0*(right-last2)*(right-last2)/2;
          else last2 = nPoints;
          /* update intermediate bins */
          right -= 0.5;
          for (idx=first2+1; idx<last2; idx++)
            spec[idx] += f0*(right-idx);
        }
      }
    
    /* Both left and right are zero-width "triangles" */
    if ((Width1==0)&&(Width2==0)) {
      first1 = (INT32)left;
      if ((first1>=0)&&(first1<nPoints))
        spec[first1] += Amplitude*Area/delta;
    }
    
  } /* for iTri */
}


void mexFunction(int nlhs, mxArray *plhs[],
         int nrhs, const mxArray *prhs[])
{
  INT32 *Triangulation;
  double *TriAreas, *Function, *Amplitudes, *x;

  INT32 nTri, isoIntensity, nPoints, nLines, nAmps;

  double *Spectrum;

  if (nrhs!=5)
    mexErrMsgTxt("Wrong number of input arguments!");
  if (nlhs!=1)
    mexErrMsgTxt("Wrong number of output arguments!");

  /* get all input parameters */
  if (mxGetM(prhs[0])!=3)
    mexErrMsgTxt("Triangulation array must have 3 rows!");
  if (!mxIsUint32(prhs[0]))
    mexErrMsgTxt("Triangulation array must be uint32!");    
  Triangulation = (INT32*)mxGetData(prhs[0]);
  nTri = (INT32)mxGetN(prhs[0]);
  if (mxGetNumberOfElements(prhs[1])!=nTri)
    mexErrMsgTxt("Number of triangles and areas must match!");
  TriAreas = mxGetPr(prhs[1]);
  Function = mxGetPr(prhs[2]);
  nLines = mxGetNumberOfElements(prhs[2]);

  nAmps = mxGetNumberOfElements(prhs[3]);
  isoIntensity = (nAmps==1);
  if (!isoIntensity)
    if (nAmps!=nLines)
      mexErrMsgTxt("Number of line positions and amplitudes must match!");
  Amplitudes = mxGetPr(prhs[3]);
  
  x = mxGetPr(prhs[4]);
  nPoints = (INT32)(mxGetNumberOfElements(prhs[4]));

  /* allocate result array */
  plhs[0] = mxCreateDoubleMatrix(1,nPoints,mxREAL);
  Spectrum = mxGetPr(plhs[0]);

  /* call computing function */
  triproject(Spectrum,Triangulation,TriAreas,Function,Amplitudes,
             x,nPoints,nTri,nLines,isoIntensity);
  
}
