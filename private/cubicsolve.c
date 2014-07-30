/*
cubicsolvec.c      solve cubic equation in interval [0 1]

  [r,d1,d2] = cubicsolvec(C,MultipleRoots)

  Returns the roots in the interval [0,1] of the
  cubic polynomial with real coefficients in the
  4-element vector C. Polynomial:

     y = C(1)*t^3 + C(2)*t^2 + C(3)*t + C(4)

  If MultipleRoots is set to false, at most one root
  will be returned.

  r contains the roots (up to 3). d1 and d2 are the
  first and second derivative at t=r.

This is an EasySpin function.
Author: Stefan Stoll
 */

#include <math.h>
#include <mex.h>

const double PI = 3.141592653589793238462643383279502884;
const double eps = 1e-11;
const int maxNewtonSteps = 11;

int exactSolver(const double C1, const double C2, const double C3,
                const double C4, double *Roots)
/* Computes all real roots r with 0<=r<=1 of the cubic polynomial
 y = C1*t^3 + C2*t^2 + C3*t + C4 using exact formulas. It assumes
 that at least one coefficient is non-zero. */
{
  double a, b, c, Q, R, A, B, Q2, QQQ, theta;
  double s1, s2, s3;
  int N;

  N = 0;

  if (C1!=0) {
    /* Cubic equation: y = t^3 + a*t^2 + b*t + c   */
    a = C2/C1;
    b = C3/C1;
    c = C4/C1;

    Q = (a*a-3*b)/9;
    R = (2*a*a*a-9*a*b+27*c)/54;
    QQQ = Q*Q*Q;

    if (R*R<QQQ) {
      /* three real roots  */
      theta = acos(R/sqrt(QQQ));
      Q2 = -2*sqrt(Q);

      s1 = Q2*cos(theta/3) - a/3;
      s2 = Q2*cos((theta+2*PI)/3) - a/3;
      s3 = Q2*cos((theta-2*PI)/3) - a/3;

      if (s1>=0 && s1<=1) Roots[N++] = s1;
      if (s2>=0 && s2<=1) Roots[N++] = s2;
      if (s3>=0 && s3<=1) Roots[N++] = s3;
    }
    else {
      /* one real root   */
      A = -pow(fabs(R)+sqrt(R*R-QQQ),1/3.0);
      if (R<0) A = -A;
      B = A==0 ? 0 : Q/A;

      s1 = A + B - a/3;
      if (s1>=0 && s1<=1) Roots[N++] = s1;
    }
  }
  else
  if (C2!=0) {
    /* Quadratic equation: y = C2*t^2 + C3*t + C4    */
    R = C3*C3 - 4*C2*C4; /* discriminant    */
    if (R>=0) { /* negative -> 2 complex solutions   */
      Q = -(C3 + (C3<0 ? -1 : 1)*sqrt(R))/2;

      s1 = Q/C2;
      s2 = C4/Q;

      if (s1>=0 && s1<=1) Roots[N++] = s1;
      if (s2>=0 && s2<=1) Roots[N++] = s2;
    }
  }
  else
  if (C3!=0) {
    /* Linear equation: y = C3*t + C4   */
    s1 = -C4/C3;
    if (s1>=0 && s1<=1) Roots[N++] = s1;
  }
  else
  if (C4!=0) {
    /* Constant, no solution     */
  }
  else {
    /* Zero, infinitely many solutions             */
    mexErrMsgTxt("All coefficients of the cubic spline are zero!");
  }

  return N;

}

int cubicsolve(const double a, const double b, const double c,
               const double d, bool SingleRoot, double *Solution)
/* SingleRoot = true: returns either 1 or 0 roots Solution
   SingleRoot = false: returns 0-3 roots in Solution
*/
{
  bool isMonotonic = false, SameSign;
  double yex;
  int iStep;
  double delta, t;
  int i, nSolutions = 0;

  /* For some cases it's easy to detect monotonicity */
  /* ------------------------------------------------------- */
  if (SingleRoot && c*(3.0*a+2.0*b+c)>0) {
    if (c*a<0)
      isMonotonic = true;
    /*
    else {
       The following test cost too much...
      p = (a!=0) ? p = -b/a/3.0 : -100;
      Discriminant = b*b-3.0*a*c;
      isMonotonic = (p<=0) || (p>=1) || (Discriminant<=0);
    }*/
  }

  /* Test for some cases with no root between 0 and 1 */
  /* -------------------------------------------------------- */
  SameSign = (a+b+c+d)*d > 0;
  if (isMonotonic) {
    if (SameSign) return 0;
  } else {
    yex = d + c*(2*a+b)/(3*a+2*b);
    if (SameSign && yex*d>0) return 0;
  }

  /* Test if root is at 0 */
  /* -------------------------------------------------------- */
  if (d==0) {
    nSolutions = 1;
    Solution[0] = 0;
    return nSolutions;
  }
  
  /* Try Newton-Raphson for a monotonic spline */
  /* --------------------------------------------------------- */
  if (isMonotonic) {
    t = -d/(a+b+c); /* starting guess: assume straight line  */

    for (iStep=0; iStep<maxNewtonSteps; iStep++) {
      delta = (((a*t+b)*t+c)*t+d)/((3*a*t+2*b)*t+c);
      t -= delta;
      if (t<0 || t>1)
        break;
      if (fabs(delta)<eps) {
        nSolutions = 1;
        Solution[0] = t;
        break;
      }
    } /* for (Step...)    */
  } /* if (isMonotonic)   */

  /* Use full cubic solver, if Newton-Raphson failed or was not tried  */
  /*------------------------------------------------------------------*/
  if (nSolutions==0) {
    nSolutions = exactSolver(a,b,c,d,Solution);
  }

  /* Compute first and second derivative at root positions */
  /*------------------------------------------------------- */
  for (i=0;i<nSolutions;i++) {
    t = Solution[i];
    Solution[i+3] = (3*a*t+2*b)*t + c;
    Solution[i+6] = 6*a*t + 2*b;
  }

  return nSolutions;

}

/*
****************************************************************
           MEX gateway function for MATLAB
****************************************************************
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* input parameters */
  double *Coeffs;
  bool LoopFields;
  /* output parameters */
  double *Solution, *roots, *deri1, *deri2, t;
  int nSolutions, i, j;

  /* Argument checks. */
  /*
  if (nlhs>3)
    mexErrMsgTxt("Wrong number of output parameters!");
  if (nrhs<1 || nrhs>2)
    mexErrMsgTxt("Wrong number of input parameters!");
  if (mxGetNumberOfElements(prhs[0])!=4 || mxIsComplex(prhs[0]))
    mexErrMsgTxt("Coefficient vector must contain 4 real elements!");
  */
  /* Get all input arguments.   */
  Coeffs = mxGetPr(prhs[0]);
  if (nrhs==2) {
    if (mxGetNumberOfElements(prhs[1])!=1)
      mexErrMsgTxt("LoopFields must be a scalar (boolean)!");
    LoopFields = mxGetScalar(prhs[1])!=0;
  }
  else {
    LoopFields = true;
  }

  /* Allocate result array. */
  Solution = (double *)mxMalloc(sizeof(double)*9);

  /* Call solver  */
  nSolutions = cubicsolve(Coeffs[0],Coeffs[1],Coeffs[2],Coeffs[3],!LoopFields,Solution);

  /* Sort solutions (bubble sorting)  */
  for (i=0;i<nSolutions;i++) {
    t = Solution[i];
    j = i;
    while (j>0 && Solution[j-1]>t) {
      Solution[j] = Solution[j-1];
      j--;
    }
    Solution[j] = t;
  }

  /* Allocate and assign output   */
  plhs[0] = mxCreateDoubleMatrix(1,nSolutions,mxREAL);
  
  roots = mxGetPr(plhs[0]);
  for (i=0;i<nSolutions;i++) {
    roots[i] = Solution[i];
  }
  
  if (nlhs>1) {
    plhs[1] = mxCreateDoubleMatrix(1,nSolutions,mxREAL);
    deri1 = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1,nSolutions,mxREAL);
    deri2 = mxGetPr(plhs[2]);
    for (i=0;i<nSolutions;i++) {
      deri1[i] = Solution[i+3];
      deri2[i] = Solution[i+6];
    }
  }

  mxFree(Solution);

} /* void mexFunction  */
