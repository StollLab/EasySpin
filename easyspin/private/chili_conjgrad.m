% conjgrad  Conjugate gradients algorithm
%
%   [x,T0,T1] = conjgrad(A,b,shift);
%
%   Solves A*x = b using the conjugate gradient
%   algorithm. Computes the Lanczos tridiagonal
%   matrix of A using the starting vector b. The
%   diagonal is returned in T0, the upper
%   diagonal in T1.

function [x,T0,T1,err,StepsDone] = conjgrad(A,b,shift);

d = b;
x = 0;

r = b;
rr = r.'*r;
beta = 0;

iStep = 1;

MaxSteps = length(b);
Tolerance = 1e-8;
while 1
  
  if (iStep>1)
    rr_old = rr;
    rr = r.'*r;
    beta = rr/rr_old;
    d = r + beta*d;
  end
  bb(iStep) = beta;

  q = A*d + shift*d;
  alpha = rr/(d.'*q);
  aa(iStep) = alpha;

  err = sum(abs(r));
  if (err<Tolerance), break; end

  x = x + alpha*d;
  r = r - alpha*q;
  
  if (iStep==MaxSteps), break; end
  
  iStep = iStep + 1;
end
StepsDone = iStep;

%-------------------------------------------------------
% Compute tridiagonal Lanczos matrix from CG scalars
%-------------------------------------------------------
% Expressions:
%
%   T(k,k)     =   1/alpha(k) + beta(k)/alpha(k-1)
%   T(k-1,k)   =   -sqrt(beta(k))/alpha(k-1)
%   T(k,k-1)   =   T(k-1,k)
%
% T is the Lanczos tridiagonal matrix, alpha and beta
% are quantities from the conjugate gradient algorithm.

if (nargout>1)
  for k = 2:StepsDone
    t1 = sqrt(bb(k))/aa(k-1);
    T1(k-1) = t1*sign(real(t1));
    %T1(k-1) = t1;
    T0(k)  = 1/aa(k) + bb(k)/aa(k-1);
  end
  T0(1) = 1/aa(1);
  T0 = T0 - shift;
  T1(StepsDone) = 0;
end
