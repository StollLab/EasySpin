% conjgrad  Conjugate gradients algorithm
%
%   [x,T0,T1] = conjgrad(A,b,shift)
%   
%   Solves A*x = b using the conjugate gradient algorithm. Computes the Lanczos
%   tridiagonal matrix of A using the starting vector b and shifting the diagonal
%   of A by shift.
%
%   The diagonal of the Lanczos tridiagonal matrix is returned in T0, and the
%   upper diagonal in T1.

function [x,T0,T1,err,stepsDone] = chili_conjgrad(A,b,shift)

if nargin<3
  shift = 0;
end

Tolerance = 1e-8;

% Initial guess: all zeros
x = zeros(size(b));


% Run Conjugate Gradient algorithm and save CG scalars alpha and beta
%-------------------------------------------------------------------------------
r = b - A*x; % residual vector
p = r;
rr = r.'*r; % square norm without complex conjugate

for iStep = 1:length(b)
  
  rr_old = rr;
  rr = r.'*r;
  if iStep==1
    beta_ = 0;
  else
    beta_ = rr/rr_old;
    p = r + beta_*p;
  end
  Ap = A*p + shift*p;
  alpha_ = rr/(p.'*Ap);
  x = x + alpha_*p;
  r = r - alpha_*Ap;
  
  % Store CG scalars
  alpha(iStep) = alpha_;
  beta(iStep) = beta_;
  
  % Terminate if residual is sufficiently small
  err = sum(abs(r));
  if err<Tolerance, break; end
  
end

stepsDone = iStep;
err = full(err);

% Compute tridiagonal Lanczos matrix T from CG scalars
%-------------------------------------------------------------------------------
% Expressions:
%
%   T(k,k) = 1/alpha(k) + beta(k)/alpha(k-1)
%   T(k-1,k)   =   -sqrt(beta(k))/alpha(k-1)
%   T(k,k-1)   =   T(k-1,k)
%
% T is the Lanczos tridiagonal matrix, alpha and beta
% are quantities from the conjugate gradient algorithm.

if nargout>1
  for k = 2:stepsDone
    t1 = sqrt(beta(k))/alpha(k-1);
    T1(k-1) = t1*sign(real(t1));
    %T1(k-1) = t1;
    T0(k)  = 1/alpha(k) + beta(k)/alpha(k-1);
  end
  T0(1) = 1/alpha(1);
  T0 = T0 - shift;
  T1(stepsDone) = 0;
end
