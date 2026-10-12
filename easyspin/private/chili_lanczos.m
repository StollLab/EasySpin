% chili_lanczos   Computes spectrum via Lanczos tridiagonalization
%
% [spec,converged,specchange,nSteps,breakdown] = chili_lanczos(A,b,z,Opt)
%
% Inputs:
%   A           complex symmetric square NxN matrix
%   b           Nx1 starting vector
%   z           variable vector over which to evaluate the spectral function
%   Opt         structure with fields
%    .Lentz     if true (default), use Lentz method for the evaluation of the
%               continued fraction expansion
%    .Threshold termination threshold for iteration-to-iteration change in
%               spectrum
%    .LanczosBreakdownTol  threshold for detecting a (near) serious breakdown,
%               i.e. a vector y with |y.'*y| << y'*y
%
% Outputs:
%   spec        calculated spectrum
%   converged   whether spectrum is converged or not
%   specchange  array of iteration-to-iteration spectral changes
%   nSteps      number of Lanczos steps performed
%   breakdown   true if the iteration stopped due to a (near) serious breakdown

% The modified Lentz method is implemented from
%    W. H. Press et al, Numerical Recipes in C, 2nd edition
%    section 5.2, p.171

function [spec,converged,specchange,nSteps,breakdown] = chili_lanczos(A,b,z,Opt)

N = length(b);
alpha = zeros(1,N);
beta = zeros(1,N);

specchange = NaN(1,N);
oldspec = Inf;

% Initialize continued-fraction evaluations
errorRecomputationInterval = min(10,ceil(N/20));
useLentzMethod = Opt.Lentz;
if useLentzMethod
  tiny = 1e-30;
  spec = tiny;
  C = spec;
  D = 0;
else
  % In the convergence tests, use horizontal axis that is shorter than the
  % typical spectral length for performance reasons
  zTest = linspace(z(1),z(end),201);
end

converged = false;
breakdown = false;

% Lanczos iterations
q = b/sqrt(b.'*b); % important: pseudonorm/rectanorm b.'*b instead of b'*b
bq = 0;
for k = 1:N
  
  % Lanczos step
  y = A*q;
  normAq = norm(y);
  alpha(k) = q.'*y;
  y = y - alpha(k)*q - bq;
  pseudonorm2 = y.'*y;
  norm2 = real(y'*y);

  % Detect termination: y is negligible, so the Krylov space is invariant and
  % the continued fraction up to step k is exact (the tolerance allows for
  % rounding errors amplified by the non-unitary Lanczos vectors)
  terminated = sqrt(norm2)<=1e-8*normAq;

  % Detect (near) serious breakdown: y is not small, but its pseudonorm is.
  % The tridiagonalization is then unstable and the spectrum unreliable.
  if ~terminated && k<N && abs(pseudonorm2)<Opt.LanczosBreakdownTol*norm2
    breakdown = true;
    break
  end

  % Continued fraction: next convergent
  if useLentzMethod
    Cold = C;
    Dold = D;
    specold = spec;
    b = alpha(k) + z; 
    if k==1
      a = 1;
    else
      a = -beta(k-1)^2;
    end
    D = b + a.*Dold;
    %D(D==0) = tiny;
    C = b + a./Cold;
    %C(C==0) = tiny;
    D = 1./D;
    Delta = C.*D;
    spec = specold.*Delta;
    if ~mod(k,errorRecomputationInterval) || k==N
      specchange(k) = norm(Delta-1,Inf);
      converged = specchange(k)<Opt.Threshold;
      %fprintf('%3d: %g\n',k,specchange(k));
    end
  else
    if ~mod(k,errorRecomputationInterval) || k==N
      spec = chili_contfracspec(zTest,alpha,beta,k);
      respec = real(spec);
      specchange(k) = max(abs(respec-oldspec))/max(respec);
      converged = specchange(k)<Opt.Threshold;
      oldspec = respec;
    end
  end

  if terminated
    converged = true;
  end
  if converged, break; end

  beta(k) = sqrt(pseudonorm2);
  bq = beta(k)*q;
  q = y/beta(k);

end

% Organize output
if ~useLentzMethod
  spec = chili_contfracspec(z,alpha,beta,k);
end
nSteps = k;
specchange = specchange(1:k);
specchange(isnan(specchange)) = [];

end
