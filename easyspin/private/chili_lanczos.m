% chili_lanczos   Computes spectrum via Lanczos tridiagonalization
%
% [spec,converged,specchange] = chili_lanczos(A,b,z,Opt)
%
% Inputs:
%   A           complex symmetrix square NxN matrix
%   b           Nx1 starting vector
%   z           variable vector over which to evaluate the spectral function
%   Opt         structure with fields
%    .Lentz     if true (default), use Lentz method for the evaluation of the
%               continued fraction expansion
%    .Threshold termination threshold for iteration-to-iteration change in
%               spectrum
%
% Outputs:
%   spec        calculated spectrum
%   converged   whether spectrum is converged or not
%   specchange  array if iteration-to-iteration spectral changes

% The modified Lentz method is implemented from
%    W. H. Press et al, Numerical Recipes in C, 2nd edition
%    section 5.2, p.171

function [spec,converged,specchange] = chili_lanczos(A,b,z,Opt)

N = length(b);
alpha = zeros(1,N);
beta = zeros(1,N);

specchange = NaN(1,N);
oldspec = inf;

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

% Lanczos iterations
q = b/sqrt(b.'*b); % important: pseudonorm/rectanorm b.'*b instead of b'*b
bq = 0;
for k = 1:N
  
  % Lanczos step
  y = A*q;
  alpha(k) = q.'*y; 
  y = y - alpha(k)*q - bq;
  beta(k) = sqrt(y.'*y);
  
  bq = beta(k)*q;
  q = y/beta(k);

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
    if ~mod(k,errorRecomputationInterval)
      specchange(k) = norm(Delta-1,inf);
      converged = specchange(k)<Opt.Threshold;
      %fprintf('%3d: %g\n',k,specchange(k));
    end
  else
    if ~mod(k,errorRecomputationInterval)
      spec = chili_contfracspec(zTest,alpha,beta,k);
      respec = real(spec);
      specchange(k) = max(abs(respec-oldspec))/max(respec);
      converged = specchange(k)<Opt.Threshold;
      oldspec = respec;
    end
  end
  
  if converged, break; end
  
  if beta(k)==0
    specchange(k) = 0;
    break
  end
  
end

% Organize output
specchange = specchange(1:k);
specchange(isnan(specchange)) = [];

return
