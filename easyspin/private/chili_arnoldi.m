% chili_arnoldi   Computes spectrum via Arnoldi iteration
%
% [spec,converged,specchange,nSteps] = chili_arnoldi(A,b,z,Opt)
%
% Computes spec(z) = b'*(A+z*I)^-1*b by projecting A onto the Krylov space of
% A and b using the Arnoldi iteration. Unlike the complex-symmetric Lanczos
% algorithm (chili_lanczos), Arnoldi does not break down, but it needs to store
% the Krylov basis (N x nSteps) and its cost grows quadratically with nSteps.
% To limit memory, the basis is capped at about 2 GB; if the spectrum has not
% converged by then, converged is returned as false.
%
% Inputs:
%   A           square NxN matrix
%   b           Nx1 starting vector
%   z           variable vector over which to evaluate the spectral function
%   Opt         structure with fields
%    .Threshold termination threshold for the change in the spectrum between
%               successive convergence checks
%
% Outputs:
%   spec        calculated spectrum
%   converged   whether spectrum is converged or not
%   specchange  array of spectral changes between successive convergence checks
%   nSteps      number of Arnoldi steps performed

function [spec,converged,specchange,nSteps] = chili_arnoldi(A,b,z,Opt)

N = length(b);
blockSize = 100;  % number of basis vectors to allocate at a time

% Maximum number of steps, limited by the memory for the Krylov basis
maxBytes = 2^31;
maxSteps = min(N,max(100,floor(maxBytes/(16*N))));

% Convergence checks get sparser as k grows (overshoot at most ~10% of steps),
% since each check costs O(k^3). Checks are at least 2 steps apart, since the
% projected spectrum can stagnate for one step (e.g. after a Lanczos breakdown).
checkInterval = max(2,min(10,ceil(N/20)));
nextCheck = checkInterval;
specchange = NaN(1,N);
oldspec = Inf;

% In the convergence tests, use horizontal axis that is shorter than the
% typical spectral length for performance reasons
if numel(z)>201
  zTest = linspace(z(1),z(end),201);
else
  zTest = z;
end

bnorm2 = real(b'*b);
V = zeros(N,min(blockSize,N));
V(:,1) = b/sqrt(bnorm2);
H = zeros(1,0);

converged = false;
for k = 1:maxSteps

  % Arnoldi step, using classical Gram-Schmidt with one reorthogonalization pass
  % (V is used in full, since V(:,1:k) would copy it; its unused columns are zero)
  w = A*V(:,k);
  h = V'*w;
  w = w - V*h;
  dh = V'*w;
  w = w - V*dh;
  H(1:k,k) = h(1:k) + dh(1:k);
  H(k+1,k) = norm(w);

  % Krylov space is invariant under A or complete: projection is exact
  exact = H(k+1,k)<=eps*norm(H(1:k,k)) || k==N;

  lastStep = exact || k==maxSteps;
  if k>=nextCheck || lastStep
    nextCheck = k + max(checkInterval,ceil(k/10));
    respec = real(arnoldispec(H(1:k,1:k),zTest));
    specchange(k) = max(abs(respec-oldspec))/max(abs(respec));
    converged = specchange(k)<Opt.Threshold || exact;
    oldspec = respec;
  end

  if converged || lastStep, break; end

  if k+1>size(V,2)
    V(:,min(maxSteps,k+blockSize)) = 0;
  end
  V(:,k+1) = w/H(k+1,k);

end

% Organize output
nSteps = k;
spec = bnorm2*arnoldispec(H(1:k,1:k),z);
specchange = specchange(1:k);
specchange(isnan(specchange)) = [];

end

%-------------------------------------------------------------------------------
% Evaluates e1'*(H+z*I)^-1*e1 for all z. The complex Schur decomposition
% H = U*T*U' (unitary U, upper triangular T) is numerically stable and reduces
% the evaluation for each z to a triangular solve.
function spec = arnoldispec(H,z)

[U,T] = schur(H,'complex');
c = U(1,:)';  % U'*e1
I = eye(size(T));
spec = zeros(size(z));
for j = 1:numel(z)
  spec(j) = c'*((T+z(j)*I)\c);
end

end
