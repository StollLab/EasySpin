% cardamom_buildhmm    Hidden Markov model (HMM) for spin label dynamics
%
%   HMM = cardamom_buildhmm(dihedrals,nStates,nLag,Opt)
%
% Input:
%   dihedrals     3D array of dihedrals, (nDims,nSteps,nTrajectories), in radians
%   nStates       number of states for the HMM
%   nLag          lag time, as a integer multiple of the MD time step
%   Opt           structure with options
%     .Verbosity  print to command window if > 0
%     .isSeeded   whether to use systematic seeds for the centroids in k-means
%     .nTrials    number of trials in k-means clustering (if not seeded)
%
% Output:
%   HMM         structure with HMM parameters
%    .transmat  transition probability matrix
%    .eqdistr   equilibrium distribution vector
%    .stateTraj state trajectory

function HMM = cardamom_buildhmm(dihedrals,nStates,nLag,Opt)

logmsg(1,'  HMM model building ----------------------------------');
logmsg(1,'  data: %d dihedrals; %d steps; %d trajectories',...
  size(dihedrals,1),size(dihedrals,3),size(dihedrals,2));

if ~isfield(Opt,'isSeeded')
  Opt.isSeeded = false;
end

if ~isfield(Opt,'nTrials')
  Opt.nTrials = 10;
end
if Opt.isSeeded
  Opt.nTrials = 1;
end

if abs(nLag-round(nLag))>1e-5 || nLag < 1
  error('nLag must be an integer >= 1.');
end
nLag = round(nLag);

% Set up initial cluster centroids if wanted
%-------------------------------------------------------------------------------
chiStart = [];
if Opt.isSeeded
  logmsg(1,'  using provided seeds');

  nDims = size(dihedrals,1);
  if (nDims==4) || (nDims==5)
    % Empirically determined values for polyala
    %chi = {[-60,65,180],[75,180],[-90,90],[75,8,-100],[180,77]};
    % Empirically determined values for T4L V131R1
    %chi = {[-60,180],[-55,55,180],[-90,90],[-170,-100,60],[-100,-20,100]};
    % Theoretical values
    chi = {[-60,60,180],[-60,60,180],[-90,90],[-60,60,180],[-90,90]};
    
    % Remove chi3 if it is absent from the dihedrals trajectories
    if size(dihedrals,1)==4
      chi(3) = [];
    end
    
    % Convert from degrees to radians
    chi = cellfun(@(x)x*pi/180,chi, 'UniformOutput', false);
    
    % Create array with all combinations
    idx = cellfun(@(a)1:numel(a),chi,'UniformOutput',false);
    [idx{:}] = ndgrid(idx{:});
    chiStart = [];
    for k = numel(chi):-1:1
      chiStart(:,k) = reshape(chi{k}(idx{k}),[],1);
    end

  end
else
  logmsg(1,'  using random seeds');
end

% Use k-means clustering etc to get initial estimates of HMM parameters
%-------------------------------------------------------------------------------

% Reorder from (nDims,nTraj,nSteps) to (nSteps,nDims,nTraj),
% for input to clustering function.
dihedrals = permute(dihedrals,[3,1,2]);

% Perform k-means clustering, return centroids mu0 and spreads Sigma0
logmsg(1,'  clustering into %d clusters using k-means (%d repeats)',nStates,Opt.nTrials);
[stateTraj, mu0, Sigma0] = ...
  initializehmm(dihedrals, chiStart, nStates, Opt.nTrials, Opt.Verbosity);

% Print results of clustering
if Opt.Verbosity >= 1
  fprintf('    cluster population  max(stddev)/deg  mu0/deg\n');
  for k = 1:nStates
    pop(k) = sum(stateTraj==k)/numel(stateTraj);
    stddev(k) = sqrt(max(eig(Sigma0(:,:,k))));
  end
  for k = 1:nStates
    fprintf('     %3d      %0.4f       %6.2f',k,pop(k),stddev(k)*180/pi);
    fprintf('        (')
    for d = 1:size(mu0,1)
      fprintf('%4.0f ',mu0(d,k)*180/pi);
    end
    fprintf(')\n');
  end
end

logmsg(1,'  MSM parameter estimation');

% Downsample dihedrals trajectory to the desired lag time
dihedrals = dihedrals(1:nLag:end,:,:);
stateTraj = stateTraj(1:nLag:end,:);

% Estimate transition probability matrix and initial distribution
[transmat0,eqdistr0] = estimatemarkovparameters(stateTraj);

% Reorder (nSteps,nDims,nTraj) to (nDims,nSteps,nTraj), for EM function
dihedrals = permute(dihedrals,[2,1,3]);

% Optimize HMM parameters
%-------------------------------------------------------------------------------
% Determine/estimate HMM model parameters using expectation maximization
[HMM.eqdistr,HMM.transmat,HMM.mu,HMM.Sigma] = ...
  cardamom_emghmm(dihedrals,eqdistr0,transmat0,mu0,Sigma0,[],Opt.Verbosity);

% Calculate Viterbi state trajectory
%-------------------------------------------------------------------------------
logmsg(1,'  Viterbi state trajectories calculation');
% Determine most probable hidden-state trajectory
HMM.viterbiTraj = viterbitrajectory(dihedrals,HMM.transmat,HMM.eqdistr,HMM.mu,HMM.Sigma);

HMM.nLag = nLag;

end

% Helper functions 
%-------------------------------------------------------------------------------
function vTraj = viterbitrajectory(dihedrals,transmat,eqdistr,mu,Sigma)
nStates = size(transmat,1);
[nDims,nSteps,nTraj] = size(dihedrals);
vTraj = zeros(nSteps,nTraj);
for iTraj = 1:nTraj
%   [obslikelihood, ~] = mixgauss_prob(dihedrals(:,:,iTraj), mu, Sigma);
  obs = dihedrals(:,:,iTraj);
  obslikelihood = zeros(nStates,nSteps);
  for iState=1:nStates
    obslikelihood(iState,:) = gaussian_prob_pbc(obs, mu(:,iState), Sigma(:,:,iState), 2*pi);
  end
  vTraj(:,iTraj) = viterbi_path(eqdistr, transmat, obslikelihood).';
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = cov_pbc(x, mu, W)

[nPoints,nDims] = size(x);
    
% remove the centers
xc = zeros(nPoints,nDims);
for iPoint = 1:nPoints
  xc(iPoint,:) = dist_pbc(x(iPoint,:)-mu, W);
end

c = (xc' * xc) ./ nPoints;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist = dist_pbc(dist,W)

w = W/2;

% dist = w - abs(pi - abs(dist));

idx1 = dist > w;
idx2 = dist < -w;

dist(idx1) = dist(idx1) - W;
dist(idx2) = dist(idx2) + W;

% for iDim=1:size(dist,2)
%   if dist(iDim) > w
%     dist(iDim) = dist(iDim)-W;
%   elseif dist(iDim) < -w
%     dist(iDim) = dist(iDim)+W;
%   end
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stateTraj,mu0,Sigma0] = initializehmm(dihedrals,chiStart,nStates,nRepeats,verbosity)

[nSteps,nDims,nTraj] = size(dihedrals);

% if more than one trajectory, collapse 3rd dim (traj) onto 1st dim (time)
if nTraj > 1
  dihedralsTemp = dihedrals;
  dihedrals = [];
  for iTraj = 1:nTraj
    dihedrals = cat(1, dihedrals, dihedralsTemp(:,:,iTraj));
  end
end

% Do k-means clustering
[stateTraj,centroids] = ...
  cardamom_kmeans(dihedrals, nStates, nRepeats, chiStart, verbosity);

% initialize the means and covariance matrices for the HMM
mu0 = centroids.';
Sigma0 = zeros(nDims,nDims,nStates);
for iState = 1:nStates
  idxState = stateTraj==iState;
  Sigma0(:,:,iState) = cov_pbc(dihedrals(idxState,:), mu0(:,iState).', 2*pi);
%   Sigma0(:,:,iState) = cov(dihedrals(idxState,:));
end

stateTraj = reshape(stateTraj,[nSteps,nTraj]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tpmat,distr] = estimatemarkovparameters(stateTraj)
% Estimate transition probability matrix and initial probability distribution
% from a set of state trajectories.
% Input:
%    stateTraj  array (nSteps,nTraj) of state indices (1..nStates)
% Output:
%    tpmat      transition probability matrix (nStates,nStates)
%    distr      initial probability density (nStates,1)

% Determine count matrix N
%-------------------------------------------------------------------------------
% N(i,j) = number of times X(t)==i and X(t+tau)==j along all the trajectories
nStates = max(max(stateTraj));
N = sparse(nStates,nStates);

% Time step for which to determine the count matrix
tau = 1;

% Calculate (and accumulate) count matrix N
nTraj = size(stateTraj,2);
for iTraj = 1:nTraj
  state1 = stateTraj(1:end-tau,iTraj);
  state2 = stateTraj(1+tau:end,iTraj);
  N = N + sparse(state1,state2,1,nStates,nStates);
end
N = full(N);

% Calculate transition probability matrix and initial distribution
%-------------------------------------------------------------------------------
[tpmat,distr] = msmtransitionmatrix(N);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = gaussian_prob_pbc(x, m, C, W)
% GAUSSIAN_PROB Evaluate a multivariate Gaussian density.
% p = gaussian_prob(X, m, C)
% p(i) = N(X(:,i), m, C) where C = covariance matrix and each COLUMN of x is a datavector

% p = gaussian_prob(X, m, C, 1) returns log N(X(:,i), m, C) (to prevents underflow).
%
% If X has size dxN, then p has size Nx1, where N = number of examples

if length(m)==1 % scalar
  x = x(:)';
end
[d, N] = size(x);
%assert(length(m)==d); % slow
m = m(:);
M = m*ones(1,N); % replicate the mean across columns
denom = (2*pi)^(d/2)*sqrt(abs(det(C)));
dist = dist_pbc(x - M, W);
mahal = sum((dist'*inv(C)).*dist',2);   % Chris Bregler's trick
% mahal = sum(((x-M)'*inv(C)).*(x-M)',2);   % Chris Bregler's trick
if any(mahal<0)
  warning('mahal < 0 => C is not psd')
end
p = exp(-0.5*mahal) / (denom+eps);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B, B2] = mixgauss_prob(data, mu, Sigma, mixmat, unit_norm)
% EVAL_PDF_COND_MOG Evaluate the pdf of a conditional mixture of Gaussians
% function [B, B2] = eval_pdf_cond_mog(data, mu, Sigma, mixmat, unit_norm)
%
% Notation: Y is observation, M is mixture component, and both may be conditioned on Q.
% If Q does not exist, ignore references to Q=j below.
% Alternatively, you may ignore M if this is a conditional Gaussian.
%
% INPUTS:
% data(:,t) = t'th observation vector 
%
% mu(:,k) = E[Y(t) | M(t)=k] 
% or mu(:,j,k) = E[Y(t) | Q(t)=j, M(t)=k]
%
% Sigma(:,:,j,k) = Cov[Y(t) | Q(t)=j, M(t)=k]
% or there are various faster, special cases:
%   Sigma() - scalar, spherical covariance independent of M,Q.
%   Sigma(:,:) diag or full, tied params independent of M,Q. 
%   Sigma(:,:,j) tied params independent of M. 
%
% mixmat(k) = Pr(M(t)=k) = prior
% or mixmat(j,k) = Pr(M(t)=k | Q(t)=j) 
% Not needed if M is not defined.
%
% unit_norm - optional; if 1, means data(:,i) AND mu(:,i) each have unit norm (slightly faster)
%
% OUTPUT:
% B(t) = Pr(y(t)) 
% or
% B(i,t) = Pr(y(t) | Q(t)=i) 
% B2(i,k,t) = Pr(y(t) | Q(t)=i, M(t)=k) 
%
% If the number of mixture components differs depending on Q, just set the trailing
% entries of mixmat to 0, e.g., 2 components if Q=1, 3 components if Q=2,
% then set mixmat(1,3)=0. In this case, B2(1,3,:)=1.0.




if isvector(mu) & size(mu,2)==1
  d = length(mu);
  Q = 1; M = 1;
elseif ndims(mu)==2
  [d Q] = size(mu);
  M = 1;
else
  [d Q M] = size(mu);
end
[d T] = size(data);

if nargin < 4, mixmat = ones(Q,1); end
if nargin < 5, unit_norm = 0; end

%B2 = zeros(Q,M,T); % ATB: not needed allways
%B = zeros(Q,T);

if isscalar(Sigma)
  mu = reshape(mu, [d Q*M]);
  if unit_norm % (p-q)'(p-q) = p'p + q'q - 2p'q = n+m -2p'q since p(:,i)'p(:,i)=1
    %avoid an expensive repmat
    disp('unit norm')
    %tic; D = 2 -2*(data'*mu)'; toc 
    D = 2 - 2*(mu'*data);
    tic; D2 = sqdist(data, mu)'; toc
    assert(approxeq(D,D2)) 
  else
    D = sqdist(data, mu)';
  end
  clear mu data % ATB: clear big old data
  % D(qm,t) = sq dist between data(:,t) and mu(:,qm)
  logB2 = -(d/2)*log(2*pi*Sigma) - (1/(2*Sigma))*D; % det(sigma*I) = sigma^d
  B2 = reshape(exp(logB2), [Q M T]);
  clear logB2 % ATB: clear big old data
  
elseif ndims(Sigma)==2 % tied full
  mu = reshape(mu, [d Q*M]);
  D = sqdist(data, mu, inv(Sigma))';
  % D(qm,t) = sq dist between data(:,t) and mu(:,qm)
  logB2 = -(d/2)*log(2*pi) - 0.5*logdet(Sigma) - 0.5*D;
  %denom = sqrt(det(2*pi*Sigma));
  %numer = exp(-0.5 * D);
  %B2 = numer/denom;
  B2 = reshape(exp(logB2), [Q M T]);
  
elseif ndims(Sigma)==3 % tied across M
  B2 = zeros(Q,M,T);
  for j=1:Q
    % D(m,t) = sq dist between data(:,t) and mu(:,j,m)
    if isposdef(Sigma(:,:,j))
      D = sqdist(data, permute(mu(:,j,:), [1 3 2]), inv(Sigma(:,:,j)))';
      logB2 = -(d/2)*log(2*pi) - 0.5*logdet(Sigma(:,:,j)) - 0.5*D;
%       logB2 = -(d/2)*log(2*pi) - 0.5*log(det(Sigma(:,:,j))) - 0.5*D;
      B2(j,:,:) = exp(logB2);
    else
      error(sprintf('mixgauss_prob: Sigma(:,:,q=%d) not psd\n', j));
    end
  end
  
else % general case
  B2 = zeros(Q,M,T);
  for j=1:Q
    for k=1:M
      %if mixmat(j,k) > 0
      B2(j,k,:) = gaussian_prob(data, mu(:,j,k), Sigma(:,:,j,k));
      %end
    end
  end
end

% B(j,t) = sum_k B2(j,k,t) * Pr(M(t)=k | Q(t)=j) 

% The repmat is actually slower than the for-loop, because it uses too much memory
% (this is true even for small T).

%B = squeeze(sum(B2 .* repmat(mixmat, [1 1 T]), 2));
%B = reshape(B, [Q T]); % undo effect of squeeze in case Q = 1
  
B = zeros(Q,T);
if Q < T
  for q=1:Q
    %B(q,:) = mixmat(q,:) * squeeze(B2(q,:,:)); % squeeze chnages order if M=1
    B(q,:) = mixmat(q,:) * permute(B2(q,:,:), [2 3 1]); % vector * matrix sums over m
  end
else
  for t=1:T
    B(:,t) = sum(mixmat .* B2(:,:,t), 2); % sum over m
  end
end
%t=toc;fprintf('%5.3f\n', t)

%tic
%A = squeeze(sum(B2 .* repmat(mixmat, [1 1 T]), 2));
%t=toc;fprintf('%5.3f\n', t)
%assert(approxeq(A,B)) % may be false because of round off error

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = isposdef(a)
% ISPOSDEF   Test for positive definite matrix.
%    ISPOSDEF(A) returns 1 if A is positive definite, 0 otherwise.
%    Using chol is much more efficient than computing eigenvectors.

%  From Tom Minka's lightspeed toolbox

[R,p] = chol(a);
b = (p == 0);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = logdet(A)
% log(det(A)) where A is positive-definite.
% This is faster and more stable than using log(det(A)).

%  From Tom Minka's lightspeed toolbox

U = chol(A);
y = 2*sum(log(diag(U)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = sqdist(p, q, A)
% SQDIST      Squared Euclidean or Mahalanobis distance.
% SQDIST(p,q)   returns m(i,j) = (p(:,i) - q(:,j))'*(p(:,i) - q(:,j)).
% SQDIST(p,q,A) returns m(i,j) = (p(:,i) - q(:,j))'*A*(p(:,i) - q(:,j)).

%  From Tom Minka's lightspeed toolbox

[d, pn] = size(p);
[d, qn] = size(q);

if nargin == 2
  
  pmag = sum(p .* p, 1);
  qmag = sum(q .* q, 1);
  m = repmat(qmag, pn, 1) + repmat(pmag', 1, qn) - 2*p'*q;
  %m = ones(pn,1)*qmag + pmag'*ones(1,qn) - 2*p'*q;
  
else

  if isempty(A) | isempty(p)
    error('sqdist: empty matrices');
  end
  Ap = A*p;
  Aq = A*q;
  pmag = sum(p .* Ap, 1);
  qmag = sum(q .* Aq, 1);
  m = repmat(qmag, pn, 1) + repmat(pmag', 1, qn) - 2*p'*Aq;
  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function path = viterbi_path(prior, transmat, obslik)
% VITERBI Find the most-probable (Viterbi) path through the HMM state trellis.
% path = viterbi(prior, transmat, obslik)
%
% Inputs:
% prior(i) = Pr(Q(1) = i)
% transmat(i,j) = Pr(Q(t+1)=j | Q(t)=i)
% obslik(i,t) = Pr(y(t) | Q(t)=i)
%
% Outputs:
% path(t) = q(t), where q1 ... qT is the argmax of the above expression.


% delta(j,t) = prob. of the best sequence of length t-1 and then going to state j, and O(1:t)
% psi(j,t) = the best predecessor state, given that we ended up in state j at t

scaled = 1;

T = size(obslik, 2);
prior = prior(:);
Q = length(prior);

delta = zeros(Q,T);
psi = zeros(Q,T);
path = zeros(1,T);
scale = ones(1,T);


t=1;
delta(:,t) = prior .* obslik(:,t);
if scaled
  [delta(:,t), n] = normalise(delta(:,t));
  scale(t) = 1/n;
end
psi(:,t) = 0; % arbitrary value, since there is no predecessor to t=1
for t=2:T
  for j=1:Q
    [delta(j,t), psi(j,t)] = max(delta(:,t-1) .* transmat(:,j));
    delta(j,t) = delta(j,t) * obslik(j,t);
  end
  if scaled
    [delta(:,t), n] = normalise(delta(:,t));
    scale(t) = 1/n;
  end
end
[p, path(T)] = max(delta(:,T));
for t=T-1:-1:1
  path(t) = psi(path(t+1),t+1);
end

% If scaled==0, p = prob_path(best_path)
% If scaled==1, p = Pr(replace sum with max and proceed as in the scaled forwards algo)
% Both are different from p(data) as computed using the sum-product (forwards) algorithm

if 0
  if scaled
    loglik = -sum(log(scale));
    %loglik = prob_path(prior, transmat, obslik, path);
  else
    loglik = log(p);
  end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M, z] = normalise(A, dim)
% NORMALISE Make the entries of a (multidimensional) array sum to 1
% [M, c] = normalise(A)
% c is the normalizing constant
%
% [M, c] = normalise(A, dim)
% If dim is specified, we normalise the specified dimension only,
% otherwise we normalise the whole array.

if nargin < 2
  z = sum(A(:));
  % Set any zeros to one before dividing
  % This is valid, since c=0 => all i. A(i)=0 => the answer should be 0/1=0
  s = z + (z==0);
  M = A / s;
elseif dim==1 % normalize each column
  z = sum(A);
  s = z + (z==0);
  %M = A ./ (d'*ones(1,size(A,1)))';
  M = A ./ repmatC(s, size(A,1), 1);
else
  % Keith Battocchi - v. slow because of repmat
  z=sum(A,dim);
  s = z + (z==0);
  L=size(A,dim);
  d=length(size(A));
  v=ones(d,1);
  v(dim)=L;
  %c=repmat(s,v);
  c=repmat(s,v');
  M=A./c;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,Z] = mk_stochastic(T)
% MK_STOCHASTIC Ensure the argument is a stochastic matrix, i.e., the sum over the last dimension is 1.
% [T,Z] = mk_stochastic(T)
%
% If T is a vector, it will sum to 1.
% If T is a matrix, each row will sum to 1.
% If T is a 3D array, then sum_k T(i,j,k) = 1 for all i,j.

% Set zeros to 1 before dividing
% This is valid since S(j) = 0 iff T(i,j) = 0 for all j

if (ndims(T)==2) & (size(T,1)==1 | size(T,2)==1) % isvector
  [T,Z] = normalise(T);
elseif ndims(T)==2 % matrix
  Z = sum(T,2); 
  S = Z + (Z==0);
  norm = repmat(S, 1, size(T,2));
  T = T ./ norm;
else % multi-dimensional array
  ns = size(T);
  T = reshape(T, prod(ns(1:end-1)), ns(end));
  Z = sum(T,2);
  S = Z + (Z==0);
  norm = repmat(S, 1, ns(end));
  T = T ./ norm;
  T = reshape(T, ns);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TPM, pi_i, x] = msmtransitionmatrix(N, maxiteration)
% msmtransitionmatrix
% estimate transition probability matrix TPM from count matrix N
%
% Syntax
% [TPM, pi_i] = msmtransitionmatrix(N);
% [TPM, pi_i] = msmtransitionmatrix(N, maxiteration);
%
% Description
% this routines uses the reversible maximum likelihood estimator
%
% Adapted from Yasuhiro Matsunaga's mdtoolbox
% 

if ~exist('maxiteration', 'var')
  maxiteration = 1000;
end

% setup
nStates = size(N, 1);

N_sym = N + N.'; % symmetrize count matrix
x = N_sym;

N_i = sum(N, 2);
x_i = sum(x, 2);

% optimization by L-BFGS-B
fcn = @(x) myfunc_column(x, N, N_i, nStates);
opts.x0 = x(:);
opts.maxIts = maxiteration;
opts.maxTotalIts = 50000;
%opts.factr = 1e5;
%opts.pgtol = 1e-7;

[x, f, info] = cardamom_lbfgsb(fcn, zeros(nStates*nStates, 1), Inf(nStates*nStates, 1), opts);
x = reshape(x,nStates,nStates);

x_i = sum(x,2);
TPM = bsxfun(@rdivide,x,x_i);
TPM(isnan(TPM)) = 0;
pi_i = x_i./sum(x_i);


%-------------------------------------------------------------------------------
function [f, g] = myfunc_column(x, c, c_i, nStates)
x = reshape(x,nStates,nStates);
[f, g] = myfunc_matrix(x, c, c_i);
g = g(:);

end

%-------------------------------------------------------------------------------
function [f, g] = myfunc_matrix(x, c, c_i)
x_i = sum(x, 2);

% F
tmp = c .* log(bsxfun(@rdivide, x, x_i));
%index = ~(isnan(tmp));
index = (x > (10*eps));
f = - sum(tmp(index));

% G
t = c_i./x_i;
g = (c./x) + (c'./x') - bsxfun(@plus, t, t');
g((x_i < (10*eps)), :) = 0;
index = ((x > (10*eps)) & (x' > (10*eps)));
g(~index) = 0;
g(isnan(g)) = 0;
g = -g;

end

end