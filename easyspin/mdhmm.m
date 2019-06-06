% mdhmm    Build hidden Markov model (HMM) from MD trajectory of dihedrals
%
%   HMM = mdhmm(dihedrals,dt,nStates,nLag,Opt)
%
% Input:
%   dihedrals     3D array of spin-label side chain dihedral angles
%                   (nDims,nTrajectories,nSteps), in radians
%   dt            MD time step, in s
%   nStates       number of desired states for the HMM
%   nLag          desired lag time, as a integer multiple of the MD time step
%   Opt           structure with options
%     .Verbosity  print to command window if > 0
%     .isSeeded   whether to use systematic seeds for the centroids in k-means
%     .nTrials    number of trials in k-means clustering (if not seeded)
%
% Output:
%   HMM             structure with HMM parameters
%    .TransProb     transition probability matrix
%    .eqDistr       equilibrium distribution vector
%    .mu            center vectors of states
%    .Sigma         covariance matrices of states
%    .viterbiTraj   Viterbi state trajectory (most likely given the dihedrals)
%    .tauRelax      relaxation times of HMM

function HMM = mdhmm(dihedrals,dt,nStates,nLag,Opt)

if nargin<5, Opt = struct; end

if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 1; end

global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

logmsg(1,'-- HMM model building ----------------------------------');
nDims = size(dihedrals,1);
nTraj = size(dihedrals,2);
nSteps = size(dihedrals,3);
logmsg(1,'  data: %d dihedrals; %d steps; %d trajectories',nDims,nSteps,nTraj);

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
HMM.nLag = nLag;
HMM.dt = dt;
HMM.tLag = nLag*dt;
HMM.nStates = nStates;

% Use k-means clustering etc to get initial estimates of HMM parameters
%-------------------------------------------------------------------------------
logmsg(1,'  clustering into %d clusters using k-means (%d repeats)',nStates,Opt.nTrials);

% Set up initial cluster centroids if wanted
chiStart = [];
if Opt.isSeeded
  logmsg(1,'    using provided seeds');
  
  nDims = size(dihedrals,1);
  if nDims==5
    % Theoretical values
    chi = {[-60,60,180],[-60,60,180],[-90,90],[-60,60,180],[-90,90]};
    
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
  logmsg(1,'    using random seeds');
end

% Reorder from (nDims,nTraj,nSteps) to (nSteps,nDims,nTraj),
% for input to clustering function.
dihedrals = permute(dihedrals,[3,1,2]);

% Perform k-means clustering, return centroids mu0 and spreads Sigma0
[stateTraj,mu0,Sigma0] = ...
  initializehmm(dihedrals,chiStart,nStates,Opt.nTrials,Opt.Verbosity);

% Print results of clustering
if Opt.Verbosity >= 1
  fprintf('    cluster  population  max(stddev)/deg  mu0/deg\n');
  for k = 1:nStates
    pop(k) = sum(stateTraj(:)==k)/numel(stateTraj);
    stddev(k) = sqrt(max(eig(Sigma0(:,:,k))));
  end
  for k = 1:nStates
    fprintf('     %3d       %0.4f       %6.1f',k,pop(k),stddev(k)*180/pi);
    fprintf('        (')
    for d = 1:size(mu0,1)
      fprintf('%4.0f ',mu0(d,k)*180/pi);
    end
    fprintf(')\n');
  end
end

logmsg(1,'  estimation of transition probability matrix and initial distribution');

% Downsample dihedrals trajectory to the desired lag time
dihedrals = dihedrals(1:nLag:end,:,:);
stateTraj = stateTraj(1:nLag:end,:);

% Estimate transition probability matrix and initial distribution
[TransProb0,eqDistr0] = estimatemarkovparameters(stateTraj);
initDistr0 = eqDistr0;

% Reorder (nSteps,nDims,nTraj) to (nDims,nSteps,nTraj), for EM function
dihedrals = permute(dihedrals,[2,1,3]);

% Optimize HMM parameters
%-------------------------------------------------------------------------------
% Determine/estimate HMM model parameters using expectation maximization
logmsg(1,'  HMM optimization using EM algorithm');
[logLik, HMM.eqDistr,HMM.TransProb,HMM.mu,HMM.Sigma] = ...
  mdhmm_em(dihedrals,initDistr0,TransProb0,mu0,Sigma0,Opt.Verbosity);
HMM.logLik = logLik(end);

% Calculate Viterbi state trajectory
%-------------------------------------------------------------------------------
% Determine most probable hidden-state trajectory
logmsg(1,'  Viterbi state trajectories calculation');
HMM.viterbiTraj = viterbitrajectory(dihedrals,HMM.TransProb,HMM.eqDistr,HMM.mu,HMM.Sigma);

% Eliminate states not visited in Viterbi trajectory
%-------------------------------------------------------------------------------
visited = false(1,nStates);
visited(unique(HMM.viterbiTraj)) = true;
if any(~visited)
  logmsg(1,'  Eliminating %d unvisited states from model',sum(~visited));
  newStateNumbers = cumsum(visited);
  HMM.viterbiTraj = newStateNumbers(HMM.viterbiTraj).';
  HMM.mu = HMM.mu(:,visited);
  HMM.Sigma = HMM.Sigma(:,:,visited);
  [HMM.TransProb,HMM.eqDistr] = estimatemarkovparameters(HMM.viterbiTraj);
end
HMM.nStates = length(HMM.eqDistr);

% Calculate relaxation times for the TPM and time lag
%-------------------------------------------------------------------------------
logmsg(1,'  Calculate relaxation times');
lambda = eig(HMM.TransProb);
lambda = sort(real(lambda),1,'descend');
lambda = lambda(lambda>0);
HMM.tauRelax = -nLag*dt./log(lambda(2:end).'); % exclude constant eq. component
logmsg(1,'    max: %0.3g s, min: %0.3g s',max(HMM.tauRelax),min(HMM.tauRelax));

logmsg(1,'-- done ------------------------------------------------');

end
%===============================================================================
%===============================================================================


%===============================================================================
% Helper functions
%-------------------------------------------------------------------------------
function vTraj = viterbitrajectory(dihedrals,transmat,eqdistr,mu,Sigma)
nStates = size(transmat,1);
[~,nSteps,nTraj] = size(dihedrals);
vTraj = zeros(nSteps,nTraj);
for iTraj = 1:nTraj
  %[obslikelihood, ~] = mixgauss_prob(dihedrals(:,:,iTraj), mu, Sigma);
  obs = dihedrals(:,:,iTraj);
  obslikelihood = zeros(nStates,nSteps);
  for iState = 1:nStates
    obslikelihood(iState,:) = gaussian_prob_circular(obs,mu(:,iState),Sigma(:,:,iState));
  end
  vTraj(:,iTraj) = viterbi_path(eqdistr,transmat,obslikelihood).';
end

end

%===============================================================================
function c = cov_pbc(x, mu)
% Calculate covariance matrix, with samples x and mean mu

nPoints = size(x,1);
xc = dist_pbc(x,mu);
c = (xc' * xc) ./ nPoints;

end

%===============================================================================
function dist = dist_pbc(x1,x2)
% Distance between two angles, x1 and x2, in radians
% Takes the circular nature of angles into account

% One-line solution: (slower)
% dist = pi - abs(pi - abs(x1-x2));

dist = x1-x2;

idx1 = dist > pi;
idx2 = dist < -pi;

dist(idx1) = dist(idx1) - 2*pi;
dist(idx2) = dist(idx2) + 2*pi;

end

%===============================================================================
function [stateTraj,mu0,Sigma0] = initializehmm(dihedrals,chiStart,nStates,nRepeats,verbosity)

[nSteps,nDims,nTraj] = size(dihedrals);

% if more than one trajectory, collapse 3rd dim (traj) onto 1st dim (time)
if nTraj > 1
  dihedrals = reshape(dihedrals,[],nDims);
%   dihedralsTemp = dihedrals;
%   dihedrals = [];
%   for iTraj = 1:nTraj
%     dihedrals = cat(1, dihedrals, dihedralsTemp(:,:,iTraj));
%   end
end

% Do k-means clustering
[stateTraj,centroids] = ...
  mdhmm_kmeans(dihedrals, nStates, nRepeats, chiStart, verbosity);

% initialize the means and covariance matrices for the HMM
mu0 = centroids.';
Sigma0 = zeros(nDims,nDims,nStates);
for iState = 1:nStates
  idxState = stateTraj==iState;
  Sigma0(:,:,iState) = cov_pbc(dihedrals(idxState,:), mu0(:,iState).');
%   Sigma0(:,:,iState) = cov(dihedrals(idxState,:));
end

% Undo collapsing onto first dim for multiple trajectories
stateTraj = reshape(stateTraj,[nSteps,nTraj]);

end

%===============================================================================
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
nStates = max(stateTraj(:));
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

%===============================================================================
function p = gaussian_prob_circular(x, m, C)
% gaussian_prob_circular Evaluate a multivariate Gaussian density for angles
%
%   p = gaussian_prob_circular(x, m, C)
%
% Input:
%   x   data array, d x N
%         one d-dimensional datavector per column; N samples
%   m   mean, vector d x 1
%   C   covariance matrix, d x d
% Output:
%   p   probabilities, N x 1
%       p(i) = Normal(x(:,i), m, C)
%
% This takes the circular nature of angles into account.

if length(m)==1 % scalar
  x = x(:).'; % make row vector
end
[d,N] = size(x);
m = m(:);
M = m*ones(1,N); % replicate the mean across columns
denom = (2*pi)^(d/2)*sqrt(abs(det(C)));
dx = dist_pbc(x,M).';
mahal = sum((dx/C).*dx,2);
if any(mahal<0)
  warning('mahal < 0 => C is not positive semidefinite.')
end
p = exp(-0.5*mahal) / (denom+eps);

end

%===============================================================================
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
  n = sum(delta(:,t));
  if n~=0
    delta(:,t) = delta(:,t)/n;
  end
  scale(t) = 1/n;
end
psi(:,t) = 0; % arbitrary value, since there is no predecessor to t=1
for t=2:T
  for j=1:Q
    [delta(j,t), psi(j,t)] = max(delta(:,t-1) .* transmat(:,j));
    delta(j,t) = delta(j,t) * obslik(j,t);
  end
  if scaled
    n = sum(delta(:,t));
    if n~=0
      delta(:,t) = delta(:,t)/n;
    end
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

if false
  if scaled
    loglik = -sum(log(scale));
    %loglik = prob_path(prior, transmat, obslik, path);
  else
    loglik = log(p);
  end
end

end


%===============================================================================
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

[x, ~, ~] = mdhmm_lbfgsb(fcn, zeros(nStates*nStates, 1), Inf(nStates*nStates, 1), opts);
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