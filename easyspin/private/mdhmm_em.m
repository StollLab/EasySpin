% mdhmm_em  Calculate the maximum-likelihood parameters for a multivariate
%           Gaussian HMM using the Baum-Welch algorithm
%
%   function [eqDistr, TransProb, mu, Sigma] = ...
%     mdhmm_em(data, initDistr, TransProb, mu, Sigma, verbosity)
%
%
%     data           numeric matrix, size = (nDims,nSteps,nTraj)
%                    input data
%
%     initdistr      numeric vector, size = (1,nStates)
%                    initial starting state probability distribution
%
%     TransProb      numeric, size = (nStates,nStates)
%                    initial transition probability matrix describing 
%                    inter-state dynamics for a given time step
%
%     mu             numeric matrix, size = (nDims,nStates)
%                    initial centers of Gaussians assigned to states
%
%     Sigma          numeric matrix, size = (nDims,nDims,nStates)
%                    initial covariance matrices of Gaussians assigned to
%                    states
%
%     verbosity      0 or 1
%                    if set to 1, display information regarding the
%                    progress
%    
%   Output:
%
%     TransProb      numeric, size = (nStates,nStates)
%                    transition probability matrix describing inter-state 
%                    dynamics for a given time step
%
%     mu             numeric matrix, size = (nDims,nStates)
%                    centers of Gaussians assigned to states
%
%     Sigma          numeric matrix, size = (nDims,nDims,nStates)
%                    covariance matrices of Gaussians assigned to states
%
%     eqdistr        numeric vector, size = (1,nStates)
%                    equilibrium state probability distribution

% Code adapted from Kevin Murphy's HMM toolbox

function [logLikIter, eqDistr, TransProb, mu, Sigma] = ...
     mdhmm_em(data, initDistr, TransProb, mu, Sigma, verbosity)

iterMax = 200;
thresh = 1e-4;
  
previous_loglik = -inf;
converged = false;
iter = 1;
logLikIter = [];

if ~iscell(data)
  data = num2cell(data, [1 2]); % each elt of the 3rd dim gets its own cell
end

nTraj = length(data);
nStates = length(initDistr);

while (iter <= iterMax) && ~converged
  % E step
  transCounts = zeros(nStates,nStates);
  expNumVisits1 = zeros(nStates,1);
  weightsSummed = zeros(nStates,1);

  logLik = 0;
  for iTraj = 1:nTraj
    obs = data{iTraj};
    [nDims,nSteps] = size(obs);
    
    muUpdater = zeros(nDims,nStates);
    
    B = zeros(nStates,nSteps);
    for iState=1:nStates
      B(iState,:) = gaussian_prob_pbc(obs, mu(:,iState), Sigma(:,:,iState), 2*pi);
    end
    fwd_only = false;
    [~, ~, gamma, current_loglik, xi_summed] = fwdback(initDistr, TransProb, B, fwd_only);
    logLik = logLik +  current_loglik;

    transCounts = transCounts + xi_summed; % sum(xi,3);
    expNumVisits1 = expNumVisits1 + gamma(:,1);

    weightsSummed = weightsSummed + sum(gamma,2);
    for iState=1:nStates
      weights = gamma(iState,:);
      distances = dist_pbc(obs - mu(:,iState), 2*pi);
      muUpdater(:,iState) = muUpdater(:,iState) + sum(weights .* distances, 2);
    end
  end
  
  
  % M step
  [TransProb, eqDistr, ~] = msmtransitionmatrix(transCounts, 1000);
  initDistr = normalise(expNumVisits1);
%   transmat = mk_stochastic(exp_num_trans);
  cov_prior = repmat(0.01*eye(nDims,nDims), [1 1 nStates]);

  % Set any zero weights to one before dividing
  % This is valid because w(i)=0 => Y(:,i)=0, etc
  weightsSummed = weightsSummed + (weightsSummed==0);

  % Update means
  for iState=1:nStates
    mu(:,iState) = mu(:,iState) + muUpdater(:,iState) / weightsSummed(iState);
  end
  mu(mu>pi) = mu(mu>pi) - 2*pi;
  mu(mu<-pi) = mu(mu<-pi) + 2*pi;

  % Update covariance matrices using updated means
  Sigma = zeros(nDims,nDims,nStates);
  for iState=1:nStates
    weights = gamma(iState,:);
    distances = dist_pbc(obs - mu(:,iState), 2*pi);
    SigmaUpdater = weights .* distances * distances';
    Sigma(:,:,iState) = SigmaUpdater/weightsSummed(iState);
  end

  % Ensure that no covariance matrix is too small
  Sigma = Sigma + cov_prior;
  
  % Check convergence
  if verbosity, fprintf(1, '    iteration %d, loglik = %f\n', iter, logLik); end
  iter =  iter + 1;
  converged = em_converged(logLik, previous_loglik, thresh, verbosity);
  previous_loglik = logLik;
  logLikIter = [logLikIter logLik];
end

end


% Helper functions
% -------------------------------------------------------------------------

function [alpha, beta, gamma, loglik, xi_summed] = fwdback(prior, ...
   transmat, obslik, fwd_only)
% FWDBACK Compute the posterior probs. in an HMM using the forwards backwards algo.
%
% [alpha, beta, gamma, loglik, xi, gamma2] = fwdback(init_state_distrib, transmat, obslik, ...)
%
% Notation:
% Y(t) = observation, Q(t) = hidden state, M(t) = mixture variable (for MOG outputs)
% A(t) = discrete input (action) (for POMDP models)
%
% INPUT:
% init_state_distrib(i) = Pr(Q(1) = i)
% transmat(i,j) = Pr(Q(t) = j | Q(t-1)=i)
%  or transmat{a}(i,j) = Pr(Q(t) = j | Q(t-1)=i, A(t-1)=a) if there are discrete inputs
% obslik(i,t) = Pr(Y(t)| Q(t)=i)
%
%
% Optional arguments:
% 'fwd_only' - if 1, only do a forwards pass and set beta=[], gamma2=[]  [0]
%
% OUTPUTS:
% alpha(i,t) = p(Q(t)=i | y(1:t)) (or p(Q(t)=i, y(1:t)) if scaled=0)
% beta(i,t) = p(y(t+1:T) | Q(t)=i)*p(y(t+1:T)|y(1:t)) (or p(y(t+1:T) | Q(t)=i) if scaled=0)
% gamma(i,t) = p(Q(t)=i | y(1:T))
% loglik = log p(y(1:T))
% xi_summed(i,j) = sum_{t=}^{T-1} xi(i,j,t)  - changed made by Herbert Jaeger
%
% If fwd_only = 1, these become
% alpha(i,t) = p(Q(t)=i | y(1:t))
% beta = []
% gamma(i,t) = p(Q(t)=i | y(1:t))
% xi(i,j,t-1)  = p(Q(t-1)=i, Q(t)=j | y(1:t))

[nStates, nSteps] = size(obslik);

scale = ones(1,nSteps);

% scale(t) = Pr(O(t) | O(1:t-1)) = 1/c(t) as defined by Rabiner (1989).
% Hence prod_t scale(t) = Pr(O(1)) Pr(O(2)|O(1)) Pr(O(3) | O(1:2)) ... = Pr(O(1), ... ,O(T))
% or log P = sum_t log scale(t).
% Rabiner suggests multiplying beta(t) by scale(t), but we can instead
% normalise beta(t) - the constants will cancel when we compute gamma.

loglik = 0;

alpha = zeros(nStates,nSteps);
gamma = zeros(nStates,nSteps);
xi_summed = zeros(nStates,nStates);

%%%%%%%%% Forwards %%%%%%%%%%

[alpha(:,1), scale(1)] = normalise(prior(:) .* obslik(:,1));

for iStep=2:nSteps
 [alpha(:,iStep), scale(iStep)] = normalise(transmat' * alpha(:,iStep-1) .* obslik(:,iStep));
 if fwd_only  % useful for online EM
   xi_summed = xi_summed + normalise((alpha(:,iStep-1) * obslik(:,iStep)') .* transmat);
 end
end

if any(scale==0)
 loglik = -inf;
else
 loglik = sum(log(scale));
end

if fwd_only
 gamma = alpha;
 beta = [];
 return;
end

%%%%%%%%% Backwards %%%%%%%%%%

beta = zeros(nStates,nSteps);

beta(:,end) = ones(nStates,1);
gamma(:,end) = normalise(alpha(:,end) .* beta(:,end));

for iStep=nSteps-1:-1:1
 b = beta(:,iStep+1) .* obslik(:,iStep+1);
 beta(:,iStep) = normalise(transmat * b);
 gamma(:,iStep) = normalise(alpha(:,iStep) .* beta(:,iStep));
 xi_summed = xi_summed + normalise((transmat .* (alpha(:,iStep) * b')));
end

end

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

function dist = dist_pbc(dist,W)

w = W/2;

% dist = w - abs(pi - abs(dist));

idx1 = dist > w;
idx2 = dist < -w;

dist(idx1) = dist(idx1) - W;
dist(idx2) = dist(idx2) + W;

end

function [converged, decrease] = em_converged(loglik, previous_loglik, threshold, verbosity)
% EM_CONVERGED Has EM converged?
% [converged, decrease] = em_converged(loglik, previous_loglik, threshold)
%
% We have converged if the slope of the log-likelihood function falls below 'threshold', 
% i.e., |f(t) - f(t-1)| / avg < threshold,
% where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log lik at iteration t.
% 'threshold' defaults to 1e-4.
%
% This stopping criterion is from Numerical Recipes in C p423
%
% If we are doing MAP estimation (using priors), the likelihood can decrase,
% even though the mode of the posterior is increasing.

converged = false;
decrease = false;

if loglik - previous_loglik < -1e-3 % allow for a little imprecision
  if verbosity
    fprintf(1, '      likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik);
  end
  decrease = true;
  converged = false;
  return
end

delta_loglik = abs(loglik - previous_loglik);
avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2;
if (delta_loglik / avg_loglik) < threshold, converged = true; end

end

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

function [t, pi_i, x] = msmtransitionmatrix(c, maxiteration)
% msmtransitionmatrix
% estimate transition probability matrix from count matrix
%
% Syntax
%# [t, pi_i] = msmtransitionmatrix(c);
%# [t, pi_i] = msmtransitionmatrix(c, maxiteration);
%
% Description
% this routines uses the reversible maximum likelihood estimator
%
% Adapted from Yasuhiro Matsunaga's mdtoolbox
% 

% setup
if issparse(c)
  c = full(c);
end

nstate = size(c, 1);

c_sym  = c + c';
x      = c_sym;

c_i    = sum(c, 2);
x_i    = sum(x, 2);

if ~exist('maxiteration', 'var')
  maxiteration = 1000;
end

% optimization by L-BFGS-B
fcn = @(x) myfunc_column(x, c, c_i, nstate);
opts.x0 = x(:);
opts.maxIts = maxiteration;
opts.maxTotalIts = 50000;
%opts.factr = 1e5;
%opts.pgtol = 1e-7;

[x, f, info] = mdhmm_lbfgsb(fcn, zeros(nstate*nstate, 1), Inf(nstate*nstate, 1), opts);
x = reshape(x, nstate, nstate);

x_i = sum(x, 2);
t = bsxfun(@rdivide, x, x_i);
t(isnan(t)) = 0;
pi_i = x_i./sum(x_i);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, g] = myfunc_column(x, c, c_i, nstate)
x = reshape(x, nstate, nstate);
[f, g] = myfunc_matrix(x, c, c_i);
g = g(:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
