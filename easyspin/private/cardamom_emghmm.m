function [eqDistr, transmat, mu, Sigma, mixmat] = ...
     cardamom_emghmm(data, eqDistr, transmat, mu, Sigma, mixmat, verbose)
% LEARN_MHMM Compute the ML parameters of an HMM with (mixtures of) Gaussians output using EM.
% [prior, transmat, mu, sigma, mixmat] = learn_mhmm(data, ...
%   prior0, transmat0, mu0, sigma0, mixmat0, ...) 
%
% Notation: Q(t) = hidden state, Y(t) = observation, M(t) = mixture variable
%
% INPUTS:
% data{ex}(:,t) or data(:,t,ex) if all sequences have the same length
% prior(i) = Pr(Q(1) = i), 
% transmat(i,j) = Pr(Q(t+1)=j | Q(t)=i)
% mu(:,j,k) = E[Y(t) | Q(t)=j, M(t)=k ]
% Sigma(:,:,j,k) = Cov[Y(t) | Q(t)=j, M(t)=k]
% mixmat(j,k) = Pr(M(t)=k | Q(t)=j) : set to [] or ones(Q,1) if only one mixture component
%
% Optional parameters may be passed as 'param_name', param_value pairs.
% Parameter names are shown below; default values in [] - if none, argument is mandatory.
%
% 'max_iter' - max number of EM iterations [10]
% 'thresh' - convergence threshold [1e-4]
% 'verbose' - if 1, print out loglik at every iteration [1]

iterMax = 100;
thresh = 1e-4;
cov_type = 'full';
  
previous_loglik = -inf;
logLik = 0;
converged = 0;
iter = 1;
logLikIter = [];

if ~iscell(data)
  data = num2cell(data, [1 2]); % each elt of the 3rd dim gets its own cell
end

nTraj = length(data);
nSteps = size(data{1},1);
nStates = length(eqDistr);

while (iter <= iterMax) && ~converged
  % E step
%   [logLik, transCounts, exp_num_visits1, postmix, m, ip, op] = ...
%       ess_mhmm(eqDistr, transmat, mu, Sigma, data);
  transCounts = zeros(nStates,nStates);
  exp_num_visits1 = zeros(nStates,1);
  gamma_summed = zeros(nStates,1);

  logLik = 0;
  for iTraj=1:nTraj
    obs = data{iTraj};
    
    nSteps = size(obs,1);
    m = zeros(nSteps,nStates,1);
    op = zeros(nSteps,nSteps,nStates,1);
    ip = zeros(nStates,1);
    
    B = mixgauss_prob(obs, mu, Sigma);
    fwd_only = false;
    [alpha, beta, gamma, current_loglik, xi_summed] = fwdback(eqDistr, transmat, B, fwd_only);
    logLik = logLik +  current_loglik;

    transCounts = transCounts + xi_summed; % sum(xi,3);
    exp_num_visits1 = exp_num_visits1 + gamma(:,1);

    gamma_summed = gamma_summed + sum(gamma,2);
    for iState=1:nStates
      w = gamma(iState,:); % w(t) = w(i,k,t,l)
      wobs = obs .* repmat(w, [nSteps 1]); % wobs(:,t) = w(t) * obs(:,t)
      m(:,iState) = m(:,iState) + sum(wobs, 2); % m(:) = sum_t w(t) obs(:,t)
      op(:,:,iState) = op(:,:,iState) + wobs * obs'; % op(:,:) = sum_t w(t) * obs(:,t) * obs(:,t)'
      ip(iState,1) = ip(iState,1) + sum(sum(wobs .* obs, 2)); % ip = sum_t w(t) * obs(:,t)' * obs(:,t)
    end
  end
  
  
  % M step
  [transmat, eqDistr, ~] = msmtransitionmatrix(transCounts, 1000);
%   prior = normalise(exp_num_visits1);
%   transmat = mk_stochastic(exp_num_trans);
  [mu, Sigma] = mixgauss_Mstep(gamma_summed, m, op, ip);
  
  % Check convergence
  if verbose, fprintf(1, '    iteration %d, loglik = %f\n', iter, logLik); end
  iter =  iter + 1;
  converged = em_converged(logLik, previous_loglik, thresh);
  previous_loglik = logLik;
  logLikIter = [logLikIter logLik];
end

end


% Helper functions
% -------------------------------------------------------------------------

function [loglik, transCounts, exp_num_visits1, postmix, m, ip, op] = ...
    ess_mhmm(eqDistr, transmat, mu, Sigma, data)
% ESS_MHMM Compute the Expected Sufficient Statistics for a MOG Hidden Markov Model.
%
% Outputs:
% exp_num_trans(i,j)   = sum_l sum_{t=2}^T Pr(Q(t-1) = i, Q(t) = j| Obs(l))
% exp_num_visits1(i)   = sum_l Pr(Q(1)=i | Obs(l))
%
% Let w(i,k,t,l) = P(Q(t)=i, M(t)=k | Obs(l))
% where Obs(l) = Obs(:,:,l) = O_1 .. O_T for sequence l
% Then 
% postmix(i,k) = sum_l sum_t w(i,k,t,l) (posterior mixing weights/ responsibilities)
% m(:,i,k)   = sum_l sum_t w(i,k,t,l) * Obs(:,t,l)
% ip(i,k) = sum_l sum_t w(i,k,t,l) * Obs(:,t,l)' * Obs(:,t,l)
% op(:,:,i,k) = sum_l sum_t w(i,k,t,l) * Obs(:,t,l) * Obs(:,t,l)'


verbose = 0;

%[O T numex] = size(data);
nTraj = length(data);
nSteps = size(data{1},1);
nStates = length(eqDistr);
transCounts = zeros(nStates,nStates);
exp_num_visits1 = zeros(nStates,1);
postmix = zeros(nStates,1);
m = zeros(nSteps,nStates,1);
op = zeros(nSteps,nSteps,nStates,1);
ip = zeros(nStates,1);

loglik = 0;
for iTraj=1:nTraj
  obs = data{iTraj};
  B = mixgauss_prob(obs, mu, Sigma);
  fwd_only = false;
  [alpha, beta, gamma, current_loglik, xi_summed] = fwdback(eqDistr, transmat, B, fwd_only);
  loglik = loglik +  current_loglik;

  transCounts = transCounts + xi_summed; % sum(xi,3);
  exp_num_visits1 = exp_num_visits1 + gamma(:,1);
  
  postmix = postmix + sum(gamma,2);
  for iState=1:nStates
    w = gamma(iState,:); % w(t) = w(i,k,t,l)
    wobs = obs .* repmat(w, [nSteps 1]); % wobs(:,t) = w(t) * obs(:,t)
    m(:,iState) = m(:,iState) + sum(wobs, 2); % m(:) = sum_t w(t) obs(:,t)
    op(:,:,iState) = op(:,:,iState) + wobs * obs'; % op(:,:) = sum_t w(t) * obs(:,t) * obs(:,t)'
    ip(iState,1) = ip(iState,1) + sum(sum(wobs .* obs, 2)); % ip = sum_t w(t) * obs(:,t)' * obs(:,t)
  end
end

end

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

function p = gaussian_prob(x, m, C, use_log)
% GAUSSIAN_PROB Evaluate a multivariate Gaussian density.
% p = gaussian_prob(X, m, C)
% p(i) = N(X(:,i), m, C) where C = covariance matrix and each COLUMN of x is a datavector

% p = gaussian_prob(X, m, C, 1) returns log N(X(:,i), m, C) (to prevents underflow).
%
% If X has size dxN, then p has size Nx1, where N = number of examples

if nargin < 4, use_log = 0; end

if length(m)==1 % scalar
  x = x(:)';
end
[d N] = size(x);
%assert(length(m)==d); % slow
m = m(:);
M = m*ones(1,N); % replicate the mean across columns
denom = (2*pi)^(d/2)*sqrt(abs(det(C)));
mahal = sum(((x-M)'*inv(C)).*(x-M)',2);   % Chris Bregler's trick
if any(mahal<0)
  warning('mahal < 0 => C is not psd')
end
if use_log
  p = -0.5*mahal - log(denom);
else
  p = exp(-0.5*mahal) / (denom+eps);
end

end

function [mu, Sigma] = mixgauss_Mstep(w, Y, YY, YTY, varargin)
% MSTEP_COND_GAUSS Compute MLEs for mixture of Gaussians given expected sufficient statistics
% function [mu, Sigma] = Mstep_cond_gauss(w, Y, YY, YTY, varargin)
%
% We assume P(Y|Q=i) = N(Y; mu_i, Sigma_i)
% and w(i,t) = p(Q(t)=i|y(t)) = posterior responsibility
% See www.ai.mit.edu/~murphyk/Papers/learncg.pdf.
%
% INPUTS:
% w(i) = sum_t w(i,t) = responsibilities for each mixture component
%  If there is only one mixture component (i.e., Q does not exist),
%  then w(i) = N = nsamples,  and 
%  all references to i can be replaced by 1.
% YY(:,:,i) = sum_t w(i,t) y(:,t) y(:,t)' = weighted outer product
% Y(:,i) = sum_t w(i,t) y(:,t) = weighted observations
% YTY(i) = sum_t w(i,t) y(:,t)' y(:,t) = weighted inner product
%   You only need to pass in YTY if Sigma is to be estimated as spherical.
%
% Optional parameters may be passed as 'param_name', param_value pairs.
% Parameter names are shown below; default values in [] - if none, argument is mandatory.
%
% 'cov_type' - 'full', 'diag' or 'spherical' ['full']
% 'tied_cov' - 1 (Sigma) or 0 (Sigma_i) [0]
% 'clamped_cov' - pass in clamped value, or [] if unclamped [ [] ]
% 'clamped_mean' - pass in clamped value, or [] if unclamped [ [] ]
% 'cov_prior' - Lambda_i, added to YY(:,:,i) [0.01*eye(d,d,Q)]
%
% If covariance is tied, Sigma has size d*d.
% But diagonal and spherical covariances are represented in full size.

[nDims, nStates] = size(Y);
N = sum(w);
cov_prior = repmat(0.01*eye(nDims,nDims), [1 1 nStates]);
%YY = reshape(YY, [Ysz Ysz Q]) + cov_prior; % regularize the scatter matrix
YY = reshape(YY, [nDims nDims nStates]);

% Set any zero weights to one before dividing
% This is valid because w(i)=0 => Y(:,i)=0, etc
w = w + (w==0);
		    
% eqn 6
%mu = Y ./ repmat(w(:)', [Ysz 1]);% Y may have a funny size
mu = zeros(nDims, nStates);
for i=1:nStates
  mu(:,i) = Y(:,i) / w(i);
end

Sigma = zeros(nDims,nDims,nStates);
for i=1:nStates
  % eqn 12
  SS = YY(:,:,i)/w(i)  - mu(:,i)*mu(:,i)';
  Sigma(:,:,i) = SS;
end

Sigma = Sigma + cov_prior;

end

function B = mixgauss_prob(data, mu, Sigma)
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
%
% OUTPUT:
% B(i,t) = Pr(y(t) | iState(t)=i) 


nStates = size(mu,2);
[nDims, nSteps] = size(data);
  
B = zeros(nStates,nSteps);
% for iState=1:nStates
%   % D(m,t) = sq dist between data(:,t) and mu(:,j,m)
%   if ~isposdef(Sigma(:,:,iState))
%     Sigma(:,:,iState) = Sigma(:,:,iState) + 1e-3*eye(size(Sigma,1));
%   end
%   if isposdef(Sigma(:,:,iState))
%     D = sqdist(data, permute(mu(:,iState,:), [1 3 2]), inv(Sigma(:,:,iState)))';
%     logB2 = -(nDims/2)*log(2*pi) - 0.5*logdet(Sigma(:,:,iState)) - 0.5*D;
% %       logB2 = -(d/2)*log(2*pi) - 0.5*log(det(Sigma(:,:,j))) - 0.5*D;
%     B(iState,:) = exp(logB2);
%   else
%     error('mixgauss_prob: Sigma(:,:,q=%d) not psd\n', iState);
%   end
% end
  
% else % general case
%   B = zeros(nStates,nSteps);
  for iState=1:nStates
    B(iState,:) = gaussian_prob(data, mu(:,iState), Sigma(:,:,iState));
  end
% end

end



function [converged, decrease] = em_converged(loglik, previous_loglik, threshold, check_increased)
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

if nargin < 3, threshold = 1e-4; end
if nargin < 4, check_increased = 1; end

converged = 0;
decrease = 0;

if check_increased
  if loglik - previous_loglik < -1e-3 % allow for a little imprecision
    fprintf(1, '******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik);
    decrease = 1;
converged = 0;
return;
  end
end

delta_loglik = abs(loglik - previous_loglik);
avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2;
if (delta_loglik / avg_loglik) < threshold, converged = 1; end

end

function b = isposdef(a)
% ISPOSDEF   Test for positive definite matrix.
%    ISPOSDEF(A) returns 1 if A is positive definite, 0 otherwise.
%    Using chol is much more efficient than computing eigenvectors.

%  From Tom Minka's lightspeed toolbox

[R,p] = chol(a);
b = (p == 0);

end

function y = logdet(A)
% log(det(A)) where A is positive-definite.
% This is faster and more stable than using log(det(A)).

%  From Tom Minka's lightspeed toolbox

U = chol(A);
y = 2*sum(log(diag(U)));

end

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

function [varargout] = process_options(args, varargin)
% PROCESS_OPTIONS - Processes options passed to a Matlab function.
%                   This function provides a simple means of
%                   parsing attribute-value options.  Each option is
%                   named by a unique string and is given a default
%                   value.
%
% Usage:  [var1, var2, ..., varn[, unused]] = ...
%           process_options(args, ...
%                           str1, def1, str2, def2, ..., strn, defn)
%
% Arguments:   
%            args            - a cell array of input arguments, such
%                              as that provided by VARARGIN.  Its contents
%                              should alternate between strings and
%                              values.
%            str1, ..., strn - Strings that are associated with a 
%                              particular variable
%            def1, ..., defn - Default values returned if no option
%                              is supplied
%
% Returns:
%            var1, ..., varn - values to be assigned to variables
%            unused          - an optional cell array of those 
%                              string-value pairs that were unused;
%                              if this is not supplied, then a
%                              warning will be issued for each
%                              option in args that lacked a match.
%
% Examples:
%
% Suppose we wish to define a Matlab function 'func' that has
% required parameters x and y, and optional arguments 'u' and 'v'.
% With the definition
%
%   function y = func(x, y, varargin)
%
%     [u, v] = process_options(varargin, 'u', 0, 'v', 1);
%
% calling func(0, 1, 'v', 2) will assign 0 to x, 1 to y, 0 to u, and 2
% to v.  The parameter names are insensitive to case; calling 
% func(0, 1, 'V', 2) has the same effect.  The function call
% 
%   func(0, 1, 'u', 5, 'z', 2);
%
% will result in u having the value 5 and v having value 1, but
% will issue a warning that the 'z' option has not been used.  On
% the other hand, if func is defined as
%
%   function y = func(x, y, varargin)
%
%     [u, v, unused_args] = process_options(varargin, 'u', 0, 'v', 1);
%
% then the call func(0, 1, 'u', 5, 'z', 2) will yield no warning,
% and unused_args will have the value {'z', 2}.  This behaviour is
% useful for functions with options that invoke other functions
% with options; all options can be passed to the outer function and
% its unprocessed arguments can be passed to the inner function.

% Copyright (C) 2002 Mark A. Paskin
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
% USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the number of input arguments
n = length(varargin);
if (mod(n, 2))
  error('Each option must be a string/value pair.');
end

% Check the number of supplied output arguments
if (nargout < (n / 2))
  error('Insufficient number of output arguments given');
elseif (nargout == (n / 2))
  warn = 1;
  nout = n / 2;
else
  warn = 0;
  nout = n / 2 + 1;
end

% Set outputs to be defaults
varargout = cell(1, nout);
for i=2:2:n
  varargout{i/2} = varargin{i};
end

% Now process all arguments
nunused = 0;
for i=1:2:length(args)
  found = 0;
  for j=1:2:n
    if strcmpi(args{i}, varargin{j})
      varargout{(j + 1)/2} = args{i + 1};
      found = 1;
      break;
    end
  end
  if (~found)
    if (warn)
      warning(sprintf('Option ''%s'' not used.', args{i}));
      args{i}
    else
      nunused = nunused + 1;
      unused{2 * nunused - 1} = args{i};
      unused{2 * nunused} = args{i + 1};
    end
  end
end

% Assign the unused arguments
if (~warn)
  if (nunused)
    varargout{nout} = unused;
  else
    varargout{nout} = cell(0);
  end
end

end

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

  if isempty(A) || isempty(p)
    error('sqdist: empty matrices');
  end
  Ap = A*p;
  Aq = A*q;
  pmag = sum(p .* Ap, 1);
  qmag = sum(q .* Aq, 1);
  m = repmat(qmag, pn, 1) + repmat(pmag', 1, qn) - 2*p'*Aq;
  
end

end


function c = msmcountmatrix(indexOfCluster, tau, nstate)
% msmcountmatrix
% calculate transition count matrix from a set of binned trajectory data
%
% Syntax
%# c = msmcountmatrix(indexOfCluster);
%# c = msmcountmatrix(indexOfCluster, tau);
%# c = msmcountmatrix(indexOfCluster, tau, nstate);
%# c = msmcountmatrix(indexOfCluster, [], nstate);
%
% Description
% calculate count matrix of transition from state i to state j during time step tau
%
% Adapted from Yasuhiro Matsunaga's mdtoolbox
% 

% setup
if ~iscell(indexOfCluster)
  indexOfCluster_noncell = indexOfCluster;
  clear indexOfCluster;
  indexOfCluster{1} = indexOfCluster_noncell;
  clear indexOfCluster_noncell;
end
ntrj = numel(indexOfCluster);

if ~exist('nstate', 'var') || isempty(nstate)
  nstate = max(cellfun(@(x) max(x), indexOfCluster));
  disp(sprintf('Message: nstate = %d is used.', nstate));
end

if ~exist('tau', 'var') || isempty(tau)
  tau = 1;
  disp('Message: tau = 1 is used.');
end

% count transitions
c = sparse(nstate, nstate);

for itrj = 1:ntrj
  nframe = numel(indexOfCluster{itrj});

  index_from = 1:(nframe-tau);
  index_to   = (1+tau):nframe;
  indexOfCluster_from = indexOfCluster{itrj}(index_from);
  indexOfCluster_to   = indexOfCluster{itrj}(index_to);

  % ignore invalid indices
  nframe = numel(indexOfCluster_from);
  s = ones(nframe, 1);

  id = (indexOfCluster_from <= 0);
  s(id) = 0;
  indexOfCluster_from(id) = 1;

  id = (indexOfCluster_to   <= 0);
  s(id) = 0;
  indexOfCluster_to(id)   = 1;

  id = isnan(indexOfCluster_from);
  s(id) = 0;
  indexOfCluster_from(id) = 1;

  id = isnan(indexOfCluster_to);
  s(id) = 0;
  indexOfCluster_to(id)   = 1;

  % calc count matrix
  % count transitions and make count matrix C_ij by using a sparse
  % matrix
  c_itrj = sparse(indexOfCluster_from, indexOfCluster_to, s, nstate, nstate);
  c = c + c_itrj;
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

[x, f, info] = cardamom_lbfgsb(fcn, zeros(nstate*nstate, 1), Inf(nstate*nstate, 1), opts);
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
