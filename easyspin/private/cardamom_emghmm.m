function [LL, prior, transmat, mu, Sigma, mixmat] = ...
     mhmm_em(data, prior, transmat, mu, Sigma, mixmat, varargin);
% LEARN_MHMM Compute the ML parameters of an HMM with (mixtures of) Gaussians output using EM.
% [ll_trace, prior, transmat, mu, sigma, mixmat] = learn_mhmm(data, ...
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
% 'cov_type' - 'full', 'diag' or 'spherical' ['full']
%
% To clamp some of the parameters, so learning does not change them:
% 'adj_prior' - if 0, do not change prior [1]
% 'adj_trans' - if 0, do not change transmat [1]
% 'adj_mix' - if 0, do not change mixmat [1]
% 'adj_mu' - if 0, do not change mu [1]
% 'adj_Sigma' - if 0, do not change Sigma [1]
%
% If the number of mixture components differs depending on Q, just set  the trailing
% entries of mixmat to 0, e.g., 2 components if Q=1, 3 components if Q=2,
% then set mixmat(1,3)=0. In this case, B2(1,3,:)=1.0.

if ~isempty(varargin) & ~ischar(varargin{1}) % catch old syntax
  error('optional arguments should be passed as string/value pairs')
end

[max_iter, thresh, verbose, cov_type,  adj_prior, adj_trans, adj_mix, adj_mu, adj_Sigma] = ...
    process_options(varargin, 'max_iter', 10, 'thresh', 1e-4, 'verbose', 1, ...
		    'cov_type', 'full', 'adj_prior', 1, 'adj_trans', 1, 'adj_mix', 1, ...
		    'adj_mu', 1, 'adj_Sigma', 1);
  
previous_loglik = -inf;
loglik = 0;
converged = 0;
num_iter = 1;
LL = [];

if ~iscell(data)
  data = num2cell(data, [1 2]); % each elt of the 3rd dim gets its own cell
end
numex = length(data);


O = size(data{1},1);
Q = length(prior);
if isempty(mixmat)
  mixmat = ones(Q,1);
end
M = size(mixmat,2);
if M == 1
  adj_mix = 0;
end

while (num_iter <= max_iter) & ~converged
  % E step
  [loglik, exp_num_trans, exp_num_visits1, postmix, m, ip, op] = ...
      ess_mhmm(prior, transmat, mixmat, mu, Sigma, data);
  
  
  % M step
  [transmat, prior, ~] = msmtransitionmatrix(exp_num_trans, 1000);
%   if adj_prior
%     prior = normalise(exp_num_visits1);
%   end
%   if adj_trans 
%     transmat = mk_stochastic(exp_num_trans);
%   end
%   if adj_mix
%     mixmat = mk_stochastic(postmix);
%   end
  if adj_mu | adj_Sigma
    [mu2, Sigma2] = mixgauss_Mstep(postmix, m, op, ip, 'cov_type', cov_type);
    if adj_mu
      mu = reshape(mu2, [O Q M]);
    end
    if adj_Sigma
      Sigma = reshape(Sigma2, [O O Q M]);
    end
  end
  
  if verbose, fprintf(1, 'iteration %d, loglik = %f\n', num_iter, loglik); end
  num_iter =  num_iter + 1;
  converged = em_converged(loglik, previous_loglik, thresh);
  previous_loglik = loglik;
  LL = [LL loglik];
end

end


% Helper functions
% -------------------------------------------------------------------------

function [loglik, exp_num_trans, exp_num_visits1, postmix, m, ip, op] = ...
    ess_mhmm(prior, transmat, mixmat, mu, Sigma, data)
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
numex = length(data);
O = size(data{1},1);
Q = length(prior);
M = size(mixmat,2);
exp_num_trans = zeros(Q,Q);
exp_num_visits1 = zeros(Q,1);
postmix = zeros(Q,M);
m = zeros(O,Q,M);
op = zeros(O,O,Q,M);
ip = zeros(Q,M);

mix = (M>1);

loglik = 0;
if verbose, fprintf(1, 'forwards-backwards example # '); end
for ex=1:numex
  if verbose, fprintf(1, '%d ', ex); end
  %obs = data(:,:,ex);
  obs = data{ex};
  T = size(obs,2);
  if mix
    [B, B2] = mixgauss_prob(obs, mu, Sigma, mixmat);
    [alpha, beta, gamma,  current_loglik, xi_summed, gamma2] = ...
	fwdback(prior, transmat, B, 'obslik2', B2, 'mixmat', mixmat);
  else
    B = mixgauss_prob(obs, mu, Sigma);
%     [alpha, beta, gamma,  current_loglik, xi_summed] = fwdback(prior, transmat, B);
    [alpha, beta, gamma,  current_loglik, xi_summed] = fwdback(prior, transmat, B, 'fwd_only', 1);
  end    
  loglik = loglik +  current_loglik; 
  if verbose, fprintf(1, 'll at ex %d = %f\n', ex, loglik); end

  exp_num_trans = exp_num_trans + xi_summed; % sum(xi,3);
  exp_num_visits1 = exp_num_visits1 + gamma(:,1);
  
  if mix
    postmix = postmix + sum(gamma2,3);
  else
    postmix = postmix + sum(gamma,2); 
    gamma2 = reshape(gamma, [Q 1 T]); % gamma2(i,m,t) = gamma(i,t)
  end
  for i=1:Q
    for k=1:M
      w = reshape(gamma2(i,k,:), [1 T]); % w(t) = w(i,k,t,l)
      wobs = obs .* repmat(w, [O 1]); % wobs(:,t) = w(t) * obs(:,t)
      m(:,i,k) = m(:,i,k) + sum(wobs, 2); % m(:) = sum_t w(t) obs(:,t)
      op(:,:,i,k) = op(:,:,i,k) + wobs * obs'; % op(:,:) = sum_t w(t) * obs(:,t) * obs(:,t)'
      ip(i,k) = ip(i,k) + sum(sum(wobs .* obs, 2)); % ip = sum_t w(t) * obs(:,t)' * obs(:,t)
    end
  end
end
if verbose, fprintf(1, '\n'); end

end

function [alpha, beta, gamma, loglik, xi_summed, gamma2] = fwdback(init_state_distrib, ...
   transmat, obslik, varargin)
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
%   (Compute obslik using eval_pdf_xxx on your data sequence first.)
%
% Optional parameters may be passed as 'param_name', param_value pairs.
% Parameter names are shown below; default values in [] - if none, argument is mandatory.
%
% For HMMs with MOG outputs: if you want to compute gamma2, you must specify
% 'obslik2' - obslik(i,j,t) = Pr(Y(t)| Q(t)=i,M(t)=j)  []
% 'mixmat' - mixmat(i,j) = Pr(M(t) = j | Q(t)=i)  []
%  or mixmat{t}(m,q) if not stationary
%
% For HMMs with discrete inputs:
% 'act' - act(t) = action performed at step t
%
% Optional arguments:
% 'fwd_only' - if 1, only do a forwards pass and set beta=[], gamma2=[]  [0]
% 'scaled' - if 1,  normalize alphas and betas to prevent underflow [1]
% 'maximize' - if 1, use max-product instead of sum-product [0]
%
% OUTPUTS:
% alpha(i,t) = p(Q(t)=i | y(1:t)) (or p(Q(t)=i, y(1:t)) if scaled=0)
% beta(i,t) = p(y(t+1:T) | Q(t)=i)*p(y(t+1:T)|y(1:t)) (or p(y(t+1:T) | Q(t)=i) if scaled=0)
% gamma(i,t) = p(Q(t)=i | y(1:T))
% loglik = log p(y(1:T))
% xi(i,j,t-1)  = p(Q(t-1)=i, Q(t)=j | y(1:T))  - NO LONGER COMPUTED
% xi_summed(i,j) = sum_{t=}^{T-1} xi(i,j,t)  - changed made by Herbert Jaeger
% gamma2(j,k,t) = p(Q(t)=j, M(t)=k | y(1:T)) (only for MOG  outputs)
%
% If fwd_only = 1, these become
% alpha(i,t) = p(Q(t)=i | y(1:t))
% beta = []
% gamma(i,t) = p(Q(t)=i | y(1:t))
% xi(i,j,t-1)  = p(Q(t-1)=i, Q(t)=j | y(1:t))
% gamma2 = []
%
% Note: we only compute xi if it is requested as a return argument, since it can be very large.
% Similarly, we only compute gamma2 on request (and if using MOG outputs).
%
% Examples:
%
% [alpha, beta, gamma, loglik] = fwdback(pi, A, multinomial_prob(sequence, B));
%
% [B, B2] = mixgauss_prob(data, mu, Sigma, mixmat);
% [alpha, beta, gamma, loglik, xi, gamma2] = fwdback(pi, A, B, 'obslik2', B2, 'mixmat', mixmat);

if 0 % nargout >= 5
  warning('this now returns sum_t xi(i,j,t) not xi(i,j,t)')
end

if nargout >= 5, compute_xi = 1; else compute_xi = 0; end
if nargout >= 6, compute_gamma2 = 1; else compute_gamma2 = 0; end

[obslik2, mixmat, fwd_only, scaled, act, maximize, compute_xi, compute_gamma2] = ...
   process_options(varargin, ...
       'obslik2', [], 'mixmat', [], ...
       'fwd_only', 0, 'scaled', 1, 'act', [], 'maximize', 0, ...
                   'compute_xi', compute_xi, 'compute_gamma2', compute_gamma2);

[Q T] = size(obslik);

if isempty(obslik2)
 compute_gamma2 = 0;
end

if isempty(act)
 act = ones(1,T);
 transmat = { transmat } ;
end

scale = ones(1,T);

% scale(t) = Pr(O(t) | O(1:t-1)) = 1/c(t) as defined by Rabiner (1989).
% Hence prod_t scale(t) = Pr(O(1)) Pr(O(2)|O(1)) Pr(O(3) | O(1:2)) ... = Pr(O(1), ... ,O(T))
% or log P = sum_t log scale(t).
% Rabiner suggests multiplying beta(t) by scale(t), but we can instead
% normalise beta(t) - the constants will cancel when we compute gamma.

loglik = 0;

alpha = zeros(Q,T);
gamma = zeros(Q,T);
if compute_xi
 xi_summed = zeros(Q,Q);
else
 xi_summed = [];
end

%%%%%%%%% Forwards %%%%%%%%%%

t = 1;
alpha(:,1) = init_state_distrib(:) .* obslik(:,t);
if scaled
 [alpha(:,t), scale(t)] = normalise(alpha(:,t));
end

for t=2:T
 trans = transmat{act(t-1)};
 if maximize
   m = max_mult(trans', alpha(:,t-1));
 else
   m = trans' * alpha(:,t-1);
 end
 alpha(:,t) = m(:) .* obslik(:,t);
 if scaled
   [alpha(:,t), scale(t)] = normalise(alpha(:,t));
 end
 if compute_xi & fwd_only  % useful for online EM
   xi_summed = xi_summed + normalise((alpha(:,t-1) * obslik(:,t)') .* trans);
 end
 %assert(approxeq(sum(alpha(:,t)),1))
end
if scaled
 if any(scale==0)
   loglik = -inf;
 else
   loglik = sum(log(scale));
 end
else
 loglik = log(sum(alpha(:,T)));
end

if fwd_only
 gamma = alpha;
 beta = [];
 gamma2 = [];
 return;
end

%%%%%%%%% Backwards %%%%%%%%%%

beta = zeros(Q,T);
if compute_gamma2
  if iscell(mixmat)
    M = size(mixmat{1},2);
  else
    M = size(mixmat, 2);
  end
 gamma2 = zeros(Q,M,T);
else
 gamma2 = [];
end

beta(:,T) = ones(Q,1);
gamma(:,T) = normalise(alpha(:,T) .* beta(:,T));
t=T;
if compute_gamma2
 denom = obslik(:,t) + (obslik(:,t)==0); % replace 0s with 1s before dividing
 if iscell(mixmat)
   gamma2(:,:,t) = obslik2(:,:,t) .* mixmat{t} .* repmat(gamma(:,t), [1 M]) ./ repmat(denom, [1 M]);
 else
   gamma2(:,:,t) = obslik2(:,:,t) .* mixmat .* repmat(gamma(:,t), [1 M]) ./ repmat(denom, [1 M]);
 end
end

for t=T-1:-1:1
 b = beta(:,t+1) .* obslik(:,t+1);
 trans = transmat{act(t)};
 if maximize
   B = repmat(b(:)', Q, 1);
   beta(:,t) = max(trans .* B, [], 2);
 else
   beta(:,t) = trans * b;
 end
 if scaled
   beta(:,t) = normalise(beta(:,t));
 end
 gamma(:,t) = normalise(alpha(:,t) .* beta(:,t));
 if compute_xi
   xi_summed = xi_summed + normalise((trans .* (alpha(:,t) * b')));
 end
 if compute_gamma2
   denom = obslik(:,t) + (obslik(:,t)==0); % replace 0s with 1s before dividing
   if iscell(mixmat)
     gamma2(:,:,t) = obslik2(:,:,t) .* mixmat{t} .* repmat(gamma(:,t), [1 M]) ./ repmat(denom,  [1 M]);
   else
     gamma2(:,:,t) = obslik2(:,:,t) .* mixmat .* repmat(gamma(:,t), [1 M]) ./ repmat(denom,  [1 M]);
   end
 end
end

% We now explain the equation for gamma2
% Let zt=y(1:t-1,t+1:T) be all observations except y(t)
% gamma2(Q,M,t) = P(Qt,Mt|yt,zt) = P(yt|Qt,Mt,zt) P(Qt,Mt|zt) / P(yt|zt)
%                = P(yt|Qt,Mt) P(Mt|Qt) P(Qt|zt) / P(yt|zt)
% Now gamma(Q,t) = P(Qt|yt,zt) = P(yt|Qt) P(Qt|zt) / P(yt|zt)
% hence
% P(Qt,Mt|yt,zt) = P(yt|Qt,Mt) P(Mt|Qt) [P(Qt|yt,zt) P(yt|zt) / P(yt|Qt)] / P(yt|zt)
%                = P(yt|Qt,Mt) P(Mt|Qt) P(Qt|yt,zt) / P(yt|Qt)

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

[cov_type, tied_cov,  clamped_cov, clamped_mean, cov_prior, other] = ...
    process_options(varargin,...
		    'cov_type', 'full', 'tied_cov', 0,  'clamped_cov', [], 'clamped_mean', [], ...
		    'cov_prior', []);

[Ysz Q] = size(Y);
N = sum(w);
if isempty(cov_prior)
  %cov_prior = zeros(Ysz, Ysz, Q);
  %for q=1:Q
  %  cov_prior(:,:,q) = 0.01*cov(Y(:,q)');
  %end
  cov_prior = repmat(0.01*eye(Ysz,Ysz), [1 1 Q]);
end
%YY = reshape(YY, [Ysz Ysz Q]) + cov_prior; % regularize the scatter matrix
YY = reshape(YY, [Ysz Ysz Q]);

% Set any zero weights to one before dividing
% This is valid because w(i)=0 => Y(:,i)=0, etc
w = w + (w==0);
		    
if ~isempty(clamped_mean)
  mu = clamped_mean;
else
  % eqn 6
  %mu = Y ./ repmat(w(:)', [Ysz 1]);% Y may have a funny size
  mu = zeros(Ysz, Q);
  for i=1:Q
    mu(:,i) = Y(:,i) / w(i);
  end
end

if ~isempty(clamped_cov)
  Sigma = clamped_cov;
  return;
end

if ~tied_cov
  Sigma = zeros(Ysz,Ysz,Q);
  for i=1:Q
    if cov_type(1) == 's'
      % eqn 17
      s2 = (1/Ysz)*( (YTY(i)/w(i)) - mu(:,i)'*mu(:,i) );
      Sigma(:,:,i) = s2 * eye(Ysz);
    else
      % eqn 12
      SS = YY(:,:,i)/w(i)  - mu(:,i)*mu(:,i)';
      if cov_type(1)=='d'
	SS = diag(diag(SS));
      end
      Sigma(:,:,i) = SS;
    end
  end
else % tied cov
  if cov_type(1) == 's'
    % eqn 19
    s2 = (1/(N*Ysz))*(sum(YTY,2) + sum(diag(mu'*mu) .* w));
    Sigma = s2*eye(Ysz);
  else
    SS = zeros(Ysz, Ysz);
    % eqn 15
    for i=1:Q % probably could vectorize this...
      SS = SS + YY(:,:,i)/N - mu(:,i)*mu(:,i)';
    end
    if cov_type(1) == 'd'
      Sigma = diag(diag(SS));
    else
      Sigma = SS;
    end
  end
end

if tied_cov
  Sigma =  repmat(Sigma, [1 1 Q]);
end
Sigma = Sigma + cov_prior;

end

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
    if ~isposdef(Sigma(:,:,j))
      Sigma(:,:,j) = Sigma(:,:,j) + 1e-3*eye(size(Sigma,1));
    end
    if isposdef(Sigma(:,:,j))
      D = sqdist(data, permute(mu(:,j,:), [1 3 2]), inv(Sigma(:,:,j)))';
      logB2 = -(d/2)*log(2*pi) - 0.5*logdet(Sigma(:,:,j)) - 0.5*D;
%       logB2 = -(d/2)*log(2*pi) - 0.5*log(det(Sigma(:,:,j))) - 0.5*D;
      B2(j,:,:) = exp(logB2);
    else
      error('mixgauss_prob: Sigma(:,:,q=%d) not psd\n', j);
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

function p = approxeq(a, b, tol, rel)
% APPROXEQ Are a and b approximately equal (to within a specified tolerance)?
% p = approxeq(a, b, thresh)
% 'tol' defaults to 1e-3.
% p(i) = 1 iff abs(a(i) - b(i)) < thresh
%
% p = approxeq(a, b, thresh, 1)
% p(i) = 1 iff abs(a(i)-b(i))/abs(a(i)) < thresh

if nargin < 3, tol = 1e-2; end
if nargin < 4, rel = 0; end

a = a(:);
b = b(:);
d = abs(a-b);
if rel
  p = ~any( (d ./ (abs(a)+eps)) > tol);
else
  p = ~any(d > tol);
end

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


function c = msmcountmatrix(indexOfCluster, tau, nstate)
%% msmcountmatrix
% calculate transition count matrix from a set of binned trajectory data
%
%% Syntax
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

%% setup
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

%% count transitions
c = sparse(nstate, nstate);

for itrj = 1:ntrj
  nframe = numel(indexOfCluster{itrj});

  index_from = 1:(nframe-tau);
  index_to   = (1+tau):nframe;
  indexOfCluster_from = indexOfCluster{itrj}(index_from);
  indexOfCluster_to   = indexOfCluster{itrj}(index_to);

  %% ignore invalid indices
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

  %% calc count matrix
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

%% setup
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

%% optimization by L-BFGS-B
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
function [f, g] = myfunc_column(x, c, c_i, nstate);
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
