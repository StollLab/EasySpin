% cardamom_emghmm      Optimize a multivariate Gaussian-hidden Markov model 
%                     parameter set using the expectation-maximization 
%                     algorithm
%
%  [LL, ModelOut] = cardamom_emghmm(data, ModelIn, maxIter)
%
%   Input:
%
%     data           numeric matrix, size = (nPoints,nDims)
%                    input data
%
%     ModelIn        structure 
%                    contains initial model parameters as fields
%
%       prior        numeric matrix, size = (nStates,)
%                    equilibrium probability distribution for starting
%                    states
%
%       transmat     numeric matrix, size = (nStates,nStates,)
%                    transition probability matrix
%
%       mu           numeric matrix, size = (nDims,nStates)
%                    means of the Gaussians representing the emmission
%                    probability distributions
%
%       Sigma        numeric matrix, size = (nDims,nDims,nStates)
%                    covariance matrices of the Gaussians
%
%     maxIter        integer
%                    maximum number of iterations to run the optimization
%                    algorithm
%
%
%   Output:
%
%     logL           double
%                    log-likelihood of the HMM model parameter values
%
%     ModelOut       structure 
%                    contains optimized model parameters as fields

%  Implementation based on:
%    Sezer, Freed, Roux, J. Phys. Chem. B 112, 11014 (2008)
%    Rabiner, Proc. IEEE 77, 257 (1989)
%    HMM Toolbox by Kevin Murphy
%      (http://www.ai.mit.edu/~murphyk/Software/hmm.html)

function [logL, ModelOut] = cardamom_emghmm(data, ModelIn, maxIter)
   
prior = ModelIn.prior;
transmat = ModelIn.transmat;
mu = ModelIn.mu;
Sigma = ModelIn.Sigma;

% check prior
if ~iscolumn(prior)
  error('Prior must be a column vector of size (nStates,1).')
end
[nStates] = numel(prior);

% check transmat
if ~ismatrix(transmat) || size(transmat,1)~=nStates || size(transmat,1)~=size(transmat,2)
  error('Transmat must be a square matrix of size (nStates,nStates).')
end

% check Sigma
if size(Sigma,1)~=size(Sigma,2)
  error('Sigma must be of size (nDims,nDims,nStates).')
end
[nDims,nDims,nStates] = size(Sigma);

% check mu
if size(mu,1)~=nDims || size(mu,2)~=nStates
  error('mu must be of size (nDims,nStates).')
end

% check data
if size(data,2)~=nDims
  error('data must be of size (nPoints,nStates).')
end
nSteps = size(data,1);

% O = 3;     % Number of coefficients in a vector (number of dimensions)
% T = 500;    % Number of vectors in a sequence (length of chains)
% nex = 5;  % Number of sequences (number of chains)
% M = 1;     % Number of mixtures (number of Gaussians per state)
% Q = 2;     % Number of states;

converged = 0;
nIter = 1;

thresh = 1e-4;

% exp_num_trans = 0;
% exp_num_visits1 = 0;
logL_old = inf;
transmat_old = inf*ones(nStates);

while ~converged && nIter<maxIter+1
  
  % expectation step ------------------------------------------------------
  % -----------------------------------------------------------------------
  
  % calculate emmission probabilities
%   b = multivargauss(data, mu, Sigma);  
  b = mixgauss_prob(data.', mu, Sigma).';
  
  [logL, alpha, beta, gamma, sumxi, scale] = forwardbackward(data, transmat, b, prior);
%   [logL, alpha, beta, gamma, xi] = forwardbackward(data, transmat, b, prior);
  
  gamma = gamma.';

  % maximization step -----------------------------------------------------
  % -----------------------------------------------------------------------
  
  % update transmat
  % Eq. 51 from Ref
  prior = gamma(1,:).';
  
  transmat = sumxi;
%   for i = 1:nStates
%     for j = 1:nStates
%       transmat(i,j) = sum(xi(i,j,1:end-1)) / sum(gamma(i,1:end-1));
%     end
%   end
  
  transmat = mk_stochastic(transmat);
%   [transmat, prior] = msmtransitionmatrix(transmat, 1000);
  
  % update mu
  mubar = zeros(nDims,nStates);
  for i = 1:nStates
    mubar(:,i) = sum(gamma(:,i).*(data.' - mu(:,i)).') / sum(gamma(:,i));
  end
  mu = mu + mubar;
  
  % update Sigma
  SigmaBar = zeros(nDims,nDims,nStates);
  for i = 1:nStates
    diffdatamu = (data.' - mu(:,i));
    dmu = zeros(nDims,nDims,nSteps);
    for iStep = 1:nSteps
      dmu(:,:,iStep) = diffdatamu(:,iStep)*diffdatamu(:,iStep).';
    end
    SigmaBar(:,:,i) = sum(permute(gamma(:,i),[2,3,1]).*dmu,3) / sum(gamma(:,i));
  end
    
%   exp_num_trans = exp_num_trans + sumxi; % sum(xi,3);
%   exp_num_visits1 = exp_num_visits1 + gamma(1,:); % will set to prior later
  
%  if adj_prior
%     prior = normalise(exp_num_visits1);
%   end
%   if adj_trans 
%     transmat = mk_stochastic(exp_num_trans);
%   end
%   if adj_mix
%     mixmat = mk_stochastic(postmix);
%   end
%   if adj_mu | adj_Sigma
%     [mu2, Sigma2] = mixgauss_Mstep(postmix, m, op, ip, 'cov_type', cov_type);
%     if adj_mu
%       mu = reshape(mu2, [O Q M]);
%     end
%     if adj_Sigma
%       Sigma = reshape(Sigma2, [O O Q M]);
%     end
%   end

%   emission = emission0;
% 
%   % non-reversible transmat
%   transmat = zeros(nStates, nStates);
%   for iStep = 2:nSteps
%     log_xi = log_alpha(iStep-1, :)' + log_beta(iStep, :);
%     transmat = transmat + exp(log_xi + log_emission0(:, data(iStep))' + log_T0)./factor(iStep);
%   end
%   transmat = transmat./sum(transmat, 2);
%   transmat(isnan(transmat)) = 0;

%   transmat = 
  
  % reversible transmat
%   [transmat, prior] = msmtransitionmatrix(exp_num_trans, 1000);
%   [transmat, prior] = msmtransitionmatrix(transmat, 1000);

  fprintf(1, 'iteration %d, loglik = %f\n', nIter, logL)
%   fprintf(1, '(%0.3f,%0.3f,%0.3f)\n', mu(1,1), mu(2,1), mu(3,1))
%   fprintf(1, '%0.3f,%0.3f,\n%0.3f,%0.3f\n\n', transmat(1,1), transmat(1,2), transmat(2,1), transmat(2,2))

  % check convergence
%   converged = abs(logL_old - logL) < thresh;
  converged = all(abs(transmat_old(:) - transmat(:)) < thresh);
  
  transmat_old = transmat;
  logL_old = logL;
  nIter = nIter + 1;
  
end

ModelOut.prior = prior;
ModelOut.transmat = transmat;
ModelOut.mu = mu;
ModelOut.Sigma = Sigma;

end

function b = multivargauss(x, mu, Sigma)
% x: points on the grid to be evaluated, size = (nPoints,nDims)
% mu: means of the distributions, size = (nDims,nStates)
% Sigma: covariance matrix, size = (nDims,nDims,nStates)
% b: emmission probabilities, size = (nPoints,nStates)

if size(Sigma,1)~=size(Sigma,2)
  error('Sigma must be of size (nDims,nDims,nStates).')
end

[nDims,nDims,nStates] = size(Sigma);

if size(mu,1)~=nDims || size(mu,2)~=nStates
  error('mu must be of size (nDims,nStates).')
end

if size(x,2)~=nDims
  error('x must be of size (nPoints,nStates).')
end

nPoints = size(x,1);

b = zeros(nPoints,nStates);

for iState = 1:nStates
  Sigma0 = Sigma(:,:,iState);
  mu0 = mu(:,iState).';
  p = 1/(2*pi*det(Sigma0)^(1/2))*exp(-0.5.*(x-mu0)*inv(Sigma0)*(x-mu0).');
  p = diag(p);
  b(:,iState) = p;
end

end

function [logL,alpha,beta,gamma,sumxi,scale] = forwardbackward(data,transmat,b,prior)
% function [logL,alpha,beta,gamma,xi] = forwardbackward(data,transmat,b,prior)

[nSteps, nDims] = size(data);
nStates = size(transmat,1);

% b = b.';

% forward algorithm -------------------------------------------------------

% initialize forward parameter alpha
alpha = zeros(nStates, nSteps);
% alpha = zeros(nSteps, nStates);
% factor = zeros(nSteps, 1);

% prior = prior.';
b = b.';

[alpha(:, 1), scale(1)] = normalize(prior.*b(:,1));
% alpha(1, :) = normalize((prior.').*b(1,:)));

% factor(1) = sum(alpha(1, :));
% alpha(1, :) = alpha(1, :)./factor(1);

sumxi = zeros(nStates,nStates);
for iStep = 2:nSteps
  sumalphatrans = transmat.' * alpha(:,iStep-1);
  [alpha(:, iStep),scale(iStep)] = normalize(sumalphatrans(:).*b(:,iStep));
%   sumalphatrans = sum(alpha(iStep-1,:) .* transmat);
%   [alpha(iStep, :),scale(iStep)] = normalize((sumalphatrans.*b(iStep,:)).');

  sumxi = sumxi + normalize((alpha(:,iStep-1) * b(:,iStep)') .* transmat);
%   sumxi = sumxi + normalise((alpha(:,iStep-1)*b(iStep,:).') .* transmat);
%   factor(iStep) = sum(alpha(iStep, :));
%   alpha(iStep, :) = alpha(iStep, :)./factor(iStep);
end

% logL = sum(log(factor));

% backward algorithm ------------------------------------------------------

% initialize beta
beta = zeros(nStates, nSteps);
beta(:, end) = 1;
% beta = zeros(nSteps, nStates);
% beta(end, :) = 1;

for iStep = (nSteps-1):-1:1
  [beta(:, iStep),trash] = normalize( transmat * (b(:,iStep+1) .* beta(:,iStep+1)) );
%   [beta(iStep, :),trash] = normalize(sum(transmat .* b(iStep+1,:).' ...
%                            .* beta(iStep+1,:).'));
end

condProb = sum(alpha.*beta,1);

% compute gamma and xi
gamma = normalize(alpha.*beta);
% gamma = alpha.*beta./condProb;

% xi = zeros(nStates,nStates,nSteps);
% for iStep = 1:nSteps-1
%   for i = 1:nStates
%     for j = 1:nStates
%       xi(i,j,iStep) = ( (alpha(i,iStep)*transmat(i,j)) ...
%                        *(b(j,iStep+1).*beta(j,iStep+1)).') ./ condProb(iStep); 
% %       xi(i,j,iStep) = ( (alpha(iStep,i)*transmat(i,j)) ...
% %                        *(b(iStep+1,j).*beta(iStep+1,j)).') ./ condProb(iStep); 
%     end
%   end
% end

% gamma = zeros(nSteps,nStates);
% for iStep = 1:nSteps
%   gamma(iStep,:) = normalize(alpha(iStep,:).*beta(iStep,:));
% end

% compute xi
% xi = zeros(nStates,nStates,nSteps);

% logL = log(sum(alpha(end,:)));
logL = sum(log(scale));

% log_alpha = log(alpha);
% log_beta = log(beta);
% log_T0 = log(transmat);
% log_b = log(b);

end

function [B,scale] = normalize(A)
% rescale an array so that all entries sum to one

scale = sum(A(:));
if abs(scale) <1e-14
  error('Zero encountered during normalization.')
end

B = A/scale;

end