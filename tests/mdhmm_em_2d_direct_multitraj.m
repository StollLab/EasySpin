function [err,data] = test(opt,olddata)
% Direct test for multivariate Gaussian HMM expectation maximization 
% function using multiple trajectories

rng(1)

% Simulate state trajectory using dynamic parameters
% -------------------------------------------------------------------------

nDims = 2;
nStates = 3;
nRepeats = 5;
verbosity = 0;

Sys.TransProb = [  0.7, 0.1, 0.2; 
                   0.1, 0.8, 0.1;
                   0.1, 0.3, 0.6 ];
Par.nTraj = 3;
Par.dt = 1e-9;
Par.nSteps = 10000;
Opt.statesOnly = true;
[~,stateTraj] = stochtraj_jump(Sys,Par,Opt);

testData = zeros(Par.nSteps, nDims, Par.nTraj);

% Create observation trajectory using state trajectory and normal 
% distributions assigned to each state
% -------------------------------------------------------------------------

% Width of normal distributions 
width = 20/180*pi;

% Centers of normal distributions
muTrue = [   0,   0; 
              150, 150; 
             -150,   0 ]/180*pi;

for iTraj = 1:Par.nTraj
  traj = stateTraj(iTraj,:);
  for iStep = 1:Par.nSteps
    state = traj(iStep);
    testData(iStep,:,iTraj) = width*randn(1,2) + muTrue(state,:);
  end
end

% Wrap data according to periodic boundary conditions
testDataWrapped = testData;
testDataWrapped(testData>pi) = testData(testData>pi) - 2*pi;
testDataWrapped(testData<-pi) = testData(testData<-pi) + 2*pi;

% Perform clustering using random seeds
% -------------------------------------------------------------------------

nTraj = size(testDataWrapped,3);

% if more than one trajectory, collapse 3rd dim (traj) onto 1st dim (time)
if nTraj > 1
  testDataWrappedCollapsed = [];
  for iTraj = 1:nTraj
    testDataWrappedCollapsed = cat(1, testDataWrappedCollapsed, testDataWrapped(:,:,iTraj));
  end
end

centroids0 = [];
[idxBest, centroidsBest] = runprivate('mdhmm_kmeans',...
                                      testDataWrappedCollapsed,nStates,nRepeats,centroids0,verbosity);

                                    
% Train HMM parameters against observation trajectory
% -------------------------------------------------------------------------

% Initialize the means and covariance matrices for the HMM
mu0 = centroidsBest.';
Sigma0 = zeros(nDims,nDims,nStates);
for iState = 1:nStates
  idxState = idxBest==iState;
%   stateData = testDataWrappedCollapsed(idxState,:);
  Sigma0(:,:,iState) = cov_pbc(testDataWrappedCollapsed(idxState,:), mu0(:,iState).', 2*pi);
end

% Initialize TPM and equilibrium distribution using uniform distributions
TransProb0 = rand(nStates);
TransProb0 = TransProb0./sum(TransProb0,2);

eqDistr0 = rand(nStates,1);
eqDistr0 = eqDistr0./sum(eqDistr0,1);

% EM function wants an array of shape (nDims,nSteps,nTraj)
testDataWrapped = permute(testDataWrapped,[2,1,3]);

% Run expectation maximization algorithm
[~,eqDistr1,TransProb1,mu1,Sigma1] = ...
    runprivate('mdhmm_em',testDataWrapped,eqDistr0,TransProb0,mu0,Sigma0,verbosity);

% Compare
% -------------------------------------------------------------------------

err = 0;

% Assignment of output state indices might not correspond to input state
% indices, so check which ordering gives the best agreement and re-assign
% indices accordingly
idxPerm = perms([1,2,3]);
bestAgreement = inf;
for k = 1:size(idxPerm,1)
  idx = idxPerm(k,:);
  thisAgreement = max(max(abs(mu1(:,idx)-muTrue')));
  if thisAgreement < bestAgreement
    bestAgreement = thisAgreement;
    bestidx = idx;
  end
end

% Reorder HMM parameters based on best agreement
mu1 = mu1(:,bestidx);
temp = TransProb1(:,bestidx);
TransProb1 = temp(bestidx,:);

muDiff = abs(mu1 - muTrue');
transmatChange = abs(Sys.TransProb - TransProb1)./TransProb1;

muThreshold = 5/180*pi;
transmatThreshold = 0.5;

if    (max(muDiff(:)) > muThreshold) ...
   || (max(transmatChange(:)) > transmatThreshold)
  err = 1;
end

data = [];

if opt.Display
  
  hold on
  % scatter(testDataWrapped(1,:), testDataWrapped(2,:), 10, 'filled')
  for iState = 1:nStates
    idxState = idxBest==iState;
    scatter(testDataWrapped(1,idxState), testDataWrapped(2,idxState), 10, 'filled')
  end
%   scatter(mu0(1,:), mu0(2,:), 40, 'filled')
  scatter(muTrue(:,1), muTrue(:,2), 40, 'filled')
  scatter(mu1(1,:), mu1(2,:), 40, 'green')
  hold off
  
end

end

% Helper functions
% -------------------------------------------------------------------------

function c = cov_pbc(x, mu, W)

nPoints = size(x,1);
    
% remove the centers
xc = dist_pbc(x-mu, W);

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

end
