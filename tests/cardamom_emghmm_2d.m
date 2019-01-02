function [err,data] = test(opt,olddata)
% Regression test for expectation maximization for GHMM function

rng(1)

% Simulate trajectory from HMM
% -------------------------------------------------------------------------

% load('.\mdfiles\MTSSL_polyAla_traj.mat')
% testData = squeeze(Traj.dihedrals(1:2,:,:)).';

nDims = 2;

Sys.TransProb = [  0.8, 0.1, 0.1; 
                   0.1, 0.8, 0.1;
                   0.1, 0.1, 0.8 ];
Par.nTraj = 1;
Par.dt = 1e-9;
Par.nSteps = 1000;
Opt.statesOnly = true;
[~,stateTraj] = stochtraj_jump(Sys,Par,Opt);

testData = zeros(Par.nSteps, nDims);

for iStep = 1:Par.nSteps
  state = stateTraj(iStep);
  switch state
    case 1
      testData(iStep,:) = randn(1,2) + 2;
    case 2
      testData(iStep,:) = randn(1,2) - 2;
    case 3
      testData(iStep,:) = randn(1,2);
  end
end

nStates = 3;
nRepeats = 5;
verbosity = 0;

centroids0 = [];

% Perform clustering
% -------------------------------------------------------------------------

[idxBest, centroidsBest] = runprivate('cardamom_kmeans',...
                                      testData,nStates,nRepeats,centroids0,verbosity);

                                    
% Train HMM parameters
% -------------------------------------------------------------------------

% initialize the means and covariance matrices for the HMM
mu0 = centroidsBest.';
Sigma0 = zeros(nDims,nDims,nStates);
for iState = 1:nStates
  idxState = idxBest==iState;
  Sigma0(:,:,iState) = cov(testData(idxState,:));
end

mixmat0 = ones(nStates,1);

transmat0 = rand(nStates);
transmat0 = transmat0./sum(transmat0,2);

prior0 = rand(nStates,1);
prior0 = prior0./sum(prior0,1);

% EM function wants an array of shape (nDims,nSteps,nTraj)
testData = permute(testData,[2,1,3]);

[prior1,transmat1,mu1,Sigma1,~] = ...
    runprivate('cardamom_emghmm',testData,prior0,transmat0,mu0,Sigma0,[],verbosity);

% Compare
% -------------------------------------------------------------------------

data.transmat = transmat1;

if ~isempty(olddata)
  err = any(abs(olddata.transmat-transmat1)>1e-10);
else
  err = [];
end

end
