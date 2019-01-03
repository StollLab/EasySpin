function [err,data] = test(opt,olddata)
% Regression test for expectation maximization for GHMM function

rng(1)

% Simulate trajectory from HMM
% -------------------------------------------------------------------------

% load('.\mdfiles\MTSSL_polyAla_traj.mat')
% testData = squeeze(Traj.dihedrals(1:2,:,:)).';

nDims = 2;
nStates = 3;
nRepeats = 5;
verbosity = 0;

Sys.TransProb = [  0.8, 0.1, 0.1; 
                   0.1, 0.8, 0.1;
                   0.1, 0.1, 0.8 ];
% Sys.TransProb = [  0.8, 0.2; 
%                    0.2, 0.8 ];
Par.nTraj = 1;
Par.dt = 1e-9;
Par.nSteps = 1000;
Opt.statesOnly = true;
[~,stateTraj] = stochtraj_jump(Sys,Par,Opt);

testData = zeros(Par.nSteps, nDims);

width = 40;

for iStep = 1:Par.nSteps
  state = stateTraj(iStep);
  switch state
    case 1
      testData(iStep,:) = [width*randn() - 150, width*randn()];
    case 2
      testData(iStep,:) = width*randn(1,2);
    case 3
      testData(iStep,:) = [width*randn() + 150, width*randn() + 150];
  end
end

testData = testData/180*pi;

testDataWrapped = testData;
testDataWrapped(testData>pi) = testData(testData>pi) - 2*pi;
testDataWrapped(testData<-pi) = testData(testData<-pi) + 2*pi;

centroids0 = [];

% Perform clustering
% -------------------------------------------------------------------------

[idxBest, centroidsBest] = runprivate('cardamom_kmeans',...
                                      testDataWrapped,nStates,nRepeats,centroids0,verbosity);

                                    
% Train HMM parameters
% -------------------------------------------------------------------------

% initialize the means and covariance matrices for the HMM
mu0 = centroidsBest.';
Sigma0 = zeros(nDims,nDims,nStates);
for iState = 1:nStates
  idxState = idxBest==iState;
  Sigma0(:,:,iState) = cov_pbc(testData(idxState,:), mu0(:,iState).', 2*pi);
%   Sigma0(:,:,iState) = cov(testData(idxState,:));
end

mixmat0 = ones(nStates,1);

transmat0 = rand(nStates);
transmat0 = transmat0./sum(transmat0,2);

eqdistr0 = rand(nStates,1);
eqdistr0 = eqdistr0./sum(eqdistr0,1);

% EM function wants an array of shape (nDims,nSteps,nTraj)
testDataWrapped = permute(testDataWrapped,[2,1,3]);

[eqdistr1,transmat1,mu1,Sigma1,~] = ...
    runprivate('cardamom_emghmm',testDataWrapped,eqdistr0,transmat0,mu0,Sigma0,[],verbosity);

% Compare
% -------------------------------------------------------------------------

data.transmat = transmat1;
data.eqdistr = eqdistr1;

% hold on
% % scatter(testDataWrapped(1,:), testDataWrapped(2,:), 10, 'filled')
% for iState = 1:nStates
%   idxState = idxBest==iState;
%   scatter(testDataWrapped(1,idxState), testDataWrapped(2,idxState), 10, 'filled')
% end
% scatter(mu0(1,:), mu0(2,:), 40, 'filled')
% scatter(mu1(1,:), mu1(2,:), 40, 'green')
% hold off


if ~isempty(olddata)
  err = any(any(abs(olddata.transmat-transmat1)>1e-10)) ...
       || any(abs(olddata.eqdistr-eqdistr1)>1e-10);
else
  err = [];
end

end



% Helper functions
% -------------------------------------------------------------------------

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
