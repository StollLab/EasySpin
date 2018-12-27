function [err,data] = test(opt,olddata)
% Regression test for k-means clustering function

rng(1)

% Simulate clustered data
% -------------------------------------------------------------------------

% load('.\mdfiles\MTSSL_polyAla_traj.mat')
% data = squeeze(Traj.dihedrals(1:2,:,:)).';

nPoints = 100;
nDims = 2;

testData1 = randn(nPoints,nDims);
testData2 = randn(nPoints,nDims)+3;
testData3 = randn(nPoints,nDims)-3;

testData = [testData1; testData2; testData3];

nClusters = 3;
nRepeats = 5;
verbosity = 0;

centroids0 = [];

% Perform clustering
% -------------------------------------------------------------------------

[idxBest, centroidsBest] = runprivate('cardamom_kmeans',...
                                      testData,nClusters,nRepeats,centroids0,verbosity);


data.centroidsBest = centroidsBest;

if ~isempty(olddata)
  err = any(abs(olddata.centroidsBest-centroidsBest)>1e-10);
else
  err = [];
end

end
