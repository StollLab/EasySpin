function [err,data] = test(opt,olddata)
% Regression test for k-means clustering function

rng(1)

% Simulate clustered data
% -------------------------------------------------------------------------

% load('.\mdfiles\MTSSL_polyAla_traj.mat')
% data = squeeze(Traj.dihedrals(1:2,:,:)).';

nPoints = 100;
nDims = 2;

width = 30;

% testData1 = randn(nPoints,nDims);
% testData2 = randn(nPoints,nDims)+3;
% testData3 = randn(nPoints,nDims)-3;
testData1 = width*randn(nPoints,nDims);
testData2 = width*randn(nPoints,nDims)+150;
testData3 = [20*randn(nPoints,1)-150, 20*randn(nPoints,1)];

% testData = [testData1; testData2];
testData = [testData1; testData2; testData3];
testData = testData/180*pi;
testData(testData>pi) = testData(testData>pi) - 2*pi;
testData(testData<-pi) = testData(testData<-pi) + 2*pi;
% testData = mod(testData,pi);


nClusters = 3;
nRepeats = 5;
verbosity = 0;

centroids0 = [];

% Perform clustering
% -------------------------------------------------------------------------

[idxBest, centroidsBest] = runprivate('mdhmm_kmeans',...
                                      testData,nClusters,nRepeats,centroids0,verbosity);


data.centroidsBest = centroidsBest;

if ~isempty(olddata)
  err = any(abs(olddata.centroidsBest-centroidsBest)>1e-10);
else
  err = [];
end

end
