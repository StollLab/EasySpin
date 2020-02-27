function ok = test(opt)

% Direct test on cluster indices for k-means clustering function in 2D for 
% Gaussian clusters and periodic boundary conditions

rng(1)

% Simulate clustered data
% -------------------------------------------------------------------------

nPoints = 100;
nDims = 2;

% Width of normal distributions
width = 20;

% Centers of normal distributions
centroidsTrue = [   0,   0; 
              150, 150; 
             -150,   0 ]/180*pi;

% Simulate clusters of data
testData = [];
clusteridx = [];
for k = 1:size(centroidsTrue,1)
  randData = width*randn(nPoints,nDims);
  testData = cat(1, testData, randData/180*pi + centroidsTrue(k,:));
  clusteridx = cat(1, clusteridx, ones(nPoints,1)*k);
end

% Wrap data to periodic boundary conditions
testDataWrapped = testData;
testDataWrapped(testData>pi) = testDataWrapped(testData>pi) - 2*pi;
testDataWrapped(testData<-pi) = testDataWrapped(testData<-pi) + 2*pi;

% Scrample order of datapoints
idx = randperm(size(testDataWrapped,1));
testDataWrapped = testDataWrapped(idx,:);
clusteridx = clusteridx(idx);

% Perform clustering with random seeds
nClusters = 3;
nRepeats = 10;
verbosity = 0;
centroids0 = [];

[idxBest, centroidsBest] = runprivate('mdhmm_kmeans',...
                                      testDataWrapped,nClusters,nRepeats,centroids0,verbosity);

% Check whether k-means produced correct clustering
for k = 1:nClusters
  idx_ = idxBest(clusteridx == k);
  ok(k) = all(idx_==idx_(1));
end

% Plot if wanted
if opt.Display
  clf
  hold on
  for k = 1:nClusters
    scatter(testDataWrapped(idxBest==k,1), testDataWrapped(idxBest==k,2), 10, 'filled')
  end
  scatter(centroidsTrue(:,1), centroidsTrue(:,2), 40, 'filled')
  scatter(centroidsBest(:,1), centroidsBest(:,2), 70)
  hold off
  axis equal
  axis([-pi pi -pi pi]);
end

end
