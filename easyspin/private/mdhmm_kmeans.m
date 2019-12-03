% mdhmm_kmeans  Perform k-means clustering on input data using random
%                  initialization
%
%   [idxBest,centroidsBest] = mdhmm_kmeans(data,nClusters,nRepeats)
%
%
%     data           numeric matrix, size = (nPoints,nDims)
%                    input data
%
%     nClusters      integer
%                    number of clusters to assign
%
%     nRepeats       integer
%                    number of times to repeat the clustering, where the
%                    final result is chosen as the one the gave the least
%                    sum of distances
%
%     centroids0     numeric matrix, size = (nClusters,nDims)
%                    (optional) initial cluster centroids
%
%     verbosity      0 or 1
%                    if set to 1, display information regarding the
%                    progress
%    
%   Output:
%
%     idxBest        numeric matrix, size = (nPoints,1)
%                    cluster assignments using integers in range
%                    (1,nClusters)
%
%     centroidsBest  numeric matrix, size = (nClusters,nDims)
%                    coordinates of centroids

function [idxBest, centroidsBest] = mdhmm_kmeans(data,nClusters,nRepeats,centroids0,verbosity)

if ~exist('nRepeats','var'), nRepeats = 1; end
if ~exist('centroids0','var'), centroids0 = []; end
if ~exist('verbosity','var'), verbosity = 0; end

maxIterations = 500;
centroidsChangeThreshold = 1e-4; % radians

[nPoints, nDims] = size(data);

if ~isempty(centroids0)
  if size(centroids0,2) ~= nDims
    error('centroids0 must have exactly %d columns.',nDims)
  end
  if size(centroids0,1) > nClusters
    error('centroids0 must have at most %d rows.',nClusters)
  end
  nRepeats = 1;
end

% center and scale data
% data = data - mean(data,1);
% data = data./max(data,[],1);

idxClustersAll = {1,nRepeats};
centroidsAll = {1,nRepeats};
sumdist2 = zeros(1,nRepeats);

% Use implicit singleton expansion if available, i.e. in R2016b (9.1) or later
implicitSingletonExpansion = ~verLessThan('matlab','9.1');

for iRepeat = 1:nRepeats
  
  centroids = zeros(nClusters,nDims);

  % Initialize centroids and cluster indices
  if ~isempty(centroids0)
    nClusters0 = size(centroids0,1);
    centroids(1:nClusters0,:) = centroids0;
    if nClusters0 < nClusters
      % generate additional clusters using random samples from data
      randidx = randperm(nPoints,nClusters-nClusters0);
      centroids(nClusters0+1:nClusters,:) = data(randidx,:);
    end
  else
%     randidx = randperm(nPoints,nClusters);
%     centroids = data(randidx,:);
    centroids = kppseed(data, nClusters);
  end

%   idxClusters = randi(nClusters,[nPoints,1]);
%   idxClustersLast = zeros(nPoints,1);t

  for iIter = 1:maxIterations

    centroidsLast = centroids;

    if iIter==1
      % check for centroids that have no data points assigned to them
      allClustersAssigned = false;
      clustersAssigned = false(1,nClusters);
      
      while ~allClustersAssigned
        % calculate distances to centroids
        distances2 = calcdistances(data, centroids, implicitSingletonExpansion);
        % assignment step
        [minDistances2, idxClusters] = min(distances2,[],2);
        clustersAssigned(unique(idxClusters)) = true;
        allClustersAssigned = all(clustersAssigned);
        
        if ~allClustersAssigned
          % move the centroids that have no assigned data points
          newCentroids = pi*(2*rand(sum(~clustersAssigned), nDims)-1);
          centroids(~clustersAssigned, :) = newCentroids;
        end
        
      end
      
    else
      % calculate distances to centroids
      distances2 = calcdistances(data, centroids, implicitSingletonExpansion);
      % assignment step
      [minDistances2, idxClusters] = min(distances2,[],2);
      
    end

    % update step
    for iCluster = 1:nClusters
      centroids(iCluster,:) = calc_centroids_pbc(data(idxClusters==iCluster,:));
    end
    
    centroidsChange = dist_pbc(centroids,centroidsLast);
    %disp([iRepeat iIter max(abs(centroidsChange(:)))]);
    
    if iIter > 1
      if max(abs(centroidsChange(:))) < centroidsChangeThreshold, break, end
%       if all(idxClustersLast==idxClusters), break, end
    end
    
%     idxClustersLast = idxClusters;

  end
  
  if iIter == maxIterations
    warning('Centroid changes did not converge within %d iterations.', maxIterations)
  end
  
  centroids = centroidsLast;
  
  sumdist2(1,iRepeat) = sum(minDistances2);
  
  if verbosity
    fprintf('    repeat %d: best sum of squared distances is %0.3f after %d iterations\n',iRepeat,sumdist2(1,iRepeat),iIter);
  end
  
  idxClustersAll{iRepeat} = idxClusters;
  centroidsAll{iRepeat} = centroids;
  
end

[~,iBest] = min(sumdist2,[],2);
idxBest = idxClustersAll{iBest};
centroidsBest = centroidsAll{iBest};

end

% Helper functions
% -------------------------------------------------------------------------

function centroids0 = kppseed(data, nSeeds)
% Seed k-means clustering using k-means++ algorithm, where each successive
% seed is chosen randomly with probability proportional to the distance to
% the nearest center

% Adapted from code by Mo Chen (sth4nth@gmail.com).
[nPoints,nDims] = size(data);
D = inf(nPoints,1);

centroids0 = zeros(nSeeds,nDims);
centroids0(1,:) = data(ceil(nPoints*rand),:);
for k = 2:nSeeds
  dist = dist_pbc(data,centroids0(k-1,:));
  D = min(D,sum(dist.^2,2));
  D = D./sum(D);
  centroids0(k,:) = data(find(rand<cumsum(D),1),:);
end

end

function distances2 = calcdistances(data, centroids, implicitSingletonExpansion)

nPoints = size(data, 1);
nClusters = size(centroids, 1);

distances2 = zeros(nPoints,nClusters);
for k = 1:nClusters
  if implicitSingletonExpansion
    differences = pi - abs(pi - abs(data-centroids(k,:))); % R2016b and later
  else
    differences = pi - abs(pi - abs(bsxfun(@minus,data,centroids(k,:)))); % pre-R2016b
  end
  distances2(:,k) = sum(differences.^2, 2);

end

end

function dist = dist_pbc(x1,x2)
% Distance between two angles x1 and x2, taking circular nature into account
dist = pi - abs(pi - abs(x1-x2));
end

function centroids = calc_centroids_pbc(x)
% Calculate centroid coordinates while accounting for periodic boundary
% conditions

refPoint = x(1,:);

% Unwrap points with respect to the reference point
distances = x - refPoint;
distances(distances>pi) = distances(distances>pi) - 2*pi;
distances(distances<-pi) = distances(distances<-pi) + 2*pi;
xUnwrapped = refPoint + distances;

% Wrap centroids as needed
centroids = mean(xUnwrapped,1);
centroids(centroids>pi) = mean(centroids(centroids>pi) - 2*pi,1);
centroids(centroids<-pi) = mean(centroids(centroids<-pi) + 2*pi,1);

end