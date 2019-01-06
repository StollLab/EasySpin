% cardamom_kmeans  Perform k-means clustering on input data using random
%                  initialization
%
%   [idxBest,centroidsBest] = cardamom_kmeans(data,nClusters,nRepeats)
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

function [idxBest, centroidsBest] = cardamom_kmeans(data,nClusters,nRepeats,centroids0,verbosity)

if ~exist('nRepeats','var'), nRepeats = 1; end
if ~exist('centroids0','var'), centroids0 = []; end
if ~exist('verbosity','var'), verbosity = 0; end

[nPoints, nDims] = size(data);

nIterations = 200;
centroidsChangeThreshold = 1e-4;

if ~isempty(centroids0)
  if size(centroids0,2) ~= nDims || size(centroids0,1) > nClusters
    error('centroids0 must have at most nClusters rows and have exactly nDims columns.')
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

  % Initialize centroids and cluster indices
  if ~isempty(centroids0)
    nClusters0 = size(centroids0,1);
    centroids(1:nClusters,:) = centroids0;
    if nClusters0 < nClusters
      % generate additional clusters using random samples from data
      randidx = randperm(nPoints,nClusters-nClusters0);
      centroids(nClusters:nClusters0,:) = data(randidx,:);
    end
  else
%     randidx = randperm(nPoints,nClusters);
%     centroids = data(randidx,:);
    centroids = kppseed(data, nClusters);
  end

%   idxClusters = randi(nClusters,[nPoints,1]);
%   idxClustersLast = zeros(nPoints,1);t

  for iIter = 1:nIterations

    centroidsLast = centroids;

    % calculate distances to centroids
    distances2 = zeros(nPoints,nClusters);
    for k = 1:nClusters
      if implicitSingletonExpansion
        differences = abs(data-centroids(k,:)); % R2016b and later
      else
        differences = abs(bsxfun(@minus,data,centroids(k,:))); % pre-R2016b
      end
      differences = pi - abs(differences-pi);
      distances2(:,k) = sum(differences.^2, 2);
    end
    
    % assignment step
    [minDistances2, idxClusters] = min(distances2,[],2);

    % update step
    for iCluster = 1:nClusters
%       centroids(iCluster,:) = mean(data(idxClusters==iCluster,:), 1);
      centroids(iCluster,:) = calc_centroids_pbc(data(idxClusters==iCluster,:),2*pi);
    end
    
    centroidsChange = centroids - centroidsLast;
    
    if iIter > 1
      if max(abs(centroidsChange(:))) < centroidsChangeThreshold, break, end
%       if all(idxClustersLast==idxClusters), break, end
    end
    
    idxClustersLast = idxClusters;

  end
  
  if iIter == nIterations
    warning('Centroid changes did not converge within %d iterations.', nIterations)
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

% if ~isempty(centroids0)
%   nRepeats = 1;
%   if nRepeats > 1
%     centroids0 = repmat(centroids0,[1,1,nRepeats]);
%   end
%   opts = statset('Display','final','MaxIter',nIters);
%   [idxBest,centroidsBest] = kmeans(data,size(centroids0,1),...
%                                   'Distance','sqeuclidean',...
%                                   'Replicates',nRepeats,...
%                                   'Start',centroids0,...
%                                   'Options',opts);
% else
%   opts = statset('Display','final','MaxIter',nIters);
%   [idxBest,centroidsBest] = kmeans(data,nClusters,...
%                                   'Distance','sqeuclidean',...
%                                   'Replicates',nRepeats,...
%                                   'Options',opts);
% end

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
  dist = dist_pbc(data-centroids0(k-1,:), 2*pi);
  D = min(D,sum(dist.^2,2));
  D = D./sum(D);
  centroids0(k,:) = data(find(rand<cumsum(D),1),:);
end

end

function dist = dist_pbc(dist,W)

w = W/2;

% dist = w - abs(pi - abs(dist));

idx1 = dist > w;
idx2 = dist < -w;

dist(idx1) = dist(idx1) - W;
dist(idx2) = dist(idx2) + W;

end

function centroids = calc_centroids_pbc(x,W)
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