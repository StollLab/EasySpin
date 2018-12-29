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

[nPoints, nDims] = size(data);

nIters = 200;

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
sumdist = zeros(1,nRepeats);

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
    randidx = randperm(nPoints,nClusters);
    centroids = data(randidx,:);
  end

%   idxClusters = randi(nClusters,[nPoints,1]);
%   idxClustersLast = zeros(nPoints,1);t

  for iIter = 1:nIters

    centroidsLast = centroids;

    % calculate distances to centroids
    distances = zeros(nPoints,nClusters);
    for k = 1:size(centroids,1)
      % loop over centroids
      % TODO: needs correct distance metric for each angle dimension
      % (modulo 2pi, subtract 2pi if >pi, abs?)
      distances(:,k) = sqrt(sum((data-centroids(k,:)).^2, 2));
    end
    
    % assignment step
    [minDistances, idxClusters] = min(distances,[],2);

    % update step
    for iCluster = 1:nClusters
      centroids(iCluster,:) = mean(data(idxClusters==iCluster,:), 1);
    end
    
    centroidsChange = centroids - centroidsLast;
    
    if iIter > 1
%       if all(abs(centroidsChange(:)) < 1e-4), break, end
      if all(idxClustersLast==idxClusters), break, end
    end
    
    idxClustersLast = idxClusters;

  end
  
  if iIter == nIters
    warning('Centroid changes did not converge within %d iterations.', nIters)
  end
  
  centroids = centroidsLast;
  
  sumdist(1,iRepeat) = sum(minDistances);
  
  if verbosity
    fprintf('Best sum of distances is %0.3f on repeat %d.\n',sumdist(1,iRepeat),iRepeat)
  end
  
  idxClustersAll{iRepeat} = idxClusters;
  centroidsAll{iRepeat} = centroids;
  
end

[~,iBest] = min(sumdist,[],2);
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