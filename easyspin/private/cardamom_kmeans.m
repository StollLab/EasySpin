% cardamom_kmeans  Perform k-means clustering on input data using random
%                  initialization
%
%   [idxStates,centroids] = cardamom_kmeans(data,nClusters,nRepeats)
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
%   Output:
%
%     idxClusters    numeric matrix, size = (nPoints,1)
%                    cluster assignments using integers in range
%                    (1,nClusters)
%
%     centroids      numeric matrix, size = (nClusters,nDims)
%                    coordinates of centroids

function [idxClusters, centroids] = cardamom_kmeans(data,nClusters,nRepeats,verbosity)

[nPoints, nDims] = size(data);

% % center and scale data
% data = data - mean(data,1);
% data = data./max(data,[],1);

idxClusters = {1,nRepeats};
centroids = {1,nRepeats};
sumdist = zeros(1,nRepeats);

for iRepeat = 1:nRepeats

  % initialize centroids and cluster indices
  randidx = randperm(nPoints,nClusters);
  centroidsNew = data(randidx,:);

  idxClustersNew = randi(nClusters,[nPoints,1]);
  idxClustersLast = zeros(nPoints,1);

  while any(idxClustersNew ~= idxClustersLast)

    idxClustersLast = idxClustersNew;

    % assignment step
    [dist, idxClustersNew] = calcdistances(data, centroidsNew, nPoints, nClusters);

    % update step
    for iCluster = nClusters
      centroidsNew(iCluster,:) = mean(data(idxClustersNew==iCluster,:),1);
    end

  end
  
  sumdist(1,iRepeat) = sum(dist);
  
  if verbosity
    fprintf('Best sum of distances is %0.3f on repeat %d.\n',sumdist(1,iRepeat),iRepeat)
  end
  
  idxClusters{iRepeat} = idxClustersNew;
  centroids{iRepeat} = centroidsNew;
  
end

[~,iBest] = min(sumdist,[],2);
idxClusters = idxClusters{iBest};
centroids = centroids{iBest};

end

function [distances,idx] = calcdistances(data, centroids, nPoints, nClusters)

distances = zeros(nPoints,nClusters);
for k = 1:size(centroids,1)
  % loop over centroids
  distances(:,k) = sqrt(sum((data - centroids(k,:)).^2,2));
end

[distances,idx] = min(distances,[],2);

end