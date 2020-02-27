function ok = test(opt)

% Test k-means clustering on 1D angular data

rng(1);

% Generate clustered angular 1D data
% (in degrees)
chicluster{1} = [0:5:20 -20:5:5];
chicluster{2} = 70:5:110;
chicluster{3} = [160:5:180 -175:5:-160];
chicluster{4} = -110:5:-70;
% (convert to radians)
chicluster = cellfun(@(x)x*pi/180,chicluster,'UniformOutput',false);

% Combine clustered data into a single list
chi = [];
clusteridx = [];
nClusters = numel(chicluster);
for k = 1:nClusters
  n = numel(chicluster{k});
  chi = [chi; chicluster{k}(:)];
  clusteridx = [clusteridx; ones(n,1)*k];
end

% Scrample order of datapoints
idx = randperm(numel(chi));
chi = chi(idx);
clusteridx = clusteridx(idx);

% Run k-means clustering
seeds = [-60; 30; 100; 160]*pi/180;
nRepeats = 1;
[idx_km,~] = runprivate('mdhmm_kmeans',chi,nClusters,nRepeats,seeds);

% Check whether k-means produced correct clustering
for k = 1:nClusters
  idx_ = idx_km(clusteridx == k);
  ok(k) = all(idx_==idx_(1));
end

% Plot if wanted
if opt.Display
  clf
  plot(cos(chi),sin(chi),'.')
  hold on
  axis equal
  axis([-1 1 -1 1]);
  for k = 1:4
    chi_ = chi(idx_km==k);
    plot(cos(chi_),sin(chi_),'o');
  end
end
