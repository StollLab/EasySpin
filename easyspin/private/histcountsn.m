function [n,bin] = histcountsn(X,edges)

nPoints = size(X,1);
nDims = size(X,2);

% Discretize data along each axis
inrange = true(nPoints,1);
for d = 1:nDims
  Xd{d} = discretize(X(:,d),edges{d});
  inrange(isnan(Xd{d})) = false;
  siz(d) = numel(edges{d})-1;
end

% Eliminate out-of-range data points
for d = 1:nDims
  Xd{d} = Xd{d}(inrange);
end

% Calculate histogram counts for linearized index
Xlin = sub2ind(siz,Xd{:});
n = histcounts(Xlin,1:prod(siz)+1);

% Reshape final size
n = reshape(n,siz);
bin = cat(2,Xd{:});

end
