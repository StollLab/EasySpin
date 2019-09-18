% given a time step dt and an averaging window tWindow, performs a moving 
% average of a trajectory of a scalar or vector quantity

function vAvg = trajavg(v,dt,tWindow)

if tWindow<dt
  error('The averaging window tWindow must be greater than or equal to the time step dt.')
end

% assumes that size(v,1) is the time dimension and that size(v,2) is the
% number of components

% size of averaging window
nWindow = ceil(tWindow/dt);

nSteps = size(v,1);

% size of trajectory after averaging
M = floor(nSteps/nWindow);

nComps = size(v,2);

vAvg = zeros(M,nComps);

for iComp=1:nComps
  temp = conv(v(:,iComp), 1/nWindow*ones(nWindow,1), 'valid');
  vAvg(:,iComp) = temp(1:nWindow:end);
end

end

