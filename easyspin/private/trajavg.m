% performs a moving average over nWindow time steps of a single, long 
% trajectory of a scalar or vector quantity

function vAvg = trajavg(v,MD,Par)

% assumes that size(v,1) is the time dimension and that size(v,2) is the
% number of components

if ~isfield(MD,'dt')
  error('MD.dt, the time step for the trajectory, is missing.')
end

if ~isfield(Par,'dt')
  error('Par.dt, the time step for the trajectory, is missing.')
end

% size of averaging window
nWindow = ceil(Par.dt/MD.dt);

% size of trajectory after averaging
M = floor(MD.nSteps/nWindow);

nComps = size(v,2);

vAvg = zeros(M,nComps);

for iComp=1:nComps
  temp = conv(v(:,iComp), 1/nWindow*ones(nWindow,1), 'valid');
  vAvg(:,iComp) = temp(1:nWindow:end);
end

end

