% esintpol  Interpolates data of a given symmetry. 
%
%   yi = esintpol(y,gridParams,Factor,phi,the,InterpType)
%   
%   y:            Data, in row vectors. Matrices get interpolated along rows.
%   gridParams:   [nKnots, closedPhi, nOctants, maxPhi], as provided by gridparam()
%   Factor:       interpolation factor (only used for Dinfh)
%   InterpType:   interpolation type, 'G3' (default),'L3','L1'
%   phi,theta:    interpolation points

function yi = esintpol(y,gridParams,Factor,phi,the,InterpType)

if nargin<5 || nargin>6
  error('5 or 6 input arguments ar required.')
end

% Set default or user-defined parameters
if nargin<6 || isempty(InterpType)
  InterpType = 'G3';
end

% parse parameters
globalInterpolation = upper(InterpType(1)) == 'G';
cubicInterpolation = InterpType(2) == '3';

% List of unique phi intervals as set by sphgrid(),
% octant numbers and end conditions (p is periodic,
% z is zero endslopes).

% Ci      4p   [0,2*pi)
% C2h     2p   [0,pi)
% S6      2p   [0,2*pi/3)
% C4h     1p   [0,pi/2)
% C6h     1p   [0,pi/3)

% D2h     1z   [0,pi/2]
% Th      1z   [0,pi/2]
% D3d     1z   [0,pi/3]
% D4h     1z   [0,pi/4]
% Oh      1z   [0,pi/4]
% D6h     1z   [0,pi/6]
% Dinfh   0     0
% O3     -1     0

% Convert symmetry to octant number and boundary condition
% flag, error if invalid point group is encountered.

if ischar(gridParams)
  error('Don''t supply symmetry string!');
else
  % maxPhi = Symmetry(1); %not needed
  nKnots = gridParams(1);
  periodic = ~gridParams(2);
  nOctants = gridParams(3);
  maxPhi = gridParams(4);
end
 
if nOctants>0 && nargin<3
  error('Not enough input parameters!');
end

% Compute number of slices
if nOctants>0
  if nOctants==8
    nExpectedData = nKnots*(nKnots-1)/2*4 +1; % upper hemisphere + equator
    nExpectedData = nExpectedData + (nKnots-1)*(nKnots-2)/2*4 + 1; % lower hemisphere
  else
    nExpectedData = nKnots*(nKnots-1)/2*nOctants + 1;
  end
  if ~periodic
    nExpectedData = nExpectedData + nKnots - 1;
  end
else % Dinfh
  nExpectedData = nKnots;
end

% If number of slices is not an integer, throw error.
if nExpectedData~=length(y)
  error('Wrong number of points. Not a valid data set!');
end

switch nOctants
case 0 % Dinfh
  %============================================================
  % Dinfh
  %============================================================
  % Linear interpolation, vectorized
  % Cubic interpolation, local and global, vectorized
  % Slope estimators: simple and Fritsch-Carlson
  %============================================================
  [m,n] = size(y);
  
  if ~cubicInterpolation % linear interpolation
    yi = interp1(y.',1:1/Factor:n);
    
  elseif globalInterpolation % global cubic interpolation
    yi = esspline1d(y,1,1:1/Factor:n);
    
  else % local cubic interpolation
    % prepare [x^3 x^2 x 1] for each point
    x = (0:1/Factor:1-1/Factor).';
    x = [x.^3 x.^2 x ones(Factor,1)];
    % coefficient matrix for Hermite interpolation
    H = [2 -2 1 1; -3 3 -2 -1; 0 0 1 0; 1 0 0 0];
    
    Estimator = 'FC';
    % run over all row vectors in y
    for r = 1:m
      switch Estimator
        case 'avg' % simple slope average
          Tangents = [0 (y(r,3:end)-y(r,1:end-2))/2 0];
        case 'FC' % Fritsch-Carlson monotone slopes
          del = diff(y(r,:));
          k = find(sign(del(1:n-2)).*sign(del(2:n-1))>0);
          dmax = max(abs(del(k)), abs(del(k+1)));
          dmin = min(abs(del(k)), abs(del(k+1)));
          Tangents = zeros(1,n);
          Tangents(k+1) = 2*dmin.*dmax./(del(k)+del(k+1));
      end
      
      % control vectors for all intervals
      C = [y(r,1:end-1); y(r,2:end); Tangents(1:end-1); Tangents(2:end)];
      
      % evaluate interpolating Hermite splines
      yii = x*H*C;
      yi(r,:) = [yii(:).' y(r,n)];
    end
  end
  %============================================================
  %============================================================
  %============================================================
  
otherwise
  
  % Convert triangular to rectangular grid.
  z = rectify(y,nKnots,nOctants,periodic,cubicInterpolation);
  
  % Attention: max(iphi) and max(ithe) must be integers, otherwise
  % interp2 returns NaNs.
  if nOctants<8
    iphi = 1 + nOctants*(nKnots-1) * (phi/maxPhi);
    ithe = 1 + (nKnots-1) * (the/(pi/2));
  else
    iphi = 1 + 4*(nKnots-1) * (phi/maxPhi);
    ithe = 1 + 2*(nKnots-1) * (the/pi);
  end
  
  if ~cubicInterpolation
    yi = interp2(z,iphi,ithe,'linear');
  elseif globalInterpolation
    yi = esspline2d(z,ithe,iphi);
  else
    error('Local cubic interpolator for %d (oct) not available!',nOctants);
  end
  
end

return


%--------------------------------------------------------------------------
function z = rectify(y,nKnots,nOctants,periodic,cubicInterp)
%--------------------------------------------------------------------------
%
%  y(iKnots) --> interpolation along
%                constant theta  --> z(theta,phi)
%
% Converts triangularly gridded data y to a rectangular grid z.
%
% nKnots    number of knots between north pole and equator
% periodic  flag telling if data are periodic or not
% nOctants  number of octants the triangular grid covers (see sphgrid)
%
% theta is the the first index z(theta,phi).

fullSphere = nOctants==8;

% nr = number of rows (theta) in output array
% nc = number of columns (phi) in output array (includes end point)
% dlen = change in number of points from one theta row to the next
if fullSphere
  nr = 2*nKnots-1;
  nc = 4*(nKnots-1) + 1;
  dlen = 4;
else
  nr = nKnots;
  nc = nOctants*(nKnots-1) + 1;
  dlen = nOctants;
end

z = zeros(nr,nc);

% add north pole
z(1,:) = y(1);

% process rest of upper hemisphere
idx = 2;
len = dlen + 1 - periodic;
for ir = 2:nKnots-1
  if cubicInterp
    if periodic
      z(ir,:) = esspline1d(y([idx:idx+len-1 idx]),2,nc);
    else
      z(ir,:) = esspline1d(y(idx:idx+len-1),1,nc);
    end
  else
    if periodic
      z(ir,:) = fastlinearinterp1d(y([idx:idx+len-1 idx]),linspace(1,len,nc));
    else
      z(ir,:) = fastlinearinterp1d(y(idx:idx+len-1),linspace(1,len,nc));
    end
  end
  idx = idx + len;
  len = len + dlen;
end

% add equator (no interpolation needed)
if periodic
  z(nKnots,:) = y([idx:idx+len-1 idx]);
else
  z(nKnots,:) = y(idx:idx+len-1);
end

% process lower hemisphere
if fullSphere
  idx = idx + len;
  len = len - dlen;
  for ir = nKnots+1:2*nKnots-2
    if cubicInterp
      if periodic
        z(ir,:) = esspline1d(y([idx:idx+len-1 idx]),2,nc);
      else
        z(ir,:) = esspline1d(y(idx:idx+len-1),1,nc);
      end
    else
      if periodic
        z(ir,:) = fastlinearinterp1d(y([idx:idx+len-1 idx]),linspace(1,len,nc));
      else
        z(ir,:) = fastlinearinterp1d(y(idx:idx+len-1),linspace(1,len,nc));
      end
    end
    idx = idx + len;
    len = len - dlen;
  end
  
  % add south pole
  z(2*nKnots-1,:) = y(idx);
end

return

%--------------------------------------------------------------------------
function yy = fastlinearinterp1d(y,xx)
%--------------------------------------------------------------------------
% Linearly interpolates 1D vector y defined over 1:length(y) at values xx.
k = min(max(1+floor(xx-1),1),length(y)-1);
yy = y(k) + (xx-k).*(y(k+1)-y(k));
return

