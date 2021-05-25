% gridinterp  Interpolates spherical-grid data
%
%   yi = gridinterp(y,gridParams,phi,theta,InterpType)
%   
% Inputs:
%   y            data, in row vectors. Matrices get interpolated along rows.
%   gridParams   [nOrientations, GridSize, closedPhi, nOctants, maxPhi]
%   phi,theta    grid points for interpolation
%   InterpType   interpolation type, 'G3' (default),'L3','L1'
% Outputs:
%   yi           y interpolated over interpolation points

% List of unique phi intervals as set by sphgrid(), octant numbers, and end
% conditions (p is periodic, z is zero endslopes).

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

function yi = gridinterp(y,gridParams,phi,theta,InterpType)

if nargin<4 || nargin>5
  error('4 or 5 input arguments ar required.')
end

% Set default or user-defined parameters
if nargin<5 || isempty(InterpType)
  InterpType = 'G3';
end

% Parse parameters
globalInterpolation = upper(InterpType(1)) == 'G';
cubicInterpolation = InterpType(2) == '3';

% Extract grid parameters
nOrientations = gridParams(1);
GridSize = gridParams(2);
periodic = ~gridParams(3);
nOctants = gridParams(4);
maxPhi = gridParams(5);
 
if nOctants>0 && nargin<3
  error('Not enough input parameters!');
end

% Make sure data has expected number of orientations
if length(y)~=nOrientations
  error('Wrong number of points. Not a valid data set!');
end

switch nOctants
case 0
  %============================================================
  % Dinfh
  %============================================================
  % Linear interpolation, vectorized
  % Cubic interpolation, local and global, vectorized
  % Slope estimators: simple and Fritsch-Carlson
  %============================================================
  [m,n] = size(y);
  Factor = (numel(phi)-1)/(n-1);
  iphi = 1:1/Factor:n;
  
  if ~cubicInterpolation % linear interpolation
    yi = interp1(y.',iphi);
    
  elseif globalInterpolation % global cubic interpolation
    yi = esspline1d(y,1,iphi);
    
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
  
otherwise
  %============================================================
  % nOctants >=1 (D2h, etc)
  %============================================================
  % Convert triangular to rectangular grid.
  z = rectify(y,GridSize,nOctants,periodic,cubicInterpolation);
  
  % Attention: max(iphi) and max(ithe) must be integers, otherwise
  % interp2 returns NaNs.
  if nOctants<8
    iphi = 1 + nOctants*(GridSize-1) * (phi/maxPhi);
    ithe = 1 + (GridSize-1) * (theta/(pi/2));
  else
    iphi = 1 + 4*(GridSize-1) * (phi/maxPhi);
    ithe = 1 + 2*(GridSize-1) * (theta/pi);
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
function z = rectify(y,GridSize,nOctants,periodic,cubicInterp)
%--------------------------------------------------------------------------
%
%  y(iKnots) --> interpolation along
%                constant theta  --> z(theta,phi)
%
% Converts triangularly gridded data y to a rectangular grid z.
%
% GridSize    number of knots between north pole and equator
% periodic  flag telling if data are periodic or not
% nOctants  number of octants the triangular grid covers (see sphgrid)
%
% theta is the the first index z(theta,phi).

fullSphere = nOctants==8;

% nr = number of rows (theta) in output array
% nc = number of columns (phi) in output array (includes end point)
% dlen = change in number of points from one theta row to the next
if fullSphere
  nr = 2*GridSize-1;
  nc = 4*(GridSize-1) + 1;
  dlen = 4;
else
  nr = GridSize;
  nc = nOctants*(GridSize-1) + 1;
  dlen = nOctants;
end

z = zeros(nr,nc);

% Add north pole
z(1,:) = y(1);

% Process rest of upper hemisphere
idx = 2;
len = dlen + 1 - periodic;
for ir = 2:GridSize-1
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

% Add equator (no interpolation needed)
if periodic
  z(GridSize,:) = y([idx:idx+len-1 idx]);
else
  z(GridSize,:) = y(idx:idx+len-1);
end

% Process lower hemisphere
if fullSphere
  idx = idx + len;
  len = len - dlen;
  for ir = GridSize+1:2*GridSize-2
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
  z(2*GridSize-1,:) = y(idx);
end

return

%--------------------------------------------------------------------------
function yy = fastlinearinterp1d(y,xx)
%--------------------------------------------------------------------------
% Linearly interpolates 1D vector y defined over 1:length(y) at values xx.
k = min(max(1+floor(xx-1),1),length(y)-1);
yy = y(k) + (xx-k).*(y(k+1)-y(k));
return
