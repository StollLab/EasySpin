% esintpol  Interpolates data of a given symmetry. 
%
%   yi = esintpol(y,Sym,Factor,Opt)
% 
%   y:            Data, in row vectors. Matrices get interpolated along rows.
%   Sym:          [nKnots, periodic, nOctants], as computed by symparam()
%   Opt:         'g3' (default),'l3','l1'

function yi = esintpol(y,Symmetry,Factor,Options,phi,the)

% special case
%if Factor==1, error('Interpolation factor must be bigger than 1!'); end
%if Factor==1, yi=y; end

% Set default or user-defined parameters
if nargin<4 || isempty(Options)
  Options = 'G3';
end

% parse parameters
Global = upper(Options(1)) == 'G';
Cubic = Options(2) == '3';

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

if ischar(Symmetry)
  error('Don''t supply symmetry string!');
else
  % maxPhi = Symmetry(1); %not needed
  nKnots = Symmetry(1);
  periodic = Symmetry(2);
  nOctants = Symmetry(3);
end
 
if nOctants>0 && nargin<3
  error('Not enough input parameters!');
end

% Compute number of slices
if nOctants>0
  nExpectedData = nKnots*(nKnots-1)/2*nOctants + 1;
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

% If interpolation points are not given, compute them.
if nargin<5
  if nOctants<4, closedPhi = 'f'; else, closedPhi = []; end
  grid = sphgrid(Symmetry,(nKnots-1)*Factor+1,closedPhi);
  phi = grid.phi.';
  the = grid.the.';
end

switch nOctants
case 0, % Dinfh
  %============================================================
  % Dinfh
  %============================================================
  % Linear interpolation, vectorized
  % Cubic interpolation, local and global, vectorized
  % Slope estimators: simple and Fritsch-Carlson
  %============================================================
  [m,n] = size(y);
  
  if (~Cubic) % linear interpolation
    yi = interp1(y.',1:1/Factor:n);
    
  elseif (Global) % global cubic interpolation
    yi = esspline1d(y,1,1:1/Factor:n);
    
  else % local cubic interpolation
    % prepare [x^3 x^2 x 1] for each point
    x = (0:1/Factor:1-1/Factor).';
    x = [x.^3 x.^2 x ones(Factor,1)];
    % coefficient matrix for Hermite interpolation
    H = [2 -2 1 1; -3 3 -2 -1; 0 0 1 0; 1 0 0 0];
    
    Estimator = 2;
    % run over all row vectors in y
    for r = 1:m
      switch Estimator
      case 1 % simple slope average
        Tangents = [0 (y(r,3:end)-y(r,1:end-2))/2 0];
      case 2 % Fritsch-Carlson monotone slopes
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
  z = rectify(y,nKnots,nOctants,periodic,~Cubic);
  
  % Attention: max(iphi) and max(ithe) must be integers, otherwise
  % interp2 returns NaNs.
  iphi = 1 + nOctants*(nKnots-1) * (phi/phi(end));
  ithe = 1 + (nKnots-1) * (the/the(end));
  
  if ~Cubic
    yi = interp2(z,iphi,ithe,'*linear');
  elseif Global
    yi = esspline2d(z,ithe,iphi);
  else
    error('Local cubic interpolator for %d (oct) not available!',nOctants);
  end
  
end

return


%----------------------------------------------------------
function z = rectify(y,nr,nOctants,periodic,linear)
%----------------------------------------------------------
%
%  y(iKnots) --> interpolation along
%                constant theta  --> z(theta,phi)
%
% Converts triangularly gridded data y to rectangular z.
% nr is the number of theta slices, periodic is a flag
% telling if data are periodic or not. nOctants gives the number
% of octants the triangular grid covers (see sphgrid).
% theta is the the first index z(theta,phi).

nc = nOctants*(nr-1) + 1 - periodic;
pos = 2;
len = nOctants + 1 - periodic;
if periodic
  z = zeros(nr,nc+1);
  z(1,:) = y(1);
  for i = 2:nr-1
    if linear
      z(i,:) = fastlinearinterp1d(y([pos:pos+len-1 pos]),linspace(1,len,nc+1));
    else
      z(i,:) = esspline1d(y([pos:pos+len-1 pos]),2,nc+1);
    end
    pos = pos + len;
    len = len + nOctants;
  end
  z(nr,:) = y([pos:end pos]);
else
  z = zeros(nr,nc);
  z(1,:) = y(1);
  for i = 2:nr-1
    if linear
      z(i,:) = fastlinearinterp1d(y(pos:pos+len-1),linspace(1,len,nc));
    else
      z(i,:) = esspline1d(y(pos:pos+len-1),1,nc);
    end
    pos = pos + len;
    len = len + nOctants;
  end
  z(nr,:) = y(pos:end);
end

return

%-------------------------------------------------------------------------------
function yy = fastlinearinterp1d(y,xx)
%-------------------------------------------------------------------------------
% Lienarly interpolates 1D vector y defined over 1:length(y) at values xx.
k = min(max(1+floor(xx-1),1),length(y)-1);
yy = y(k) + (xx-k).*(y(k+1)-y(k));
return


% esintpol test area
%===============================================================================

% Dinfh symmetry

opt.Scope = 'local';
opt.Order = 3;
the = linspace(0,pi/2,10);
y = [sin(the).^2; cos(the).^2+.1];
for k = 1:10
  yy = esintpol(y,'Dinfh',k);
  fthe = linspace(0,pi/2,length(yy));
  plot(fthe,yy,'.-b',the,y,'or');
  pause
end
close

% C2h symmetry

Symmetry = 'C2h';

opt.Scope = 'global';
opt.Order = 3;
% y = [1 2 3 3 3.5 4]; % D2h
y = [1 1 2 1 2 3 4];
for k = 1:10
  grid = sphgrid(Symmetry,2*k+1);
  p = grid.phi;
  t = grid.theta;
  yy = esintpol(y,Symmetry,k,opt,p,t);
  showdata(yy,Symmetry,2*k+1);
  pause
end
close

% C4h symmetry

Symmetry = 'C4h';

opt.Scope = 'global';
opt.Order = 3;
y = [1 1 1 2 1 2 3];
n = 4;
for k = 1:10
  kk = (n-1)*k+1;
  grid = sphgrid(Symmetry,kk);
  p = grid.phi;
  t = grid.theta;
  yy = esintpol(y,Symmetry,k,opt,p,t);
  showdata(yy,Symmetry,kk);
  pause
end
close

