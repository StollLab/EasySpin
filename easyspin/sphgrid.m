% sphgrid  Spherical triangular grids
%
%   vecs = sphgrid(Symmetry,nKnots)
%   [vecs,Weights] = sphgrid(Symmetry,nKnots,'c')
%   [phi,theta] = sphgrid(Symmetry,nKnots)
%   [phi,theta,Weights] = sphgrid(Symmetry,nKnots)
%
%   Returns a set of unique orientations together with
%   the covered solid angle by each, for a given
%   symmetry.
%
%   Input:
%   - Symmetry: the point group in Schoenflies notation
%     ('Ci', 'Dinfh', 'C2h', etc.)
%   - nKnots: number of knots along a quarter of a
%     meridian, must be at least 2
%
%   Output:
%   - vecs: an 3xm array of orientations (unit column vectors)
%   - Weights: associated solid angles, 1xm vector, sum is 4*pi
%   - phi,theta: polar angles in radians
%
%   Only centrosymmetric point groups and C1 are supported.

% This grid is often called SOPHE grid, after
%   Wang/Hanson JMRA 117 1-8.
% However, it was first published in
%   Kurihara, Monthly Weather Review 93(7), 399-415 (July 1965)

%------------------------------------------------
% Undocumented: 'f' option, and Symmetry as nOct.
%------------------------------------------------

% Weights are approximate.

function varargout = sphgrid(Symmetry,nKnots,Options)

if nargin==0, help(mfilename); return; end

if nargin<2 || nargin>3, error('Wrong number of input arguments!'); end

% Determine computation options
%-------------------------------------------------------------------------------
% Cartesian: return vectors instead of angles
% explicitClosedPhi: force closed phi interval computation
% ComputeWeights: compute weights or not

if nargin<3
  Cartesian = 0;
  explicitClosedPhi = false;
else
  Cartesian = any(strfind(Options,'c'));
  explicitClosedPhi = any(strfind(Options,'f'));
end
computeWeights = nargout>2 || (Cartesian && (nargout==2));

% Initializations
%-------------------------------------------------------------------------------

% Determine maximum phi and open/closed interval from symmetry
% specification.

if ischar(Symmetry)
  [maxPhi,openPhi,nOctants] = symparam(Symmetry);
  FullSphere = strcmp(Symmetry,'C1');
else
  nOctants = Symmetry;
  FullSphere = false;
  switch nOctants
    case {-1, 0}, maxPhi = 0;
    case 1, maxPhi = pi/2;
    case 2, maxPhi = pi;
    case 4, maxPhi = 2*pi;
    case 8, maxPhi = 2*pi; FullSphere = true;
    otherwise
      error('Wrong number of octants, must be -1, 0, 1, 2, 4 or 8.');
  end
  openPhi = 1;
end

% Force open phi interval if input parameter is set.
if explicitClosedPhi
  openPhi = 0;
end

dtheta = (pi/2)/(nKnots-1);

% Now we have maxPhi, openPhi, dtheta, nOctants and nKnots,
% completely specifying the grid we want.

if nOctants > 0 % if not Dinfh or O3 symmetry
  
  % Initializations
  nOct = ceil(maxPhi/(pi/2));
  nPoints = nKnots + nOct*nKnots*(nKnots-1)/2;
  theta = zeros(1,nPoints);
  phi = zeros(1,nPoints);
  Weights = zeros(1,nPoints);
  
  if computeWeights
    sindth2 = sin(dtheta/2);
    w1 = 0.5 + openPhi/2;
  end
  
  % North pole
  if computeWeights
    Weights(1) = maxPhi*(1-cos(dtheta/2));
  end
  
  % All but equatorial slice
  Start = 2;
  for iSlice = 2:nKnots-1
    nPhi = nOct*(iSlice-1)+1;
    dPhi = maxPhi/(nPhi-1);
    idx = Start+(0:nPhi-1);
    if computeWeights
      Weights(idx) = 2*sin((iSlice-1)*dtheta)*sindth2*dPhi*...
        [w1 ones(1,nPhi-2) .5];
    end
    theta(idx) = (iSlice-1)*dtheta;
    phi(idx) = linspace(0,maxPhi,nPhi);
    Start = Start + nPhi;
  end
  
  % Equatorial slice
  nPhi = nOct*(nKnots-1)+1;
  dPhi = maxPhi/(nPhi-1);
  idx = Start+(0:nPhi-1);
  theta(idx) = pi/2;
  phi(idx) = linspace(0,maxPhi,nPhi);
  if (computeWeights)
    Weights(idx) = sindth2*dPhi*[w1 ones(1,nPhi-2) .5];
  end
  
  % Border removal
  if openPhi
    openRemove = cumsum(nOct*(1:nKnots-1)+1)+1;
    phi(openRemove) = [];
    theta(openRemove) = [];
    Weights(openRemove) = [];
  end
  
  % Add lower hemisphere for C1
  if FullSphere
    idx = length(theta)-nPhi+openPhi:-1:1;
    phi = [phi phi(idx)];
    theta = [theta pi-theta(idx)];
    if (computeWeights)
      Weights(idx) = Weights(idx)/2; % half of real value
      Weights = [Weights Weights(idx)];
    end
  end
  
  if computeWeights
    Weights = 2*(2*pi/maxPhi)*Weights;
  end
  
elseif nOctants==0 % Dinfh symmetry
  
  phi = zeros(1,nKnots);
  theta = linspace(0,pi/2,nKnots);
  Weights = -2*(2*pi)*diff(cos([0 dtheta/2:dtheta:pi/2 pi/2]));
  
else % O3 symmetry
  
  phi = 0;
  theta = 0;
  Weights = 4*pi;
  
end


% output
%-------------------------------------------------------------------------------
switch nargout
  case 1
    varargout = {ang2vec(phi,theta)};
  case 2
    if Cartesian
      varargout = {ang2vec(phi,theta),Weights};
    else
      varargout = {phi,theta};
    end
  case 3
    varargout = {phi,theta,Weights};
  case 0
    v = ang2vec(phi,theta);
    plot3(v(1,:),v(2,:),v(3,:),'.')
    axis equal
    xlabel('x');
    ylabel('y');
    zlabel('z');
  otherwise
    error('Wrong number of output arguments!');
end

return
