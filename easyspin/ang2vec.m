% ang2vec  Polar angles to cartesian coordinates 
%
%   V = ang2vec(phi,theta)
%   [x,y,z] = ang2vec(phi,theta)
%
%   Converts angles to cartesian unit vector.
%
%   Inputs:
%     theta is the angle down from the z axis.
%     phi is the angle in the xy plane between the x axis and
%       the vectors xy projection (counterclockwise).
%     phi and theta must be arrays with the same number of elements.
%
%   Outputs:
%     Angles have to be in radians. V is a 3xN matrix containing the
%     cartesian vectors along columns.
%     x, y, and z are the separate Cartesian components of V.

function varargout = ang2vec(phi,theta)

if nargin==0, help(mfilename); return; end

if nargout==2 || nargout>3
  error('One or three output arguments are required.');
end

nphi = numel(phi);
ntheta = numel(theta);

if nphi~=1 && ntheta~=1 && nphi~=ntheta
  error('phi and theta must have the same number of elements!');
end

if nphi==1
  phi = phi(ones(1,ntheta));
elseif ntheta==1
  theta = theta(ones(1,nphi));
end

% Convert both angles to row vectors.
phi = phi(:).';
theta = theta(:).';

% Compute cartesian column vectors.
sintheta = sin(theta);
x = sintheta.*cos(phi);
y = sintheta.*sin(phi);
z = cos(theta);

switch nargout
  case 1
    varargout = {[x;y;z]};
  case 3
    varargout = {x,y,z};
end

return
