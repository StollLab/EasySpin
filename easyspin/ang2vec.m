% ang2vec  Polar angles to cartesian coordinates 
%
%   V = ang2vec(phi,theta)
%   [x,y,z] = ang2vec(phi,theta)
%
%   Converts angles to cartesian unit vector.
%   Angles have to be in radians. V is a 3xN
%   matrix containing the cartesian vectors
%   along columns. theta is the angle down
%   from the z axis, phi is the angle in the xy
%   plane between the x axis and the vectors xy
%   projection (counterclockwise). phi and
%   theta must be arrays with the same number of
%   elements.

function varargout = ang2vec(phi,theta)

if (nargin==0), help(mfilename); return; end

nphi = numel(phi);
ntheta = numel(theta);

if (nphi~=1) && (ntheta~=1) && (nphi~=ntheta)
  error('phi and theta must have the same number of elements!');
end

if (nphi==1)
  phi = phi(ones(1,ntheta));
elseif (ntheta==1)
  theta = theta(ones(1,nphi));
end

% Convert both angles to row vectors.
p = phi(:).';
t = theta(:).';

% Compute cartesian column vectors.
sintheta = sin(t);
x = sintheta.*cos(p);
y = sintheta.*sin(p);
z = cos(t);
switch nargout
case 1
  varargout = {[x;y;z]};
case 3
  varargout = {x,y,z};
end

return
