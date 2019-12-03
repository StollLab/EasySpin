% rotplane   Set of vectors/orientations in a plane
%
%     v = rotplane(n,chi)
%     v = rotplane(n,chi12,nChi)
%     [phi,theta] = ...
%
%     Computes vectors in a plane perpendicular to
%     a given direction n.
%
%     n      rotation axis, either 3-element vector, or shortcut
%              'x','y','z','xy','xz','yz','xyz',
%              '-x','-y','-z','-xy','-xz','-yz','-xyz'
%     chi    array of plane angles, in radians
%     chi12  [chi1 chi2], start and end angle
%            in plane, in radians
%            rotplane(n,[chi1 chi2],nChi) corresponds to
%            rotplane(n,linspace(chi1,chi2,nChi))
%    
%     v      array of vectors giving the orientations
%            in the plane; one vector per column
%     phi,theta   angles describing the orientations

function varargout = rotplane(n,chi,nChi)

if nargin==0, help(mfilename); return; end

switch nargin
  case 1
    error('Second input argument (chi) is missing.');
  case 2
    % nothing to be done
  case 3
    chi = linspace(chi(1),chi(end),nChi);
  otherwise
    error('Needs 2 or 3 input parameters.');
end

if ischar(n)
  n = letter2vec(n);
else
  if numel(n)~=3
    error('n must be a 3-element vector or a shortcut such as ''x'', ''y'', etc.');
  end
end

z = n(:);

normz = norm(z);
if normz==0
  error('n must be a vector with nonzero length.');
end

Sign = -1;

z = z(:)/normz;
y = [-z(2); z(1); 0];
if norm(y)<eps
  y = [0; 1; 0];
else
  y = y/norm(y);
end
x = cross(y,z);
v = x*cos(chi) + Sign*y*sin(chi);

switch nargout
  case 0
    varargout = {v};
  case 1
    varargout = {v};
  case 2
    [phi,theta] = vec2ang(v);
    varargout = {phi,theta};
end

return
