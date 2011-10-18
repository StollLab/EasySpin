% sphrand  Generate random points over unit sphere
%
%   vecs = sphrand(N)
%   vecs = sphrand(N,k)
%   [phi,theta] = sphrand(...)
%   [x,y,z] = sphrand(...)
%
%   Generates N randomly distributed point over k
%   octants of the unit sphere. The underlying
%   distribution is uniform.
%
%   Input
%   - N: number of points
%   - k: number of octants, can be 1, 2, 4 or 8.
%          1   x>0,y>0,z>0
%          2   y>0,z>0
%          4   z>0 (upper hemisphere)
%          8   full sphere
%        4 is the default.
%
%    Output
%   - vecs: 3xN array with column vectors
%   - phi, theta: polar angles in radians

function varargout = sphrand(N,nOctants)

if (nargin==0), help(mfilename); return; end

% N: number of points
% nOctants: number of octants

if (nargin==1), nOctants = 4; end

switch nOctants
  case {1,2,4}
    z = rand(1,N);
    phi = nOctants*pi/2*rand(1,N);
  case 8
    z = 2*rand(1,N) - 1;
    phi = 2*pi*rand(1,N);
  otherwise
    error('Number of octants is %d, but must be 1, 2, 4 or 8!',nOctants);
end

switch nargout
  case 2
    varargout = {phi,acos(z)};
  case 1
    r = sqrt(1-z.^2);
    varargout = {[r.*cos(phi); r.*sin(phi); z]};
  case 3
    r = sqrt(1-z.^2);
    varargout = {r.*cos(phi), r.*sin(phi), z};
  case 0
    r = sqrt(1-z.^2);
    x = r.*cos(phi);
    y = r.*sin(phi);
    plot3(x,y,z,'.');
    axis equal
    xlabel('x');
    ylabel('y');
    zlabel('z');
  otherwise
    error('Wrong number of output arguments!');
end

return
