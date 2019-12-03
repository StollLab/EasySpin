% vec2ang  Convert cartesian vectors to polar angles 
%
%   [phi,theta] = vec2ang(v)
%   ang = vec2ang(v)
%
%   Converts cartesian vectors in v to polar angles phi and theta.
%
%   Inputs:
%     v can be a 3xN array (list of column vectors) or Nx3 array (list of row
%        vectors). If 3x3, column vectors are assumed.
%     v can be a vector shorthand such as 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'.
%
%   Outputs:
%     theta is the angle down from the z axis to the vector v
%     phi is the angle between the x axis and the projection of the vector v
%        onto the xy plane, anticlockwise.
%
%   If only one output is requested, then an array with the phi values in the
%   first row and the theta values in the secon row is returned (ang = [phi; theta]).
%
%   Angles are in radians. The vectors in v don't have to be normalized.

function [phi,theta] = vec2ang(v)

if nargin==0, help(mfilename); return; end

if ischar(v)
  v = letter2vec(v);
end

% Flip array if it is obviously a list of row vectors.
if size(v,1)~=3 && size(v,2)==3
  v = v.';
end

% Error if the input array has not the correct form.
if size(v,1)~=3
  error('At least one dimension of the input array must be 3.');
end

% Error if the input array doesn't contain real numbers.
if ~isnumeric(v) || ~isreal(v)
  error('Argument must contain real numbers!');
end

% Compute angles using atan2() which correctly sets the signs.
% Vectors DON'T have to be normalized to unit length.
x = v(1,:);
y = v(2,:);
z = v(3,:);
phi = atan2(y,x);
theta = pi/2 - atan2(z,sqrt(x.^2+y.^2));

switch nargout
  case {0,1}, phi = [phi; theta];
  case 2 % phi and theta separate
end

return
