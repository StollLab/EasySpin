% rotaxi2mat   Compute rotation matrix from axis and angle
%
%    R = rotaxi2mat(n,phi)
%
%    Generates the matrix representing the rotation
%    around the axis given by the 3-element vector n
%    by the angle phi (in radians).
%
%    There are three shortcuts for x, y, and z axis:
%      n=1 implies n=[1;0;0]
%      n=2 means   n=[0;1;0]
%      n=3 gives   n=[0;0;1].
%
%    Example:
%      % A rotation by 2*pi/3 around the axis [1;1;1]
%      % permutes the x, y and z axes.
%      R = rotaxi2mat([1;1;1],2*pi/3)

function R = rotaxi2mat(n,phi)

if nargin==0, help(mfilename); return; end

if ~isreal(phi) | numel(phi)~=1
  error('phi must be a real number!');
end

if ~isreal(n)
  error('Rotation axis vector n must have real elements!');
end

switch numel(n)
case 1,
  switch n
  case 1, n = [1;0;0];
  case 2, n = [0;1;0];
  case 3, n = [0;0;1];
  otherwise
    error('Rotation axis shortcut is invalid!');
  end
case 3,
  n = n/norm(n);
otherwise
  error('Rotation axis must be a vector or a scalar shortcut!');
end

N = zeros(3,3);

% N_{ij} = -e_{ijk} n_k, e_{ijk} = +1 if ijk even permutation of 123,
% -1 if odd permutation, and 0 otherwise (if two indices are equal)
N([8 3 4]) = n;
N = N - N.';

R = eye(3) + N*sin(phi) + N^2*(1-cos(phi));

% Alternative:
%R = expm(phi*N);

% Remove numerical errors for entries with 0, +1 and -1
thresh = 1e-10;
R(abs(R)<thresh) = 0;
R(abs(R-1)<thresh) = +1;
R(abs(R+1)<thresh) = -1;

return
