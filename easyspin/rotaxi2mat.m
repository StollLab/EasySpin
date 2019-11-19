% rotaxi2mat   Compute rotation matrix from rotation axis and rotation angle
%
%    R = rotaxi2mat(n,rho)
%
%    Generates the matrix representing the rotation around the axis given by
%    the 3-element vector n by the angle rho (in radians). n does not have to
%    be normalized.
%
%    There are shortcuts for a few specific directions of the rotation axis:
%      n='x' implies n=[1;0;0]
%      n='y'         n=[0;1;0]
%      n='z'         n=[0;0;1]
%      n='xy'        n=[1;1;0]
%      n='xz'        n=[1;0;1]
%      n='yz'        n=[0;1;1]
%      n='xyz'       n=[1;1;1]
%
%    Example:
%      % A rotation by 2*pi/3 (120 degrees) around the axis [1;1;1]
%      % permutes the x, y and z axes.
%      R = rotaxi2mat('xyz',2*pi/3)

function R = rotaxi2mat(n,rho)

if nargin==0, help(mfilename); return; end

if nargin<2
  error('Provide a second input argument (rotation angle rho).')
end

if ~isreal(rho) || numel(rho)~=1
  error('Second input (angle rho) must be a real number.');
end

if ischar(n)
  n = letter2vec(n);
else
  if numel(n)~=3 || ~isreal(n)
    error('Rotation axis must be a 3-element vector or a shortcut such as ''x'', ''y'', etc.');
  end
end

n = n/norm(n);

N = zeros(3,3);

% N_{ij} = -e_{ijk} n_k, with
%                / +1 if ijk is an even permutation of 123
%      e_{ijk} = | -1 if ijk is an odd permutation of 123
%                \  0 otherwise (if two indices are equal)
N([8 3 4]) = n;
N = N - N.';

R = eye(3) + N*sin(rho) + N^2*(1-cos(rho));
%R = expm(phi*N); % alternative

% Remove numerical errors for entries with 0, +1 and -1
thresh = 1e-10;
R(abs(R)<thresh) = 0;
R(abs(R-1)<thresh) = +1;
R(abs(R+1)<thresh) = -1;

return
