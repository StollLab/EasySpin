% rotmat2axi   Rotation axis and angle form rotation matrix
%
%    [n,rho] = rotmat2axi(R)
%
%    Determines the axis n (a 3-element vector) and the angle
%    rho (in  radians) of the rotation described by the 3x3
%    rotation matrix R.
%
%    Example:
%      % A rotation matrix describing a rotation around y
%      R = erot(0,pi/3,0);
%      % rotmax2axi yields axis y and angle pi/3
%      [n,rho] = rotmat2axi(R)

function [n,rho] = rotmat2axi(R)

if nargin==0, help(mfilename); return; end

if ~isnumeric(R) || ~isreal(R) || any(size(R)~=[3 3])
  error('R must be a 3x3 real matrix');
end

% trace(R) corresponds to (cos(beta) + 1) * (cos(alpha+gamma) - 1)

if trace(R)==3
  % no rotation, rho==0
  rho = 0;
  n = [0;0;1];
else
  rho = real(acos((trace(R)-1)/2));
  n = [R(2,3)-R(3,2); R(3,1)-R(1,3); R(1,2)-R(2,1)];
  n = n/norm(n);
end

return
