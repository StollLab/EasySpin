% wignerd  Compute Wigner D-matrix
%
%   D = wignerd(J,angles);
%   D = wignerd(J,angles,'-');
%   D = wignerd(J,angles,'+');
%
%   This function computes the Wigner rotation matrix D^J_(m1,m2)
%   where m1 and m2 run from J to -J.
%
%   angles ... [alpha beta gamma], Euler angles in radians, or
%              beta, second Euler angle in radians
%   J      ... rank
%   D      ... Wigner rotation matrix, size (2J+1)x(2J+1)
%
%   If angles contains only one number, the Wigner small-d matrix
%   d^J_(m1,m2) is computed.
%
%   If '+' is given (or the third argument is omitted), wignerd computes
%   the (2J+1)x(2J+1) matrix
%     <Jm1|expm(1i*alpha*Jz)*expm(1i*beta*Jy)*expm(1i*gamma*Jz)|Jm2>
%     (sign convention as in Edmonds, Mathematica)
%   If '-' is given, wignerd returns the (2J+1)x(2J+1) matrix representation of
%     <Jm1|expm(-1i*alpha*Jz)*expm(-1i*beta*Jy)*expm(-1i*gamma*Jz)|Jm2>
%     (sign convention as in Brink/Satcher, Zare, Sakurai, Varshalovich)
%
%   The basis is ordered +J..-J from left to right and from top to bottom,
%   so that e.g. the output matrix element D(1,2) corresponds to D^J_(J,J-1).

function D = wignerd(J,angles,phase)

switch nargin
  case 0, help(mfilename); return;
  case 1, error('Angles are missing');
  case 2, phase = '-';
end

switch numel(angles)
  case 1, % only beta given
  case 3, % [alpha beta gamma]
  otherwise
    error('Angles must be either [alpha beta gamma] or beta.')
end

if numel(J)~=1 || ~isreal(J) || mod(J,0.5) || (J<0) || ~isnumeric(J)
  error('J must be one of 0, 1/2, 1, 3/2, 2, etc');
end

if ~ischar(phase) || (numel(phase)~=1)
  error('Third argument must be either ''+'' or ''-''.');
end

if (phase=='+')
elseif (phase=='-')
  angles = -angles;
else
  error('Third argument must be either ''+'' or ''-''.');
end

if ~any(angles)
  D = eye(2*J+1);
  return;
end

v = sqrt((1:2*J).*(2*J:-1:1))/2; % off-diagonals of Jy matrix
Jy = diag(v,+1) - diag(v,-1); % assemble Jy matrix
if numel(angles)==1
  D = expm(angles*Jy);
else
  mz = J:-1:-J; % diagonal of Jz matrix
  D = exp(1i*angles(1)*mz.') * exp(1i*angles(3)*mz);
  D = D.*expm(angles(2)*Jy);
end
