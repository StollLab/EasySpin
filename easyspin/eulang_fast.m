% eulang_fast  Calculate Euler angles from rotation matrix
%
%   [alpha,beta,gamma] = eulang_fast(R)
%
%   Returns the three Euler angles alpha, beta, and gamma
%   (in radians) of the rotation matrix R, which must be
%   a 3x3 real matrix with determinant very close to +1.
%   No validity checks are performed on R.
%
%   For a definition of the Euler angles, see erot().
%
%   [alpha,beta,gamma] and [alpha+-pi,-beta,gamma+-pi]
%   give the same rotation matrix. eulang_fast() returns the
%   set with beta>=0.

function [alpha,beta,gamma] = eulang_fast(R)

% No validity checks on R

% Calculate Euler angles using analytical expressions
%-------------------------------------------------------------------------------
% Degenerate cases:  R(3,3) = cos(beta) = +-1 -> beta=n*pi. In these
% cases, alpha and gamma are not distinguishable. We collect the entire
% z rotation angle in alpha and set gamma to zero.
degenerateCaseLimit = 1e-8;
if abs(R(3,3)-1)<=degenerateCaseLimit
  alpha = atan2(R(1,2),R(2,2));
  beta = 0;
  gamma = 0;
elseif abs(R(3,3)+1)<=degenerateCaseLimit
  alpha = atan2(-R(1,2),R(2,2));
  beta = pi;
  gamma = 0;
else
  alpha = atan2(R(3,2),R(3,1));
  beta = atan2(sqrt(R(3,1)^2+R(3,2)^2),R(3,3));
  gamma = atan2(R(2,3),-R(1,3));
end

% Assure alpha and gamma are positive (unless numerically close to zero).
%-------------------------------------------------------------------------------
negativeThreshold = -1e-8;
if alpha<negativeThreshold, alpha = alpha + 2*pi; end
if gamma<negativeThreshold, gamma = gamma + 2*pi; end

% Assure beta is positive
%-------------------------------------------------------------------------------
if beta<0
  beta = -beta;
  if alpha<0
    alpha = alpha + pi;
    gamma = gamma + pi;
  else
    alpha = alpha - pi;
    gamma = gamma - pi;
  end
end

end
