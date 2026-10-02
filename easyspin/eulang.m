% eulang  Euler angles from rotation matrix
%
%   angles = eulang(R)
%   [alpha,beta,gamma] = eulang(R)
%
%   Returns the three Euler angles alpha, beta and gamma (in radians) of the
%   rotation matrix R, which must be a real orthogonal 3x3 matrix with
%   determinant +1.
%
%   alpha and gamma are in the range [0,2*pi), and beta is in [0,pi].
%   If beta is 0 or pi, alpha and gamma cannot be separated; the entire
%   z rotation is then returned in alpha, and gamma is set to zero.
%
%   If the matrix is close to orthogonal, a neighboring orthogonal matrix is
%   calculated using singular-value decomposition.

% The second input, nocheck, is undocumented. If set to true, all validity
% checks are bypassed. This improves performance, since the orthogonality
% check is expensive.

function varargout = eulang(R,nocheck)

if nargin<1, help(mfilename); return; end
if nargin<2, nocheck = false; end

% Thresholds
%-------------------------------------------------------------------------------
orthErrorLimit = 1e-2;        % orthogonality error above which eulang errors
orthogonalizeLimit = 1e-6;    % orthogonality error above which R is orthogonalized
degenerateCaseLimit = 1e-14;  % sin(beta) below which beta is taken as 0 or pi
negativeAngleLimit = 1e-8;    % alpha, gamma above -negativeAngleLimit are not shifted by 2*pi

if ~nocheck

  % Check size and real-valuedness
  %-------------------------------------------------------------------------------
  if any(size(R)~=3) || ~isreal(R) || ~all(isfinite(R(:)))
    error('eulang: Rotation matrix must be a real-valued 3x3 matrix.');
  end

  % Check orthogonality
  %-------------------------------------------------------------------------------
  % (The determinant of an orthogonal matrix is +-1, but the converse is not true.)
  % Check orthonormality of columns
  orthogonalityError = norm(R.'*R-eye(3));
  if orthogonalityError>orthErrorLimit
    error('eulang: Rotation matrix is not orthogonal, deviation is %g.',orthogonalityError);
  end

  % Check sign of determinant
  %-------------------------------------------------------------------------------
  if det(R)<0
    error('eulang: Rotation matrix has negative determinant. Change the signs in one column or row.');
  end

  % Orthogonalize if needed
  %-------------------------------------------------------------------------------
  % Construct the closest orthogonal rotation matrix using
  % singular-value decomposition.
  if orthogonalityError>orthogonalizeLimit
    fprintf('eulang: Rotation matrix is not orthogonal, deviation is %g.\n',orthogonalityError);
    fprintf('eulang: Orthogonalizing using singular-value decomposition (SVD).\n');
    [U,~,V] = svd(R);
    R = U*V.';  % has same sign of determinant as R
  end

end


% Calculate Euler angles using analytical expressions
%-------------------------------------------------------------------------------
% Degenerate case (beta = 0 or pi): alpha and gamma are not separable, so the
% entire z rotation goes into alpha and gamma is set to zero.
% Non-degenerate case: alpha from R(3,1:2) is ill-conditioned for small
% sin(beta), so gamma is computed from alpha+gamma (beta<pi/2) or gamma-alpha
% (beta>pi/2), which are well-conditioned in R(1:2,1:2). Errors in alpha then
% do not affect the reconstructed rotation matrix.
sinbeta = hypot(R(3,1),R(3,2));
if sinbeta<=degenerateCaseLimit
  if R(3,3)>0
    alpha = atan2(R(1,2)-R(2,1),R(1,1)+R(2,2));  % alpha+gamma, with gamma = 0
    beta = 0;
  else
    alpha = atan2(-(R(1,2)+R(2,1)),R(2,2)-R(1,1));  % alpha-gamma, with gamma = 0
    beta = pi;
  end
  gamma = 0;
else
  alpha = atan2(R(3,2),R(3,1));
  beta = atan2(sinbeta,R(3,3));  % always in [0,pi]
  if R(2,3)==0
    gamma = atan2(0,-R(1,3));  % R(2,3) = sin(gamma)*sin(beta), so gamma is exactly 0 or pi
  elseif R(3,3)>0
    gamma = atan2(R(1,2)-R(2,1),R(1,1)+R(2,2)) - alpha;
  else
    gamma = atan2(R(1,2)+R(2,1),R(2,2)-R(1,1)) + alpha;
  end
  gamma = gamma - 2*pi*round(gamma/(2*pi));  % wrap to [-pi,pi]
end


% Assure alpha and gamma are positive (unless numerically close to zero).
%-------------------------------------------------------------------------------
if alpha<-negativeAngleLimit, alpha = alpha + 2*pi; end
if gamma<-negativeAngleLimit, gamma = gamma + 2*pi; end


% Collect output
%-------------------------------------------------------------------------------
switch nargout
  case {0,1}
    angles = [alpha,beta,gamma];
    varargout = {angles};
  case 3
    varargout = {alpha,beta,gamma};
  otherwise
    error('eulang: Wrong number of output arguments.')
end
