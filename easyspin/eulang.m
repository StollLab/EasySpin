% eulang  Euler angles from rotation matrix
%
%   angles = eulang(R)
%   [alpha,beta,gamma] = eulang(R)
%
%   Returns the three Euler angles alpha, beta and gamma (in radians) of the
%   rotation matrix R, which must be a 3x3 real matrix with determinant very
%   close to +1.
%
%   [alpha,beta,gamma] and [alpha+-pi,-beta,gamma+-pi]
%   give the same rotation matrix. eulang() returns the
%   set with beta>=0.
%
%   If the matrix is close to orthogonal, a neighoring orthogonal matrix is
%   calculated using singular-value decomposition.

% The second input, nocheck, is undocumented. If set to true, all validity
% checks are bypassed. This improves performance, since the orthogonality
% check is expensive.

function varargout = eulang(R,nocheck)

if nargin<1, help(mfilename); return; end
if nargin<2, nocheck = false; end

if ~nocheck

  % Check size and real-valuedness
  %-------------------------------------------------------------------------------
  if any(size(R)~=3) || ~isreal(R)
    error('eulang: Rotation matrix must be a real-valued 3x3 matrix.');
  end

  % Check orthogonality
  %-------------------------------------------------------------------------------
  % (The determinant of an orthogonal matrix is +-1, but the converse is not true.)
  % Check orthonormality of columns
  orthogonalityError = norm(R'*R-eye(3));
  if orthogonalityError>1e-2
    error('eulang: Rotation matrix is not orthogonal, deviation is %f.',orthogonalityError);
  elseif orthogonalityError>1e-8
    fprintf('eulang: Rotation matrix is not orthogonal, deviation is %g.\n',orthogonalityError);
    fprintf('eulang: Orthogonalizing using singular-value decomposition (SVD).\n');
    [U,~,V] = svd(R);
    R = U*V.';
  end

  % Check sign of determinant
  %-------------------------------------------------------------------------------
  d = det(R);
  if d<0
    error('eulang: Rotation matrix has negative determinant. Change the signs in one column or row.');
  end

end


% Calculate Euler angles using analytical expressions
%-------------------------------------------------------------------------------
% Degenerate cases:  R(3,3) = cos(beta) = +-1 -> beta=n*pi. In these
% cases, alpha and gamma rotations are not separable. We collect the entire
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
thr = -1e-8;
if alpha<thr, alpha = alpha + 2*pi; end
if gamma<thr, gamma = gamma + 2*pi; end


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


% Collect output
%-------------------------------------------------------------------------------
switch nargout
  case {0,1}
    angles = [alpha,beta,gamma];
    varargout = {angles};
  case 3
    varargout = {alpha,beta,gamma};
  otherwise
    error('Wrong number of output arguments.')
end
