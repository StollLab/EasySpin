% eulang  Euler angles from rotation matrix
%
%   Angles = eulang(R)
%   [alpha,beta,gamma] = eulang(R)
%
%   Returns the three Euler angles [alpha, beta, gamma]
%   (in radians) of the rotation matrix R, which must be
%   a 3x3 real matrix with determinant very close to +1.
%
%   For a definition of the Euler angles, see erot().
%
%   [alpha,beta,gamma] and [alpha+-pi,-beta,gamma+-pi]
%   give the same rotation matrix. eulang() returns the
%   set with beta>=0.

function varargout = eulang(RotMatrix,skipFitting)

if nargin<1, help(mfilename); return; end
if nargin<2, skipFitting = 0; end

% Validity checks
%--------------------------------------------------------
if ~isreal(RotMatrix)
  error('eulang: Rotation matrix must be real.');
end

if any(size(RotMatrix)~=3)
  error('eulang: Rotation matrix must be 3x3.');
end

normError = norm(RotMatrix'*RotMatrix-eye(3));
if normError>1e-6
  disp('eulang: Rotation matrix is not orthogonal within 1e-6.');
elseif normError>1e-2
  error('eulang: Rotation matrix is not orthogonal within 1e-2.');
end

d = det(RotMatrix);
if d<0
  error('eulang: Rotation matrix has negative determinant! Change the signs in one column or row.');
end

hardLimit = 1e-1;
softLimit = 5e-3;
detError = abs(d-1);
if detError>hardLimit
  error('eulang: Determinant of rotation matrix is %+0.3f, deviates too much from +1!\nRescale argument to R/det(R)^(1/3) if result wanted.',d);
elseif detError>softLimit
  fprintf('eulang: Determinant of rotation matrix is %+0.3f.\n',d);
  fprintf('eulang: The output of eulang() is probably inaccurate!\n');
end

% Calculate Euler angles using analytical expressions
%-------------------------------------------------------------------------------
% Degenerate cases:  R(3,3) = cos(beta) = +-1 -> beta=n*pi. In these
% cases, alpha and gamma are not distinguishable. We collect the entire
% z rotation angle in alpha and set gamma to zero.
degenerateCaseLimit = 1e-8;
if abs(RotMatrix(3,3)-1)<=degenerateCaseLimit
  alpha = atan2(RotMatrix(1,2),RotMatrix(2,2));
  beta = 0;
  gamma = 0;
elseif abs(RotMatrix(3,3)+1)<=degenerateCaseLimit
  alpha = atan2(-RotMatrix(1,2),RotMatrix(2,2));
  beta = pi;
  gamma = 0;
else
  alpha = atan2(RotMatrix(3,2),RotMatrix(3,1));
  beta = atan2(sqrt(RotMatrix(3,1)^2+RotMatrix(3,2)^2),RotMatrix(3,3));
  gamma = atan2(RotMatrix(2,3),-RotMatrix(1,3));
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

Angles = [alpha, beta, gamma];

% Refinement by least-squares fitting
%-------------------------------------------------------------------------------
if ~skipFitting
  fitOptions = optimset('TolX',0.0001);
  RotErrorFunction = @(ang) norm(erot(ang)-RotMatrix);
  Angles = fminsearch(RotErrorFunction,Angles,fitOptions);
end


switch nargout
  case 0
    varargout = {Angles};
  case 1
    varargout = {Angles};
  case 3
    varargout = {Angles(1),Angles(2),Angles(3)};
  otherwise
    error('Wrong number of output arguments.')
end
