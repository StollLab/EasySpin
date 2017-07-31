function [err,data] = test(opt,olddata)
% Check that quat2euler and euler2quat are consistent

N = 4;

% Euler->quaternion->Euler

% Generate some uniformly random Euler angles
alphaIn = pi*(2*rand(1,N)-1);
betaIn = pi*rand(1,N);
gammaIn = pi*(2*rand(1,N)-1);

% Convert to quaternion
q = euler2quat(alphaIn, betaIn, gammaIn);

% Convert back to Euler angles
[alphaOut, betaOut, gammaOut] = quat2euler(q);

diffEqE = [alphaOut-alphaIn; betaOut-betaIn; gammaOut-gammaIn];

qp = euler2quat(alphaOut, betaOut, gammaOut);

diffEqEq = q - qp;

% Check for consistency
if any(abs(diffEqE(:))>1e-13)||any(abs(diffEqEq(:))>1e-13)
  err = 1;
else  
  err = 0;
end

data = [];

end
