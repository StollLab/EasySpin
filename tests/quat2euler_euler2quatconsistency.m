function [err,data] = test(opt,olddata)
% Check that rotating vectors by quaternions is equivalent to using
% rotation matrices

% Generate some uniformly random Euler angles
alphaIn = 2*pi*rand();
betaIn = pi*rand();
gammaIn = 2*pi*rand();

% Convert to quaternion
q = euler2quat(alphaIn, betaIn, gammaIn);

% Convert back to Euler angles
[alphaOut, betaOut, gammaOut] = quat2euler(q);

diff = [alphaOut-alphaIn, betaOut-betaIn, gammaOut-gammaIn];

% Check for consistency
if any(diff(:)>1e-14)
  err = 1;
else  
  err = 0;
end

data = [];

end
