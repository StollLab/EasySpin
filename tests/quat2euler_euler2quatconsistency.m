function [err,data] = test(opt,olddata)
% Check that quat2euler and euler2quat are consistent

% Generate some uniformly random Euler angles
alphaIn = pi*(2*rand()-1);
betaIn = pi*rand();
gammaIn = pi*(2*rand()-1);

% Convert to quaternion
q = euler2quat(alphaIn, betaIn, gammaIn);

% Convert back to Euler angles
[alphaOut, betaOut, gammaOut] = quat2euler(q);

diff = [alphaOut-alphaIn, betaOut-betaIn, gammaOut-gammaIn];

% Check for consistency
if any(abs(diff(:))>1e-14)
  err = 1;
else  
  err = 0;
end

data = [];

end
