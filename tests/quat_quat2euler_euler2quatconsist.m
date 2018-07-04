function [err,data] = test(opt,olddata)
% Check that quat2euler and euler2quat are consistent

N = 5;
M = 1;

% Euler->quaternion->Euler

% Generate some uniformly random Euler angles
alphaIn = pi*(2*rand(1,N,M)-1);
betaIn = pi*rand(1,N,M);
gammaIn = pi*(2*rand(1,N,M)-1);

omegaIn = [alphaIn; betaIn; gammaIn];

% Convert to quaternion
q = euler2quat(omegaIn);
% q = euler2quat(alphaIn, betaIn, gammaIn);

% Convert back to Euler angles
[alphaOut, betaOut, gammaOut] = quat2euler(q);
omegaOut = [alphaOut; betaOut; gammaOut];

% diffEqE = [alphaOut-alphaIn; betaOut-betaIn; gammaOut-gammaIn];
diffEqE = omegaIn - omegaOut;

% q = 2*rand(4,N,M);
% q = q./sqrt(sum(q.*q,1));

u1 = rand(1,N,M);
u2 = rand(1,N,M);
u3 = rand(1,N,M);

q0 = sqrt(1-u1).*sin(2*pi*u2);
q1 = sqrt(1-u1).*cos(2*pi*u2);
q2 =   sqrt(u1).*sin(2*pi*u3);
q3 =   sqrt(u1).*cos(2*pi*u3);

q = [q0; q1; q2; q3];

Index = cell(1, ndims(q));
Index(:) = {':'};

idx = q(1,Index{2:end}) < 0;

if any(idx(:))
  q(:,idx) = -q(:,idx);
end

[alpha, beta, gamma] = quat2euler(q);
qp = euler2quat(alpha, beta, gamma);

diffqEq = q - qp;

% Check for consistency
if any(abs(diffEqE(:))>1e-10)||any(abs(diffqEq(:))>1e-10)
  err = 1;
else  
  err = 0;
end

data = [];

end
