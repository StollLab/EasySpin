function [err,data] = test(opt,olddata)
% Check that quat2euler and euler2quat are consistent

N = 5;
M = 1;

% Euler->quaternion->Euler

% Generate some uniformly random Euler angles
alphaIn = 2*pi*rand(1,N,M);
betaIn = pi*rand(1,N,M);
gammaIn = 2*pi*rand(1,N,M);
% alphaIn = linspace(0,2*pi,N);
% betaIn = linspace(0,pi,N);
% gammaIn = linspace(0,2*pi,N);

[Alpha,Beta,Gamma] = ndgrid(alphaIn,betaIn,gammaIn);
alphaIn = Alpha(:).';
betaIn = Beta(:).';
gammaIn = Gamma(:).';

omegaIn = [alphaIn; betaIn; gammaIn];

% Convert to quaternion
q = euler2quat(omegaIn);

% Convert back to Euler angles
[alphaOut, betaOut, gammaOut] = quat2euler(q);
omegaOut = [alphaOut; betaOut; gammaOut];

diffEqE = omegaIn - omegaOut;
% omegaInProb = omegaIn(abs(diffEqE)>1e-3);
% omegaInProb = omegaIn(abs(diffEqE)>1e-3);
% omegaIn(abs(diffEqE)>1e-3)/pi*180
% diffEqE(abs(diffEqE)>1e-3)/pi*180

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
err = any(abs(diffEqE(:))>1e-10)||any(abs(diffqEq(:))>1e-10);

data = [];

end
