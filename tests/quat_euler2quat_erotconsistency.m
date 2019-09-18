function [err,data] = test(opt,olddata)

% consistency with erot()
%======================================================

alpha = pi*(2*rand-1);
beta = pi*rand;
gamma = pi*(2*rand-1);

vec = 2*rand(3,1)-1;
vec = vec/norm(vec);

% test rotation by quaternion
q = euler2quat(alpha, beta, gamma);
vecp = quatvecmult(q, vec);

% test rotation by matrix
R = quat2rotmat(q);

R_erot = erot(alpha,beta,gamma);

diffq = vecp - R_erot*vec;

diffR = R - R_erot;

err1 = any(abs(diffq(:))>1e-10);

err2 = any(abs(diffR(:))>1e-10);

err = err1 + err2;

data = [];
