function [err,data] = test(opt,olddata)

% consistency with erot()
%======================================================

alpha = 2*pi*(2*rand-1);
beta = pi*(2*rand-1);
gamma = 2*pi*(2*rand-1);

% Need to convert quaternion to rotation matrix
R = quat2rotmat(euler2quat(alpha,beta,gamma));

diff = R - erot(alpha,beta,gamma);

err = any(diff(:)>1e-10);

data = [];
