function ok = diptensor_direction()

% Test whether the directionality of the dipolar tensor is calculated correctly.

rvec = [1; -2; 3];  % nm

% Calculate using dipvector
T = diptensor(gfree,gfree,rvec); % MHz

% Calculate explicitly
rvec = rvec*1e-9;  % nm -> m
r = norm(rvec);
n = rvec/r;
N = 3*(n*n')-eye(3);
pre = (mu0/4/pi)*r^-3*(bmagn*gfree)^2; % J
T_ref = -pre*N;  % J
T_ref = T_ref/planck/1e6;  % J -> MHz

ok = areequal(T,T_ref,1e-10,'rel');
