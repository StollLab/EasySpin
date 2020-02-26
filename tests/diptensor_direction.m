function ok = diptensor_direction()

% Test whether the directionality of the dipolar tensor is calculated correctly.

rvec = [1; -2; 3]; % nm

% Calculate using dipvector
T = diptensor('e','e',rvec); % MHz

% Calculate explicitly
rvec = rvec*1e-9; % nm -> m
r = norm(rvec);
n = rvec/r;
N = 3*(n*n')-eye(3);
pre = (mu0/4/pi)*r^-3*(bmagn*gfree)^2; % J
T0 = -pre*N; % J
T0 = T0/planck/1e6; % J -> MHz

ok = areequal(T,T0,1e-10,'rel');
