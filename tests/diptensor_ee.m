function err = diptensor_ee(opt)

% Test whether diptensor gives the correct result for two electrons

rvec = [0;0;1]; % nm
T = diptensor('e','e',rvec);
Tperp = T(1,1);

rvec = rvec*1e-9; % nm -> m
r = norm(rvec);
Tperp0 = (mu0/4/pi)*r^-3*(bmagn*gfree)^2; % J
Tperp0 = Tperp0/planck/1e6; % J -> MHz

err = ~areequal(Tperp,Tperp0,1e-10,'rel');
