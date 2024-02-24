function ok = diptensor_ee()

% Test whether diptensor gives the correct result for two electrons

rvec = [0;0;1];  % nm
T = diptensor(gfree,gfree,rvec);
Tperp = T(1,1);

rvec = rvec*1e-9;  % nm -> m
r = norm(rvec);
Tperp_ref = (mu0/4/pi)*r^-3*(bmagn*gfree)^2;  % J
Tperp_ref = Tperp_ref/planck/1e6;  % J -> MHz

ok = areequal(Tperp,Tperp_ref,1e-10,'rel');
