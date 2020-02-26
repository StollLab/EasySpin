function ok = test()

% Test whether diptensor gives correct results for nuclus-nucleus couplings

rvec = [0;0;1]*0.2; % nm
r = norm(rvec)*1e-9; % nm -> m

Nuc1 = '1H';
Nuc2 = '13C';
T = diptensor(Nuc1,Nuc2,rvec);
Tperp = T(1,1);

gn1 = nucgval(Nuc1);
gn2 = nucgval(Nuc2);
Tperp0 = (mu0/4/pi)*r^-3*nmagn^2*gn1*gn2; % J
Tperp0 = Tperp0/planck/1e6; % J -> MHz

ok = areequal(Tperp,Tperp0,1e-10,'rel');
