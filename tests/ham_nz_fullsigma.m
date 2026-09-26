function ok = test()

% Assert that a full (non-symmetric) shielding tensor is used as
% mu_i = pre*sum_k sigma(i,k)*I_k

Sys.Nucs = '1H';
Sys.A = 0;  % needed to pass validation
Sys.sigma = [1 0.2 0.5; 0.1 2 0.3; 0.6 0.4 3];

[mux,muy,muz] = ham_nz(Sys);

[Ix,Iy,Iz] = sop(Sys,'x2','y2','z2');
pre = +nucgval(Sys.Nucs)*nmagn/planck/1e9;  % MHz/mT
mux0 = pre*(1*Ix + 0.2*Iy + 0.5*Iz);
muy0 = pre*(0.1*Ix + 2*Iy + 0.3*Iz);
muz0 = pre*(0.6*Ix + 0.4*Iy + 3*Iz);

thr = 1e-10;
ok = areequal(mux,mux0,thr,'rel') && areequal(muy,muy0,thr,'rel') && areequal(muz,muz0,thr,'rel');
