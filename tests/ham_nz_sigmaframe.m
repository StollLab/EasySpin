function ok = test()

% Assert that the shielding tensor frame is applied correctly, by comparing
% to a hand-computed tensor in the molecular frame

Sys.Nucs = '1H';
Sys.A = 0;  % needed to pass validation
Sys.sigma = [1 2 3];
Sys.sigmaFrame = [pi/4 0 0];  % +45 degrees about z

[mux,muy,muz] = ham_nz(Sys);

% Principal axis with value 1 is along (x+y)/sqrt(2) in the molecular frame,
% giving sigma = [1.5 -0.5 0; -0.5 1.5 0; 0 0 3] in the molecular frame
[Ix,Iy,Iz] = sop(Sys,'x2','y2','z2');
pre = +nucgval(Sys.Nucs)*nmagn/planck/1e9;  % MHz/mT
mux0 = pre*(1.5*Ix - 0.5*Iy);
muy0 = pre*(-0.5*Ix + 1.5*Iy);
muz0 = pre*3*Iz;

thr = 1e-10;
ok = areequal(mux,mux0,thr,'rel') && areequal(muy,muy0,thr,'rel') && areequal(muz,muz0,thr,'rel');
