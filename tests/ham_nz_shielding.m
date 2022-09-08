function ok = test()

Sys.Nucs = '1H';
Sys.sigma = [1 2 3];
Sys.sigmaFrame = [pi/3 pi/5 pi/7];
Sys.gnscale = 0.9;
Sys.gn = nucgval(Sys.Nucs);
Sys.A = 0; % needed to pass validation

% Calculate operators using ham_nz()
[Gx,Gy,Gz] = ham_nz(Sys);

% Calculate operators explicitly
[I{1:3}] = sop(Sys,'x2','y2','z2');

pre = -nmagn/(planck*1e9)*Sys.gn.*Sys.gnscale;

sigma = diag(Sys.sigma);  % eigenframe of sigma
R = erot(Sys.sigmaFrame);
sigma = R.'*sigma*R;  % transform to molecular frame

Gx0 = 0;
Gy0 = 0;
Gz0 = 0;
for k = 1:3
  Gx0 = Gx0 + pre*sigma(1,k)*I{k};
  Gy0 = Gy0 + pre*sigma(2,k)*I{k};
  Gz0 = Gz0 + pre*sigma(3,k)*I{k};
end

thr = 1e-10;
ok(1) = areequal(Gx,Gx0,thr,'rel');
ok(2) = areequal(Gy,Gy0,thr,'rel');
ok(3) = areequal(Gz,Gz0,thr,'rel');
