function ok = test(opt)

Sys.S = 1;
Sys.g = [2.01 2.005 2.002];
Sys.D = 500;
Sys.tcorr = 10e-9;
Exp.Field = 339.4;
Exp.mwRange = [8.4 10.5];
Exp.Harmonic = 1;

% without additional broadening
[B,spc0] = chili(Sys,Exp);

% with additional broadening in Sys.lwpp
Sys.lwpp = 2;
[B,spc1] = chili(Sys,Exp);

spc0 = spc0/max(spc0);
spc1 = spc1/max(spc1);

ok = areequal(spc0,spc1,1e-2,'abs');

if opt.Display
  plot(B,spc0,'.',B,spc1);
  legend('without lwpp','with lwpp');
  axis tight
end
