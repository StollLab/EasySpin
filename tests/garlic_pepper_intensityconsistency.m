function ok = test()

% garlic: Test auto ranging for frequency sweeps

Sys.g = 2;
Sys.A = 300;
Sys.Nucs = '1H';
Sys.lwpp = 0.2;

Exp.mwFreq = 9.5;
Exp.Range = [320 360];
Exp.Harmonic = 0;

[x,yp] = pepper(Sys,Exp);
[x,yg] = garlic(Sys,Exp);

maxyp = max(yp);

ok = areequal(max(yg)/maxyp,1,1e-3,'abs');
