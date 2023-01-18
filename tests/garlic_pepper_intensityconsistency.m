function ok = test()

% garlic: Test auto ranging for frequency sweeps

Sys.g = 2;
Sys.A = 300;
Sys.Nucs = '1H';
Sys.lwpp = 1.5;

Exp.mwFreq = 9.5;
Exp.Range = [320 360];
Exp.Harmonic = 0;
Exp.nPoints = 1e4;

[x,yp] = pepper(Sys,Exp);
[x,yg] = garlic(Sys,Exp);

maxyp = max(yp);

ok = areequal(max(yg)/maxyp,1,1e-3,'rel');
