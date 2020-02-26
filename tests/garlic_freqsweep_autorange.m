function ok = test()

% garlic: Test auto ranging for frequency sweeps

Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = 100;
Sys.lw = 10;

Exp.Field = 340;

[nu,y] = garlic(Sys,Exp);

ok = min(nu)<max(nu);
