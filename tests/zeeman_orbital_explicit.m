function ok = test()

Sys.S = 1/2;
Sys.L = 1;
Sys.g = 1e-7;
Sys.soc = 1e-7;

gL = 123;
Sys.orf = gL;

H = zeeman(Sys,[0;0;1]);
E = diag(H)/bmagn*planck*1e9;

Eref = [1 0 -1 1 0 -1].'*gL;
ok = areequal(E,Eref,1e-9,'rel');
