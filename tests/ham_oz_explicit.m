function ok = test()

gL = 5/3;

Sys.S = 1/2;
Sys.L = 1;
Sys.g = 2;
Sys.gL = gL;
Sys.soc = 1000;
H = ham_oz(Sys,[0;0;1]);
E = diag(H)/bmagn*planck*1e9;

Eref = [1 0 -1 1 0 -1].'*gL;

ok = areequal(E,Eref,1e-9,'rel');
