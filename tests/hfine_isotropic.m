function [err,data] = test(opt,olddata)

%======================================================
% expansion of isotropic hfc
% system with 1 electron and 3 nuclei
%======================================================
aiso = [100 121 37];

Sys1.Nucs = '1H,1H,1H';
Sys1.A = aiso;

Sys2.Nucs = Sys1.Nucs;
Sys2.A = aiso(:)*[1 1 1];

H1 = hfine(Sys1);
H2 = hfine(Sys2);

ok = areequal(H1,H2);
err = ~ok;

data = [];
