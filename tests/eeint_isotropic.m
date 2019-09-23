function [err,data] = test(opt,olddata)

%======================================================
% expansion of isotropic ee
%======================================================
ee = [100 121 37];

Sys1.S = [1/2 1/2 1/2];
Sys1.ee = ee;

Sys2.S = Sys1.S;
Sys2.ee = ee(:)*[1 1 1];

H1 = eeint(Sys1);
H2 = eeint(Sys2);

ok = areequal(H1,H2,1e-10,'abs');
err = ~ok;

data = [];
