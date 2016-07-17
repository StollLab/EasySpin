function [err,data] = test(opt,olddata)

%======================================================
% 3 electron spins, check full vs simple
%======================================================
Sys.S = [1/2 1/2 1/2];
ee1 = [1 1 1];
ee2 = [2 2 2];
ee3 = [5 5 5];

Sys1.S = [1/2 1/2 1/2];
Sys1.ee = [ee1; ee2; ee3];
Sys2.S = Sys1.S;
Sys2.ee = [diag(ee1);diag(ee2);diag(ee3)];

H1_13 = eeint(Sys1,[1 3]);
H2_13 = eeint(Sys2,[1 3]);
H1_12 = eeint(Sys1,[1 2]);
H2_12 = eeint(Sys2,[1 2]);
H1_23 = eeint(Sys1,[2 3]);
H2_23 = eeint(Sys2,[2 3]);

ok = areequal(H1_12,H2_12) & areequal(H1_13,H2_13) & areequal(H1_23,H2_23);
err = ~ok;

data = [];
