function [err,data] = test(opt,olddata)

%=======================================================
% Fast motional regime: check whether it runs without nuclei
%=======================================================

Spins.g = [2.01 2.00 1.99];
Spins.tcorr = 2e-9;
Exp = struct('mwFreq',9.7);
Exp.CenterSweep = [346.5 10];
[x,y] = garlic(Spins,Exp);

data = [];
err = 0;
