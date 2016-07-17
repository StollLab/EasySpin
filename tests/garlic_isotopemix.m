function [err,data] = test(opt,olddata)

%=======================================================
% Isotope mixtures: spectrum should integrate to 1
%=======================================================
Sys.Nucs = 'B';
Sys.A = 10;
Exp.mwFreq = 9.668;
Exp.CenterSweep = [345 2];
Exp.Harmonic = 0;
[x,y]=garlic(Sys,Exp);

err = isequal(sum(y),1,1e-6);
data = [];
