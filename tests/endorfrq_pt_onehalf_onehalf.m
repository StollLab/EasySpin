function [err,data] = test(opt,olddata)

%-------------------------------------------------------------
%S=1/2, I=1/2
%-------------------------------------------------------------
clear
Sys.Nucs='13C';
Exp.Field=352.3;
Exp.mwFreq=9.5;
Exp.CrystalOrientation = [0 pi/4 0];

% weak coupling
Sys.A = [5 9];
a=endorfrq(Sys,Exp);
b=endorfrq_perturb(Sys,Exp);
err1 = any(abs(sort(a)-sort(b))>0.01);

% strong coupling
Sys.A = [30 40];
a=endorfrq(Sys,Exp);
b=endorfrq_perturb(Sys,Exp);
err2 = any(abs(sort(a)-sort(b))>0.1);

err = err1 || err2;
data = [];
