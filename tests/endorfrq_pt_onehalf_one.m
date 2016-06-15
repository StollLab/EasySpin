function [err,data] = test(opt,olddata)

%-------------------------------------------------------------
%S=1/2, I=1/2
%-------------------------------------------------------------
clear
Sys.Nucs='2H';
Sys.A = [5 9];
Exp.Field=3500;
Exp.CrystalOrientation = [0 pi/4 0];
Exp.mwFreq = 95;
a=endorfrq(Sys,Exp);
b=endorfrq_perturb(Sys,Exp);
err = any(abs([sort(a)-sort(b)])>0.01);

Sys.Q = 1;
a=endorfrq(Sys,Exp);
b=endorfrq_perturb(Sys,Exp);
err = err || any(abs([sort(a)-sort(b)])>0.01);

data = [];
