function [err,data] = test(opt,olddata)

%=======================================================
% garlic should give the same integral intensity as pepper
%=======================================================
Sys.g = 2.1;
Sys.Nucs = '1H';
Sys.A = 100;
Sys.lwpp = 1;

Exp.mwFreq = 9.5;
Exp.Range = [300 360];
Exp.Harmonic = 0;

[x,y1] = pepper(Sys,Exp);
[x,y2] = garlic(Sys,Exp);
dx = x(2)-x(1);

integral_garlic = sum(y1)*dx;
integral_pepper = sum(y2)*dx;

err = abs(integral_garlic-integral_pepper)>0.001;

data = [];
