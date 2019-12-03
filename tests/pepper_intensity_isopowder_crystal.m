function [err,data] = test(opt,olddata)

% Compare integrated intensity of a powder spectrum and crystal
% spectrum for an isotropic-g spin system

Sys.g = [2];
Sys.lwpp = 0.5;

Exp.mwFreq = 9.6;
Exp.Range = [341 345];
Exp.Harmonic = 0;

Exp.CrystalOrientation = [];
[x,y1] = pepper(Sys,Exp);
Int1 = sum(y1)*(x(2)-x(1));

Exp.CrystalOrientation = rand(1,3)*2*pi;
[x,y2] = pepper(Sys,Exp);
Int2 = sum(y2)*(x(2)-x(1));

data = [];

err = ~areequal(Int1,Int2,1e-4,'abs');
