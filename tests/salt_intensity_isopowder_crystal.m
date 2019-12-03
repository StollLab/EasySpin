function [err,data] = test(opt,olddata)

% Compare integrated intensity of a powder spectrum and crystal
% spectrum for an isotropic-g spin system

clear Sys Exp
Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = 2;
Sys.lwEndor = 0.5;

Exp.Range = [11 18];
Exp.Field = 350;

Exp.CrystalOrientation = [];
[x,y1] = salt(Sys,Exp);
Int1 = sum(y1)*(x(2)-x(1));

Exp.CrystalOrientation = rand(1,3)*2*pi;
[x,y2] = salt(Sys,Exp);
Int2 = sum(y2)*(x(2)-x(1));

data = [];

err = ~areequal(Int1,Int2,1e-8,'abs');
