function [err,data] = test(opt,olddata)

%=======================================================
% Simple simulation with full g tensor
%=======================================================
g = [2 2.005 2.02];
Sys.g = g;
Sys.tcorr = 4e-9;
Exp.mwFreq = 9.8;
Exp.Range = [344 354];
[x1,y1] = chili(Sys,Exp);
Sys.g = diag(g);
[x2,y2] = chili(Sys,Exp);

data = [];
err = ~areequal(y1,y2,1e-2,'abs');
