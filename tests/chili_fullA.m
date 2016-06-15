function [err,data] = test(opt,olddata)

%=============================================================
% Simple simulation with full g and A tensors (1  nucleus)
%=============================================================
g = [2.008 2.006 2.002];
A = [20 20 85];
Sys.Nucs = '15N';
Sys.tcorr = 4e-9;
Exp.mwFreq = 9.8;
Exp.Range = [344 354];

Sys.g = g;
Sys.A = A;
[x1,y1] = chili(Sys,Exp);

Sys.g = diag(g);
Sys.A = A;
[x2,y2] = chili(Sys,Exp);

Sys.g = g;
Sys.A = diag(A);
[x3,y3] = chili(Sys,Exp);

Sys.g = diag(g);
Sys.A = diag(A);
[x4,y4] = chili(Sys,Exp);

data = [];
thr = 1e-3;
err = ~areequal(y1,y2,thr) || ~areequal(y1,y3,thr) || ~areequal(y1,y4,thr);
