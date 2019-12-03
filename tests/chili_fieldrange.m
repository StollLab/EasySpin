function [err,data] = test(opt,olddata)

% Correct behavior of Exp.CenterField and Exp.Range
%-----------------------------------------------------
Sys.g = [2.008 2.0061 2.0027];
Sys.Nucs = '14N';
Sys.A = [16 16 86];
Sys.lw = 0.1;
Sys.tcorr = 1e-7;

mw = 9.7;

Exp = struct('mwFreq',mw);
[x0,y0] = chili(Sys,Exp);

Exp.CenterSweep = [mean(x0) max(x0)-min(x0)];
[x1,y1] = chili(Sys,Exp);

Exp = struct('mwFreq',mw);
Exp.Range = [min(x0) max(x0)];
[x2,y2] = chili(Sys,Exp);

thr = 1e-12;
err = ~areequal(x0,x1,thr,'rel') || ~areequal(x0,x2,thr,'rel');

data = [];
