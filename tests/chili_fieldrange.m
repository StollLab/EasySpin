function [err,data] = test(opt,olddata)

% Correct behaviour of Exp.CenterField and Exp.Range
%-----------------------------------------------------
Sys = struct('g',[2.008 2.0061 2.0027],'Nucs','14N','A',[16 16 86]);
Sys.lw = 0.1; Sys.tcorr = 1e-7;

mw = 9.7;

Exp = struct('mwFreq',mw);
[x0,y0] = chili(Sys,Exp);

Exp.CenterSweep = [mean(x0) max(x0)-min(x0)];
[x1,y1] = chili(Sys,Exp);

Exp = struct('mwFreq',mw);
Exp.Range = [min(x0) max(x0)];
[x2,y2] = chili(Sys,Exp);

err = ~areequal(x0,x1) | ~areequal(x0,x2);

data = [];
