function [err,data] = test(opt,olddata)

% Correct behavior of Exp.CenterField and Exp.Range
%-----------------------------------------------------

Sys = struct('g',2,'Nucs','1H','lw',[0,0.01],'A',50);

mw = 9.7;

Exp = struct('mwFreq',mw);
[x0,y0] = garlic(Sys,Exp);

Exp.CenterSweep = [mean(x0) max(x0)-min(x0)];
[x1,y1] = garlic(Sys,Exp);

Exp = struct('mwFreq',mw);
Exp.Range = [min(x0) max(x0)];
[x2,y2] = garlic(Sys,Exp);

err = ~areequal(x0,x1,1e-10,'rel') | ~areequal(x0,x2,1e-10,'rel');

data = [];
