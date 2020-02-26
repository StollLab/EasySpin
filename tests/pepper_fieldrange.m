function ok = test()

% Correct behaviour of Exp.CenterField and Exp.Range
%-----------------------------------------------------

Sys = struct('g',[2 2.1 2.2],'lw',1);

mw = 35;

% Automatic range
Exp = struct('mwFreq',mw);
[x0,y0] = pepper(Sys,Exp);

% Center and sweep
Exp.CenterSweep = [mean(x0) max(x0)-min(x0)];
[x1,y1] = pepper(Sys,Exp);

% Min and Max
Exp = struct('mwFreq',mw);
Exp.Range = [min(x0) max(x0)];
[x2,y2] = pepper(Sys,Exp);

ok = areequal(x0,x1,1e-10,'rel') && areequal(x0,x2,1e-10,'rel');
