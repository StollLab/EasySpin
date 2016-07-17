function [err,data] = test(opt,olddata)

% Make sure pepper calculates stick absorption spectrum when no
% line broadening is specified.

Sys.Nucs = '14N';
Sys.A = [20 20 100];

Exp.mwFreq = 0.250;
Exp.Range = [1 15];

[x,y1] = pepper(Sys,Exp);
Exp.Harmonic = 0;
[x,y2] = pepper(Sys,Exp);

if opt.Display
  plot(x,y1,x,y2,'r');
  legend('auto Harmonic','Harmonic=0');
  legend boxoff
end

err = ~areequal(y1,y2,1e-6*max(y1));
data = [];
