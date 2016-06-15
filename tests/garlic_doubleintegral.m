function [err,data] = test(opt,olddata)

% make sure doulbe integral does not depend on number of nuclei

Sys = struct('g',2,'lw',[0,0.02]);
Exp = struct('mwFreq',9.7,'Harmonic',0);
Opt = struct('Verbosity',opt.Verbosity);

Exp.Range = [337 355];
Exp.nPoints = 1e4;
for iNuc=1:6
  Sys = nucspinadd(Sys,'1H',rand*20);
  y = garlic(Sys,Exp,Opt);
  Integral(iNuc) = sum(y);
end

err = any(abs(Integral/Integral(1)-1)>1e-3);

data = [];
