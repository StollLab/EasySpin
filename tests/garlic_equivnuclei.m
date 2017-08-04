function [err,data] = test(opt,olddata)

%=======================================================
% Make sure that '1H' and n=3 is the same as '1H,1H,1H'
%=======================================================
Exp.mwFreq = 9.669;
Exp.CenterSweep = [345 5];
Exp.Harmonic = 0;
Exp.nPoints = 10e3;
Opt.Method = 'perturb2';

for nEquivNuclei = 1:10
  Sys1.n = nEquivNuclei;
  Sys1.lw = [0 0.02];
  Sys1.Nucs = '1H';
  Sys1.A = 10;
  clear Sys2;
  Sys2.lw = Sys1.lw;
  for k = 1:nEquivNuclei
    Sys2 = nucspinadd(Sys2,'1H',Sys1.A);
  end
  Spc1 = garlic(Sys1,Exp,Opt);
  Spc2 = garlic(Sys2,Exp,Opt);
  area1(nEquivNuclei) = sum(Spc1);
  area2(nEquivNuclei) = sum(Spc2);
end

data = [];
err = isequal(area1,area2,1e-4);
