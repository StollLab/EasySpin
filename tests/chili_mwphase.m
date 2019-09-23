function [err,data] = test(opt,olddata)

%=======================================================
% Assert that mw phasing is working in chili
%=======================================================
Sys.g = [2.01 2.00];
Sys.tcorr = 1e-9;
Exp.mwFreq = 9.8;
Exp.Harmonic = 0;
Exp.Range = [345 353];

Exp.mwPhase = 90*pi/180;
[x,y] = chili(Sys,Exp);

data.y = y;

if ~isempty(olddata)
  if opt.Display
    plot(x,data.y/max(data.y),x,olddata.y/max(olddata.y));
    legend('current','old');
    legend boxoff
    title(sprintf('MAD = %g',max(abs(data.y/max(data.y)-olddata.y/max(olddata.y)))));
  end
  ok = areequal(data.y/max(data.y),olddata.y/max(olddata.y),5e-4,'abs');
  err = ~ok;
else
  err = [];
end
