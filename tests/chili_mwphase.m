function [err,data] = test(opt,olddata)

%=======================================================
% Assert that mw phasing is working in chili
%=======================================================
Sys.g = [2.01 2.00];
Sys.tcorr = 1e-9;
Exp.mwFreq = 9.8;
Exp.Harmonic = 0;
Exp.Range = [345 353];

Exp.mwPhase = 30*pi/180;
[x,y] = chili(Sys,Exp);

data.y = y;

if ~isempty(olddata)
  if opt.Display
    plot(x,data.y,x,olddata.y);
    legend('current','old');
    legend boxoff
  end
  ok = areequal(data.y,olddata.y,1e-4);
  err = ~ok;
else
  err = [];
end
