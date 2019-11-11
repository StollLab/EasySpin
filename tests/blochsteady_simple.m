function [err,data] = test(options,olddata)

g = gfree;       % g value
T1 = 20;         % longitudinal relaxation time, us
T2 = 1;          % transverse relaxation time, us
deltaB0 = 0.1;   % field offset, in mT
B1 = 0.02;       % microwave field, in mT
ModAmp = 0.5;    % peak-to-peak field modulation amplitude, in mT
ModFreq = 50;    % field modulation frequency, in kHz

[t,My] = blochsteady(g,T1,T2,deltaB0,B1,ModAmp,ModFreq);

data.My = My;

if ~isempty(olddata)
  ok = areequal(data.My,olddata.My,1e-6,'rel');
  err = ~ok;
else
  err = [];
end

if options.Display
  plot(t,data.My,t,olddata.My);
  legend('now','previous');
  legend boxoff
end
