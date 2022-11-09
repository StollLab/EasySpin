function [ok,data] = test(opt,olddata)

% Simulate rapid-scan signal
g = 2;
T1 = 4;         % µs
T2 = 0.5;       % µs
B1 = 1e-3;      % mT
deltaB0 = 0;    % mT
modAmp = 1;     % mT
modFreq = 50;   % kHz

Opt.nPoints = 1000;

[t,Mx,My,Mz] = blochsteady(g,T1,T2,deltaB0,B1,modAmp,modFreq,Opt);
M = -Mx + 1i*My;

[B,spc] = rapidscan2spc(M,modAmp,modFreq,2);

data.spc = spc;

if opt.Display
  plot(B,real(spc),B,real(olddata.spc));
  legend('current','previous');
end

if ~isempty(olddata)
  ok = areequal(data.spc,olddata.spc,1e-6,'rel');
else
  ok = [];
end
