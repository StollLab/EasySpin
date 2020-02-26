function [ok,data] = rs2spc_basic(opt,olddata)

% Simulate rapid-scan signal
g = 2;
T1 = 4;    % microseconds
T2 = 0.5;  % microseconds
B1 = 1e-3;      % mT
deltaB0 = 0;    % mT
ModAmp_mT = 1;   % mT
ModFreq = 50;    % kHz  

[t,Mx,My,Mz] = blochsteady(g,T1,T2,deltaB0,B1,ModAmp_mT,ModFreq);
M = -Mx + 1i*My;


[B,spc] = rapidscan2spc(M,ModAmp_mT,ModFreq,2);

data.spc = spc;

if ~isempty(olddata)
  ok = areequal(data.spc,olddata.spc,1e-6,'rel');
else
  ok = [];
end
