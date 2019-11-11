function [err,data] = test(options,olddata)

g = gfree;       % g value
T1 = 20;         % longitudinal relaxation time, us
T2 = 1;          % transverse relaxation time, us
deltaB0 = 0.1;   % field offset, in mT
B1 = 0.02;       % microwave field, in mT
ModAmp = 0.5;    % peak-to-peak field modulation amplitude, in mT
ModFreq = 50;    % field modulation frequency, in kHz

[t,Mx,My,Mz] = blochsteady(g,T1,T2,deltaB0,B1,ModAmp,ModFreq);

data.Mx = Mx;
data.My = My;
data.Mz = Mz;

if ~isempty(olddata)
  thr = 1e-6;
  ok = areequal(data.Mx,olddata.Mx,thr,'rel') && ...
       areequal(data.My,olddata.My,thr,'rel') && ...
       areequal(data.Mz,olddata.Mz,thr,'rel');
  err = ~ok;
else
  err = [];
end

if options.Display
  subplot(3,1,1)
  plot(t,data.Mx,'.',t,olddata.Mx);
  legend('now','previous');
  legend boxoff
  ylabel('M_x');
  grid on
  subplot(3,1,2)
  plot(t,data.My,'.',t,olddata.My);
  ylabel('M_y');
  grid on
  subplot(3,1,3)
  plot(t,data.Mz,'.',t,olddata.Mz);
  ylabel('M_z');
  xlabel('time (us)');
  grid on
end
