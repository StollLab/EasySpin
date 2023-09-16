function [ok,data] = test(opt,olddata)

clear Sys Opt
% Spin System
Sys.S = 1/2;
Sys.g = [2.00906 2.0062 2.0023];
Sys.A = [11.5 11.5 95];
Sys.Nucs = '14N';
Sys.lwpp = 10;

Exp.Field = 1240;
Exp.mwFreq = 34.78;

Pulse1.Type = 'quartersin/linear';
Pulse1.trise = 0.030;
Pulse1.tp = 0.200;
Pulse1.Flip = pi/2;
Pulse1.Frequency = [-300 300];

Pulse2.Type = 'quartersin/linear';
Pulse2.trise = 0.030;
Pulse2.tp = 0.100;
Pulse2.Flip = pi;
Pulse2.Frequency = [-300 300];

Exp.Sequence = {Pulse1 0.25 Pulse2 0.25}; 
Exp.DetWindow = [-0.05 0.05] + Pulse2.tp;
Exp.DetPhase = 0;

Opt.GridSize = 7;
% Opt.SimulationMode = 'thyme';
Opt.SimFreq = 15;

[x, y] = saffron(Sys,Exp,Opt);

data.x = x;
data.y = y;

if opt.Display
  if ~isempty(olddata)
    p1 = subplot(3,1,[1 2]);
    plot(x,real(y),x,real(olddata.y));
    axis tight
    legend('new','old');
    title(mfilename);
    p2 = subplot(3,1,3);
    plot(x,real(olddata.y-y));
    axis tight
    xlabel('time (µs)');
    title('old - new')
    linkaxes([p1,p2],'x')
  end
end

if ~isempty(olddata)
  ok = areequal(y,olddata.y,1e-3,'abs');
else
  ok = [];
end
