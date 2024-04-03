function [ok,data] = test(opt,olddata)

% Spin system
Sys.S = 1/2;
Sys.g = [2.00906 2.0062 2.0023];
Sys.A = [11.5 11.5 95];  % MHz
Sys.Nucs = '14N';
Sys.lwpp = 10;  % mT

% Experiment (Q band) 
Exp.mwFreq = 34.78;  % GHz
Exp.Field = 1240;  % mT

% Define pulses
Pulse1.Type = 'quartersin/linear';
Pulse1.trise = 0.030;  % µs
Pulse1.tp = 0.200;  % µs
Pulse1.Flip = pi/2;  % rad
Pulse1.Frequency = [-300 300]; % MHz

Pulse2.Type = 'quartersin/linear';
Pulse2.trise = 0.030;  % µs
Pulse2.tp = 0.100;  % µs
Pulse2.Flip = pi;  % rad
Pulse2.Frequency = [-300 300];  % MHz

% Define pulse sequence
tau = 0.25;  % µs
Exp.Sequence = {Pulse1 tau Pulse2 tau}; 

% Define detection
Exp.DetWindow = [-0.1 0.1] + Pulse2.tp;
Exp.DetPhase = 0;

% Options
Opt.GridSize = 7;
Opt.SimulationMode = 'thyme';

[x, y] = saffron(Sys,Exp,Opt);

data.x = x;
data.y = y;

if opt.Display
  if ~isempty(olddata)
    Pulse1 = subplot(3,1,[1 2]);
    plot(x,real(y),x,real(olddata.y));
    axis tight
    legend('new','old');
    title(mfilename);
    Pulse2 = subplot(3,1,3);
    plot(x,real(olddata.y-y));
    axis tight
    xlabel('time (µs)');
    title('old - new')
    linkaxes([Pulse1,Pulse2],'x')
  end
end

if ~isempty(olddata)
  ok = areequal(y,olddata.y,1e-3,'abs');
else
  ok = [];
end

