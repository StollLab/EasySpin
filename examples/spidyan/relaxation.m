clear Exp Sys Opt Pulse
% Default Spin System
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500; % GHz
Sys.T1 = 1; % us
% Sys.T2 = 0.5; % us
Sys.eqState = -sop(Sys.S,'z');

% Pulse Definitions
Rectangular.Type = 'rectangular';

% A default Experiment/Sequence
Exp.Field = 1240; % mT 
Exp.TimeStep = 0.0001; % us
Exp.mwFreq = 33.5; % GHz
Exp.DetEvents = 1;

Exp.t = [0.2 5]; % us
Exp.Pulses = {Rectangular 0};
Exp.Frequency = 0;
Exp.Flip = pi/2;

% Options
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32; % GHz
Opt.Relaxation = [0 1];

[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

% plotting
figure(1)
clf
plot(TimeAxis*1000,real(Signal))
xlabel('t [ns]')
ylabel('<S_i>')
legend(Opt.DetOperator)