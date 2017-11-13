clear Exp Sys Opt Pulse

% Spin System
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500; % GHz

% Pulse Definitions
Pulse.Type = 'rectangular';

% Define the phase cycle with phase and detection phase 
PhaseCycle = [0   1; 
              pi -1];

% Experiment/Sequence
Exp.t = [0.025 0.4]; % us
Exp.Frequency = 0; % GHz
Exp.Pulses = {Pulse 0 Pulse};
Exp.Flip = pi;
Exp.Field = 1240; % mT
Exp.TimeStep = 0.0001; % us
Exp.mwFreq = 33.5; % GHz
Exp.DetEvents = 1;
% Add the phasecycle to the experiment structure
Exp.PhaseCycle = {PhaseCycle};

% Options
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32; % GHz

[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

% Plotting 
figure(1)
clf
plot(TimeAxis,Signal)
ylim([-1 1])
xlabel('t [ns]')
ylabel('<S_z>')