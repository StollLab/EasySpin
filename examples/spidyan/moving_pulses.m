clear Exp Sys Opt Pulse
% This script explains how to move a pulse around using Exp.Dim in
% combination with the option 'Position'. The results are shown as
% trajectory of a single spin wich is on resoance with the pulses

% Spin System
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500; % GHz

Pulse.Type = 'rectangular';

%% Creates a pi/2 - tau1 - pi - tau2 - pi - tau3 sequence
Exp.t = [0.025 0.1 0.05 0.9 0.05 1]; % us
Exp.Pulses = {Pulse 0 Pulse 0 Pulse 0};
Exp.Field = 1240; % mT
Exp.TimeStep = 0.0001; % us
Exp.Frequency = 0; % GHz
Exp.Flip = [pi/2 pi pi];
Exp.mwFreq = 33.5; % GHz
Exp.DetEvents = 1;

% Move the center pulse to the back in 100 ns steps
Exp.nPoints = 9;
Exp.Dim1 = {'p2.Position',0.1};

% Options
Opt.DetOperator = {'+1'};
Opt.FreqTranslation = -33.5; % GHz
Opt.FrameShift = 32; % GHz

% Run simulation
[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

% Plotting
figure(1)
clf
hold on
for i = 1:size(Signal,1)
  plot(TimeAxis*1000,real(squeeze(Signal(i,:,:))));
end

xlabel('t [ns]')
ylabel('<S_z>')
axis tight
ylim([-1 1])