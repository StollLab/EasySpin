clear Exp Sys Opt Pulse90 Pulse180
% This script explains how to move a pulse around using Exp.Dim in
% combination with the option 'Position'. The results are shown as
% trajectory of a single spin wich is on resoance with the pulses

% Spin System
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500; % GHz

Pulse90.Type = 'rectangular';
Pulse90.Flip = pi/2;
Pulse90.tp = 0.025;

Pulse180.Type = 'rectangular';
Pulse180.Flip = pi;
Pulse180.tp = 0.05;

%% Creates a pi/2 - tau1 - pi - tau2 - pi - tau3 sequence
Exp.Sequence = {Pulse90 0.1 Pulse180 0.9 Pulse180 1};
Exp.Field = 1240; % mT
% Exp.TimeStep = 0.00001; % us

Exp.mwFreq = 33.5; % GHz
Exp.DetEvents = 1;

% Move the center pulse to the back in 100 ns steps
Exp.nPoints = 9;
Exp.Dim1 = {'p2.Position',0.1};

% Options
Opt.DetOperator = {'z1'};
% Opt.FreqTranslation = -33.5; % GHz

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