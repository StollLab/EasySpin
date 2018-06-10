clear Exp Sys Opt Pulse
% This script uses different types of inversion pulses. Their inversion 
% behavior is shown through the trajcetory of a single spin. The spin is on
% resonance for the monochromatic pulse and in the center of the sweep
% for the frequency modulated pulses (linear chirp and hyperbolic secant).
% The monochromatic rectangular pulse is much shorter than the two
% frequency swept pules.

% Spin System
Sys.S = [1/2]; 
Sys.ZeemanFreq = [33.500];

% Options
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;
Opt.SimulationMode = 'ShiftFrame';

%% Experiment Definition for a monochromatic rectangular pulse

% Pulse definition
Monochromatic.Type = 'rectangular';

% Experiment/Sequence - A single monochromatic pulse
Exp.t = 0.06; % us
Exp.Pulses = {Monochromatic};
Exp.Field = 1240; % mT
Exp.TimeStep = 0.0001; % us
Exp.Frequency = 0; % GHz
Exp.Flip = pi;
Exp.mwFreq = 33.5; % GHz
Exp.DetEvents = 1;

% Run simulation
[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

% Plotting
figure(1)
clf
plot(TimeAxis*1000,real(Signal));
xlabel('t [ns]')
axis tight
ylim([-1 1])

%% Experiment structure for a linear chirp pulse with quartersin amplitude modulation and a bandwidth of 100 MHz
LinearChirp.Type = 'quartersin/linear';
LinearChirp.trise = 0.030;

Exp.t = 0.200; % us
Exp.Pulses = {LinearChirp};
Exp.Field = 1240; % mT
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.05 0.05]; % excitation band, GHz
Exp.Flip = pi;
Exp.mwFreq = 33.5; % GHz
Exp.DetEvents = 1;

% Run simulation
[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

% Plotting
figure(1)
hold on
plot(TimeAxis*1000,real(Signal));
xlabel('t [ns]')
axis tight
%% Experiment structure for a hyperbolic secant pulse with a bandwidth of 100 MHz
HS.Type = 'sech/tanh';
HS.beta = 10;
HS.n = 1;

Exp.t = 0.200; % us
Exp.Pulses = {HS};
Exp.Field = 1240; % mT
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.05 0.05]; % excitation band, GHz
Exp.Flip = pi;
Exp.mwFreq = 33.5; % GHz
Exp.DetEvents = 1;

% Run simulation
[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

% Plotting
figure(1)
hold on
plot(TimeAxis*1000,real(Signal));
xlabel('t [ns]')
axis tight

legend('Monochromatic','Linear Chirp','Hyperbolic Secant')