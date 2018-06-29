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

%% Experiment Definition for a monochromatic rectangular pulse

% Pulse definition
Monochromatic.Type = 'rectangular';
Monochromatic.tp = 0.06;
Monochromatic.Flip = pi/2;

% Experiment/Sequence - A single monochromatic pulse
Exp.Sequence = {Monochromatic};
Exp.Field = 1240; % mT
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
LinearChirp.tp = 0.2;
LinearChirp.Frequency = [-0.05 0.05];
LinearChirp.Flip = pi;

Exp.Sequence = {LinearChirp};
Exp.Field = 1240; % mT

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
HS.tp = 0.2;
HS.Frequency = [-0.05 0.05]; % excitation band, GHz
HS.Flip = pi;

Exp.Sequence = {HS};
Exp.Field = 1240; % mT
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