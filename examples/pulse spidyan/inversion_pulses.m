% spin trajectories with three different inversion pulses (spidyan)
%==========================================================================
% This script uses different types of inversion pulses. Their inversion 
% behavior is shown through the trajcetory of a single spin. The spin is on
% resonance for the monochromatic pulse and in the center of the sweep
% for the frequency modulated pulses (linear chirp and hyperbolic secant).
% The monochromatic rectangular pulse is much shorter than the two
% frequency-swept pules.

clear 

% Spin System
Sys.ZeemanFreq = 9.500; % GHz

% General experiment settings
Exp.DetOperator = 'z1';
Exp.mwFreq = 9.500; % GHz

% set up figure
figure(1)
clf
hold on
xlabel('t (ns)')
ylabel('<S_z>')
axis tight
ylim([-1 1])
%% Experiment Definition for a monochromatic rectangular pulse

% Pulse definition
Monochromatic.Type = 'rectangular';
Monochromatic.tp = 0.06; % pulse length in µs
Monochromatic.Flip = pi; % flip angle of pulse in rad

% Experiment/Sequence - A single monochromatic pulse
Exp.Sequence = {Monochromatic};

% Run simulation
[TimeAxis, Signal] = spidyan(Sys,Exp);

% Plotting
plot(TimeAxis*1000,real(Signal));

%% Experiment structure for a linear chirp pulse with quartersin amplitude modulation and a bandwidth of 100 MHz
LinearChirp.Type = 'quartersin/linear';
LinearChirp.trise = 0.030; % rise time (smoothing) of the edges in the time domain
LinearChirp.tp = 0.2; % pulse length in µs
LinearChirp.Frequency = [-50 50]; % frequency band
LinearChirp.Flip = pi; % flip angle in rad

% Update sequence
Exp.Sequence = {LinearChirp};

% Run simulation
[TimeAxis, Signal] = spidyan(Sys,Exp);

% Plotting
plot(TimeAxis*1000,real(Signal));

%% Experiment structure for a hyperbolic secant pulse with a bandwidth of 100 MHz
HS.Type = 'sech/tanh';
HS.beta = 10; % truncation parameter of 'sech'
HS.n = 1; % exponent of the secant function argument
HS.tp = 0.2; % pulse length
HS.Frequency = [-50 50]; % excitation band, GHz
HS.Flip = pi; % flip angle in rad

% Update sequence
Exp.Sequence = {HS};

[TimeAxis, Signal] = spidyan(Sys,Exp);

% Plotting
plot(TimeAxis*1000,real(Signal));

legend('Monochromatic','Linear Chirp','Hyperbolic Secant')