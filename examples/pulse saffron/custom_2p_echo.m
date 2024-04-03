% Two-pulse echo transient using chirp pulses (saffron)
%==========================================================================
% Simulates a chirp echo of a nitroxide using saffron

clear

Exp.Field = 324.9;  % mT

% Spin system
Sys.g = [2.009 2.006 2.002];
Sys.Nucs = '14N';
Sys.A = [11 11 95];  % MHz
Sys.lwpp = 5;  % MHz

% Simulate frequency-sweep spectrum to set pulse excitation bands
pepper(Sys,Exp);

%%

% Define basic experiment parameters
Exp.mwFreq = 9.1;   % GHz

% Define pulses
p90.Type = 'quartersin/linear';
p90.tp = 0.200;              % pulse length, µs
p90.Flip = pi/2;             % flip angle, radians
p90.Frequency = [-120 120];  % excitation band, MHz
p90.trise = 0.030;           % rise time, µs

p180 = p90;
p180.tp = 0.100;   % µs
p180.Flip = pi;    % rad

% Define pulse sequence
tau = 0.2;  % µs
Exp.Sequence = {p90 tau p180 tau+p180.tp};

% Define detection
Exp.DetWindow = [-0.15 0.15];  % µs
Exp.DetPhase = pi;             % rad, for proper phasing of the signal

% Use grid resolution sufficient to get converged echo tails
Opt.GridSize = 30;
Opt.Verbosity = true;

saffron(Sys,Exp,Opt);
