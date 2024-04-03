% Two-pulse echo using chirp pulses (saffron)
%==========================================================================
% This example shows how to use indirect dimension. A two-pulse chirp-echo
% is simulated, as a function of the flip angle of the first pulse.

clear

% Spin system
Sys.g = [2.00906 2.0062 2.0023];
Sys.A = [11.5 11.5 95]; % MHz
Sys.Nucs = '14N';

% Pulse definitions
Chirp90.Type = 'quartersin/linear';
Chirp90.tp = 0.200; % mus
Chirp90.Flip = pi/2; % rad
Chirp90.Frequency = [-120 120]; % excitation band, MHz
Chirp90.trise = 0.030; % rise time, mus

Chirp180.Type = 'quartersin/linear';
Chirp180.tp = 0.100; % mus
Chirp180.Flip = pi; % rad
Chirp180.Frequency = [-120 120]; % excitation band, MHz
Chirp180.trise = 0.030; % rise time, mus

% Experiment
Exp.mwFreq = 9.1; % GHz
Exp.Field = 324.9; % mT

% Pulse sequence and detection
tau = 0.5;
Exp.Sequence = {Chirp90 tau Chirp180 tau+Chirp180.tp};
Exp.DetWindow = [-0.02 0.02]; % detection window, mus

% Indirect dimensions
Exp.nPoints = 10; % number of data points in indirect dimension
Exp.Dim1 = {'p1.Flip' -pi/20}; % change flip angle of 1st pulse, in rad

Opt.GridSize = 10; % increase for better accuracy
Opt.Verbosity = 1; 

saffron(Sys,Exp,Opt);