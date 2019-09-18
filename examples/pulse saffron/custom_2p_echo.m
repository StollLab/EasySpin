% Two-pulse echo using chirp pulses (saffron)
%==========================================================================
% simulates a chirp echo of a nitroxide using saffron

clear

% Spin system
Sys.g = [2.009 2.006 2.002];
Sys.A = [11 11 95]; % MHz
Sys.Nucs = '14N';
Sys.lwpp = 5;

Exp.Field = 324.9; % mT

pepper(Sys,Exp); % use pepper to obtain field-sweep spectrum to set pulses

%%

% Pulse definitions
Chirp90.Type = 'quartersin/linear';
Chirp90.tp = 0.200; % mus
Chirp90.Flip = pi/2; % rad
Chirp90.Frequency = [-120 120]; % excitation band, MHz
Chirp90.trise = 0.030; % rise time, us

Chirp180.Type = 'quartersin/linear';
Chirp180.tp = 0.100; % mus
Chirp180.Flip = pi; % rad
Chirp180.Frequency = [-120 120]; % excitation band, MHz
Chirp180.trise = 0.030; % rise time, us

tau = 0.5; % mus
Exp.Sequence = {Chirp90 tau Chirp180 tau+Chirp180.tp};
Exp.mwFreq = 9.1; % GHz
Exp.Field = 324.9; % mT
Exp.DetWindow = [-0.02 0.02]; % us
Exp.DetPhase = pi; % rad, for proper phasing of the signal

Opt.nKnots = 20;

saffron(Sys,Exp,Opt);