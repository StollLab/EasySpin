% Two-pulse echo using chirp pulses (saffron)
%==========================================================================

clear, clf, clc

% Spin system
Sys.S = 1/2;
Sys.g = [2.00906 2.0062 2.0023];
Sys.A = [11.5 11.5 95];
Sys.Nucs = '14N';

% Pulse definitions
Chirp90.Type = 'quartersin/linear';
Chirp90.tp = 0.200;
Chirp90.Flip = pi/2;
Chirp90.Frequency = [-120 120]; % excitation band, MHz
Chirp90.trise = 0.030;

Chirp90.Phase = pi;

Chirp180.Type = 'quartersin/linear';
Chirp180.tp = 0.100;
Chirp180.Flip = pi;
Chirp180.Frequency = [-120 120]; % excitation band, MHz
Chirp180.trise = 0.030;

Chirp180.Phase = pi;

%
tau = 0.5;
Exp.Sequence = {Chirp90 tau Chirp180 tau+Chirp180.tp};
Exp.mwFreq = 9.1;
Exp.Field = 324.9;
Exp.DetWindow = [-0.05 0.05];


Opt.nKnots = 20;
Opt.Verbosity = true;

Exp.nPoints = 20;
Exp.Dim1 = {'p2.Flip' -pi/30};

% Exp.Dim1 = {'p1.Frequency,p2.Frequency' [+5 -5]};
saffron(Sys,Exp,Opt);
