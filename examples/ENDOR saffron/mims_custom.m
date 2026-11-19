% Mims ENDOR of a 1H, user-defined experiment (saffron)
%==========================================================================
% Powder Mims ENDOR spectrum of an S=1/2 coupled to one proton. Instead of
% using the predefined 'MimsENDOR', this example defines the experiment
% explicitly: the pulse sequence, with an RF pulse, in Exp.Sequence, and the
% sweep of the RF frequency in Exp.rfFreq, Exp.nPoints and Exp.Dim1.
% The result is the same as in the predefined experiment.

clear, clc

% Spin system: S=1/2 with one 1H
Sys.Nucs = '1H';
Sys.A = [-2 4];      % hyperfine coupling, MHz
Sys.lwEndor = 0.05;  % ENDOR line width, MHz

% Pulses
p90.Flip = pi/2;     % ideal microwave pulse, rad
pRF.Channel = 'rf';  % RF pulse
pRF.Flip = pi;       % RF pulse flip angle, rad
pRF.tp = 20;         % RF pulse length, µs

% Pulse sequence: pi/2 - tau - pi/2 - RF - pi/2 - tau - echo
tau = 0.1;           % µs
Exp.Sequence = {p90 tau p90 0 pRF 0 p90 tau};

Exp.Field = 350;     % mT

% Stimulated echo: electron coherence orders after the four pulses
Exp.Pathways = '+00-';

% Set detection phase to get a real positive echo, as in the
% predefined 'MimsENDOR'
Exp.DetPhase = -pi/2;   % rad

% Sweep of the RF frequency
Exp.rfFreq = 11;     % initial RF frequency, MHz
Exp.nPoints = 401;
Exp.Dim1 = {'rfFreq', 0.02};  % RF frequency step, MHz

Opt.GridSize = 91;

% Simulate and plot
saffron(Sys,Exp,Opt);
