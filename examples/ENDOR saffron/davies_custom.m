% Davies ENDOR of a 1H, user-defined experiment (saffron)
%==========================================================================
% Powder Davies ENDOR spectrum of an S=1/2 coupled to one proton. Instead of
% using the predefined 'DaviesENDOR', this example defines the experiment
% explicitly: the pulse sequence, with finite microwave pulses and an RF pulse,
% in Exp.Sequence, and the sweep of the RF frequency in Exp.rfFreq,
% Exp.nPoints and Exp.Dim1. The delays before and after the RF pulse do
% not change the result, since relaxation is not modeled. The result is the
% same as in the predefined experiment.

clear, clc

% Spin system: S=1/2 with one 1H
Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = [-4 8];      % hyperfine coupling, MHz
Sys.lwAnisoG = 20;   % inhomogeneous EPR line width (Gaussian FWHM), MHz;
                     % needed for Davies ENDOR, see davies_predefined.m
Sys.lwEndor = 0.2;   % ENDOR line width, MHz

% Pulses
pInv.Flip = pi;      % selective inversion pulse, rad
pInv.tp = 0.2;       % µs
p90.Flip = pi/2;     % rad
p90.tp = 0.1;        % µs
p180.Flip = pi;      % rad
p180.tp = 0.2;       % µs
pRF.Channel = 'rf';  % RF pulse
pRF.Flip = pi;       % rad
pRF.tp = 20;         % µs

% Pulse sequence: pi - T - RF - T - pi/2 - tau - pi - tau - echo
T = 1.0;             % µs
tau = 0.3;           % µs
Exp.Field = 350.1;   % mT
Exp.mwFreq = 9.8;    % GHz
Exp.Sequence = {pInv T pRF T p90 tau p180 tau};

% Echo: electron coherence orders after the four pulses (populations after the
% inversion and RF pulses)
Exp.Pathways = '00+-';
Exp.DetPhase = -pi/2;  % rad, as in the predefined 'DaviesENDOR'

% Sweep of the RF frequency
Exp.rfFreq = 8;      % initial RF frequency, MHz
Exp.nPoints = 281;
Exp.Dim1 = {'rfFreq', 0.05};  % RF frequency step, MHz

Opt.GridSize = 181;  % fine powder integration grid

% Simulate and plot
saffron(Sys,Exp,Opt);
