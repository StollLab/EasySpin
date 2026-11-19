% Mims ENDOR of a 1H, predefined experiment (saffron)
%==========================================================================
% Powder Mims ENDOR spectrum of an S=1/2 coupled to one proton, using the
% predefined experiment 'MimsENDOR'. The sequence is
%   pi/2 - tau - pi/2 - T - RF pi pulse - pi/2 - tau - echo
% and the stimulated echo is recorded as a function of the RF frequency.
% The RF pulse reduces the echo at the nuclear transition frequencies, so
% the ENDOR lines appear as dips. Blind spots in the spectrum appear at
% couplings of n/tau (n = 0, 1, 2, ...).

clear, clc

% Spin system: S=1/2 with one 1H
Sys.Nucs = '1H';
Sys.A = [-2 4];      % hyperfine coupling, MHz
Sys.lwEndor = 0.05;  % ENDOR line width, MHz

% Experiment: predefined sequence
Exp.Field = 350;     % mT
Exp.Sequence = 'MimsENDOR';
Exp.tau = 0.1;       % µs
Exp.Range = [11 19]; % RF range, MHz
Exp.nPoints = 401;

Opt.GridSize = 91;   % fine powder integration grid

% Simulate and plot
saffron(Sys,Exp,Opt);
