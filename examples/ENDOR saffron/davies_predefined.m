% Davies ENDOR of a 1H, predefined experiment (saffron)
%==========================================================================
% Powder Davies ENDOR spectrum of an S=1/2 coupled to one proton, using the
% predefined experiment 'DaviesENDOR'. The sequence is
%   pi - T - RF pi pulse - pi/2 - tau - pi - tau - echo
% The first microwave pulse is selective and inverts only a narrow part of
% the EPR line. The RF pulse reduces the inversion at the nuclear
% transition frequencies, so the ENDOR lines appear as peaks on the
% negative (inverted) echo. Davies ENDOR works best for couplings larger
% than about 1/tp of the inversion pulse.
%
% Sys.lwAnisoG gives the EPR line an inhomogeneous width, so the spins have
% a spread of resonance frequencies, and for many of them the first pulse
% inverts one hyperfine line but not the other. Without Sys.lwAnisoG, the
% pulse inverts both hyperfine lines equally, and the ENDOR spectrum is flat.

clear, clc

% Spin system: S=1/2 with one 1H
Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = [-4 8];      % hyperfine coupling, MHz
Sys.lwAnisoG = 20;   % inhomogeneous EPR line width (Gaussian FWHM), MHz
Sys.lwEndor = 0.2;   % ENDOR line width, MHz

% Experiment
Exp.mwFreq = 9.8;    % GHz
Exp.Field = 350.1;   % mT

Exp.Sequence = 'DaviesENDOR';
Exp.tau = 0.3;       % µs
Exp.tp = [0.2 0.1 0.2];  % lengths of the three microwave pulses, µs

Exp.Range = [8 22];  % RF sweep range, MHz
Exp.nPoints = 281;

Opt.GridSize = 181;  % fine powder integration grid

% Simulate and plot
saffron(Sys,Exp,Opt);
