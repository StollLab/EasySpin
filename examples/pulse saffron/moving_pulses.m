% Moving a pulse with Exp.Dim (saffron)
%==========================================================================
% This script shows how to move a pulse around using Exp.Dim1 in
% combination with the option 'Position'. The results are shown as the
% trajectory of a single spin which is on resonance with the pulses

clear

% Spin System
Sys.S = 1/2;

Pulse90.Type = 'rectangular';
Pulse90.Flip = pi/2; % rad
Pulse90.tp = 0.025; % µs

Pulse180.Type = 'rectangular';
Pulse180.Flip = pi; % rad
Pulse180.tp = 0.05; % µs

% Create a pi/2 - tau1 - pi - tau2 - pi - tau3 sequence
Exp.Sequence = {Pulse90 0.1 Pulse180 0.9 Pulse180 1};
Exp.mwFreq = 9.5; % GHz
Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree at 9.5 GHz
Exp.MolFrame = [0 0 0]; % single orientation, molecular frame = lab frame

% Move the center pulse to the back in 100 ns steps
Exp.nPoints = 9;
Exp.Dim1 = {'p2.Position',0.1};

% Detection during the whole sequence
Exp.DetWindow = 'sequence';
Exp.DetOperator = 'z1';

% Run simulation and plot
saffron(Sys,Exp);
