% moving a pulse with Exp.Dim (spidyan)
%==========================================================================
% This script show how to move a pulse around using Exp.Dim1 in
% combination with the option 'Position'. The results are shown as
% trajectory of a single spin wich is on resoance with the pulses

% Spin System
Sys.ZeemanFreq = 9.500; % GHz

Pulse90.Type = 'rectangular';
Pulse90.Flip = pi/2; % rad
Pulse90.tp = 0.025; % µs

Pulse180.Type = 'rectangular';
Pulse180.Flip = pi; % rad
Pulse180.tp = 0.05; % µs

% Create a pi/2 - tau1 - pi - tau2 - pi - tau3 sequence
Exp.Sequence = {Pulse90 0.1 Pulse180 0.9 Pulse180 1};
Exp.mwFreq = 9.5; % GHz

% Move the center pulse to the back in 100 ns steps
Exp.nPoints = 9;
Exp.Dim1 = {'p2.Position',0.1};

% Detection
Exp.DetOperator = 'z1';

% Run simulation
spidyan(Sys,Exp);