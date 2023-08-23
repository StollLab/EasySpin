% rabi oscillation (spidyan)
%==========================================================================
% creates a rabi oscillation for a single spin, using a detection window of
% length 0 and one indirect dimension
clear

% Spin System
Sys.ZeemanFreq = 9.500; % resonance frequency of spin in GHz

% eperiment setup
Pulse.Type = 'rectangular';
Pulse.tp = 0.001; % µs
Pulse.Amplitude = 30; % MHz

Exp.Sequence = {Pulse}; 
Exp.mwFreq = 9.5; % GHz

Exp.DetWindow = 0;
Exp.DetOperator = {'z1'};

Exp.nPoints = 99;
Exp.Dim1 = {'p1.tp' 0.001}; % mus

spidyan(Sys,Exp);
ylim([-1 1])