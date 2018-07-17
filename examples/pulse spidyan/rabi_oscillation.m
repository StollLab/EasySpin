% rabi oscillation (spidyan)
%==========================================================================
% this demonstrates how to simulate the spin echo of a nitroxide.
% the same result can be obtained with saffron (see example
% thyme_2p_echo.m)
clear

% Spin System
Sys.ZeemanFreq = 9.500;

% eperiment setup
Pulse.Type = 'rectangular';
Pulse.tp = 0.001;
Pulse.Amplitude = 30;

Exp.Sequence = {Pulse}; 
Exp.mwFreq = 9.5; % GHz

Exp.DetWindow = 0;
Exp.DetOperator = {'z1'};

Exp.nPoints = 99;
Exp.Dim1 = {'p1.tp' 0.001};

spidyan(Sys,Exp);
ylim([-1 1])