% Rabi oscillation (saffron)
%==========================================================================
% creates a rabi oscillation for a single spin, using a detection window of
% length 0 and one indirect dimension
clear

% Spin System
Sys.S = 1/2;

% experiment setup
Pulse.Type = 'rectangular';
Pulse.tp = 0.001; % µs
Pulse.Amplitude = 30; % MHz

Exp.Sequence = {Pulse};
Exp.mwFreq = 9.5; % GHz
Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree at 9.5 GHz
Exp.MolFrame = [0 0 0]; % single orientation, molecular frame = lab frame

Exp.DetWindow = 0;
Exp.DetOperator = {'z1'};

Exp.nPoints = 99;
Exp.Dim1 = {'p1.tp' 0.001}; % µs

saffron(Sys,Exp);
ylim([-1 1])
