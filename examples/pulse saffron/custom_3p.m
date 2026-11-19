% User-defined three-pulse ESEEM experiment (saffron)
%==========================================================================
% This example illustrates how to define a three-pulse ESEEM experiment
% yourself: the pulse sequence in Exp.Sequence, and the incrementation
% scheme in Exp.nPoints and Exp.Dim1. With ideal pulses, saffron simulates
% it with the fast subspace method. The result is the same as with the
% predefined experiment '3pESEEM'.

clear

% Spin system
Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A_ = [5 2];  % MHz

% Pulse sequence: pi/2 - tau - pi/2 - T - pi/2 - tau
tau = 0.01;  % µs
T = 0.06;    % initial T, µs
p90.Flip = pi/2;  % rad, ideal pulse
Exp.Field = 350;  % mT
Exp.Sequence = {p90 tau p90 T p90 tau};

% Stimulated echo: coherence orders after the three pulses
Exp.Pathways = '+0-';
Exp.DetPhase = -pi/2;  % rad, real positive echo as in the predefined '3pESEEM'

% Incrementation scheme: increment T (the second delay)
Exp.nPoints = 512;
Exp.Dim1 = {'d2', 0.005}; % µs

Opt.GridSize = 50;

saffron(Sys,Exp,Opt);
