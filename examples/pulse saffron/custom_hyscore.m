% User-defined HYSCORE experiment with rectangular pulses (saffron)
%==========================================================================
% This example defines a HYSCORE experiment with finite-length rectangular
% pulses yourself, instead of using the predefined 'HYSCORE': the pulse
% sequence in Exp.Sequence, and the incrementation scheme (t1 and t2) in
% Exp.nPoints, Exp.Dim1 and Exp.Dim2. Sys.lwAnisoG gives the inhomogeneous
% width of the EPR line, so that the pulses excite spins with a spread of
% resonance offsets.

clear
Sys.Nucs = '14N';
Sys.A = [3 2]; % MHz
Sys.Q = 1; % MHz
Sys.lwAnisoG = 100; % MHz, inhomogeneous EPR line width

tau = 0.1; % µs
dt = 0.08; % µs

p90.Flip = pi/2; % rad
p90.tp = 0.01; % µs

p180.Flip = pi; % rad
p180.tp = 0.01; % µs

Exp.Field = 350; % mT

% Pulse sequence: pi/2 - tau - pi/2 - t1 - pi - t2 - pi/2 - tau
Exp.Sequence = {p90 tau p90 0 p180 0 p90 tau};
Exp.Pathways = '+00-'; % echo pathways of HYSCORE, as in the predefined 'HYSCORE'
Exp.DetPhase = pi/2; % rad, real positive echo as in the predefined 'HYSCORE'

% Incrementation scheme
Exp.nPoints = [101 101]; % data points
Exp.Dim1 = {'d2' dt}; % increment of t1 in dimension 1, µs
Exp.Dim2 = {'d3' dt}; % increment of t2 in dimension 2, µs

Opt.GridSize = 91;

saffron(Sys,Exp,Opt);
