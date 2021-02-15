% HYSCORE with non-ideal pulses
%==================================================
% This example illustrates how a HYSCORE with finite-length
% pulses can be simulated using a custom sequence instead
% of the built-in.

clear
Sys.Nucs = '14N';
Sys.A = [3 2]; % MHz
Sys.Q = 1; % MHz

tau = 0.1; % mus 
dt = 0.08; % mus

p90.Flip = pi/2; % rad
p90.tp = 0.01; % mus

p180.Flip = pi; % rad
p180.tp = 0.01; % mus

Exp.Field = 350; % mT

Exp.Sequence = {p90 tau p90 0 p180 0 p90 tau};

Exp.nPoints = [101 101]; % data points
Exp.Dim1 = {'d2' dt}; % increment in dimension 1, mus
Exp.Dim2 = {'d3' dt}; % increment in dimension 2, mus

Opt.GridSize = 91;
Opt.nOffsets = 30;
Opt.lwOffset = 100; % MHz

saffron(Sys,Exp,Opt);