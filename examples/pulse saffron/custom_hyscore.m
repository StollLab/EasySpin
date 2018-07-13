% HYSCORE with non-ideal pulses
%==================================================
% This example illustrates how a HYSCORE with finite-length
% pulses can be simulated using a custom sequence instead
% of the built-in.

clear
Sys.Nucs = '14N';
Sys.A = [3 2];
Sys.Q = 1;

tau = 0.1;
dt = 0.08;

p90.Flip = pi/2;
p90.tp = 0.01;
p180.Flip = pi;
p180.tp = 0.01;

Exp.Field = 350;

Exp.Sequence = {p90 tau p90 0 p180 0 p90 tau};

Exp.nPoints = [101 101];
Exp.Dim1 = {'d2' dt};
Exp.Dim2 = {'d3' dt};

Opt.nKnots = 91;
Opt.nOffsets = 30;
Opt.lwOffset = 100;

saffron(Sys,Exp,Opt);