% HYSCORE with non-ideal pulses
%==================================================
% This example illustrates how a HYSCORE with finite-length
% pulses can be simulated using a custom sequence instead
% of the built-in.

clear
Sys.Nucs = '14N';
Sys.A = [3 2];
Sys.Q = 1;

Exp.Field = 350;
Exp.Flip = [1 1 2 1];
Exp.Inc = [0 1 2 0];
Exp.t = [0.1 0 0 0.1]*1e-3;
Exp.tp = [10 10 10 10]*1e-3;
Exp.dt = 0.08;
Exp.nPoints = 101;

Opt.nKnots = 91;
Opt.nOffsets = 30;
Opt.lwOffset = 100;

saffron(Sys,Exp,Opt);