% Three-pulse ESEEM sequence with user-defined sequence
%=======================================================
% This example illustrates how to use Exp.Flip, Exp.Inc,
% and Exp.t to create a custom pulse sequence.

clear,clc,clf

Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A_ = [5 2];

Exp.Field = 350;
Exp.Flip = [1 1 1];
Exp.Inc = [0 1 0];
tau = 0.01;
Exp.t = [tau 0 tau];
Exp.dt = 0.012;

Opt.nKnots = 91;

saffron(Sys,Exp,Opt);