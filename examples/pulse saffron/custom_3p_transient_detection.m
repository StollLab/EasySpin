% Three-pulse ESEEM sequence with user-defined sequence
%=======================================================
% This example illustrates how to use Exp.Sequence, Exp.Dim1,
% and Exp.nPoints to create a custom pulse sequence that can
% simulated using the fast method

clear,clc,clf
tic
Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A_ = [5 2];

tau = 0.4;
p90.Flip = pi/2;
p90.tp = 0.02;

Exp.Field = 350;
Exp.Sequence = {p90 tau p90 0 p90 tau};
Exp.mwFreq = 9.8;
Exp.DetWindow = [-0.1 0.1];
Exp.DetIntegrate = 1;

Exp.nPoints = 512;
Exp.Dim1 = {'d2', 0.005};

Opt.nKnots = 50;

saffron(Sys,Exp,Opt);
toc