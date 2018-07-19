% Three-pulse ESEEM sequence with user-defined sequence
%=======================================================
% This example illustrates how to use Exp.Sequence, Exp.Dim1,
% and Exp.nPoints to create a custom pulse sequence that can be
% simulated using the fast method

clear,clc,clf

Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A_ = [5 2];

tau = 0.01; % mus
p180.Flip = pi/2; % rad


Exp.Field = 350; % mT
Exp.Sequence = {p180 tau p180 0 p180 tau};

Exp.nPoints = 512;
Exp.Dim1 = {'d2', 0.005}; % mus

Opt.nKnots = 50;

saffron(Sys,Exp,Opt);