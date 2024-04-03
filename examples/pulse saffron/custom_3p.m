% Three-pulse ESEEM sequence with user-defined sequence
%=======================================================
% This example illustrates how to use Exp.Sequence, Exp.Dim1,
% and Exp.nPoints to create a custom pulse sequence that can be
% simulated using the fast method.

clear

% Spin system
Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A_ = [5 2];  % MHz

% Experiment
tau = 0.01;  % µs
p90.Flip = pi/2;  % rad
Exp.Field = 350;  % mT
Exp.Sequence = {p90 tau p90 0.06 p90 tau};
Exp.nPoints = 512;
Exp.Dim1 = {'d2', 0.005}; % µs

Opt.GridSize = 50;

saffron(Sys,Exp,Opt);
