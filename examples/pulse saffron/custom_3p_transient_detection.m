% Three-pulse ESEEM sequence with user-defined sequence
%=======================================================
% This example illustrates how to use Exp.Sequence, Exp.Dim1,
% and Exp.nPoints to create a custom pulse sequence that, for each
% acquisition point an echo transient is recorded and integrated.
% For this, saffron switches to the slower thyme method

clear

Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A_ = [5 2]; % MHz

tau = 0.4; % mus

p90.Flip = pi/2; % flip angle, rad
p90.tp = 0.02; % pulse length, mus

Exp.Field = 350; % mT
Exp.Sequence = {p90 tau p90 0.1 p90 tau};
Exp.mwFreq = 9.81; % GHz

Exp.DetWindow = [-0.05 0.05]; % mus
Exp.DetIntegrate = true;

Exp.nPoints = 512; % 512 data points in indirect dimension
Exp.Dim1 = {'d2', 0.005}; % increment of second delay by 5 ns per step

Exp.PhaseCycle{1} = [0, 1; pi, -1];
Exp.PhaseCycle{3} = [0, 1; pi, -1];

Opt.nKnots = 50;

saffron(Sys,Exp,Opt);