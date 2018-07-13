% Three-pulse echo with phasecycling (saffron)
%==========================================================================

clear, clf, clc

% Spin system
Sys.S = 1/2;
Sys.g = [2.00906 2.0062 2.0023];
Sys.A = [11.5 11.5 95];
Sys.Nucs = '14N';

% Pulse definitions
p90.tp = 0.020;
p90.Flip = pi/2;
p90.Phase = pi;

%
tau = 0.5;
T = 1;

Exp.Sequence = {p90 tau p90 T p90 tau};
Exp.mwFreq = 9.1;
Exp.Field = 324.9;
Exp.DetWindow = [-0.1 0.1];
Exp.PhaseCycle{1} = [0, 1; pi, -1];
Exp.PhaseCycle{3} = [0, 1; pi, -1];


Opt.nKnots = 30;
Opt.Verbosity = true;

saffron(Sys,Exp,Opt);
