% Three-pulse echo with phasecycling (saffron)
%==========================================================================
% creates a three pulse echo, with [+(+x)-(-x)] phase-cycles for the first
% and third pulse

clear

% Spin system
Sys.g = [2.00906 2.0062 2.0023];
Sys.A = [11.5 11.5 95]; % MHz
Sys.Nucs = '14N';

% Pulse definition
p90.tp = 0.020; % mus
p90.Flip = pi/2; % rad

tau = 0.5; % mus
T = 1; % mus

% Experiment
Exp.Sequence = {p90 tau p90 T p90 tau};
Exp.mwFreq = 9.1; % GHz
Exp.Field = 324.9; % mT
Exp.DetWindow = [-0.1 0.1]; % detection window in mus

Exp.PhaseCycle{1} = [0, 1; pi, -1]; % [+(+x)-(-x)]  phase cycle
Exp.PhaseCycle{3} = [0, 1; pi, -1]; % [+(+x)-(-x)]  phase cycle

saffron(Sys,Exp);
