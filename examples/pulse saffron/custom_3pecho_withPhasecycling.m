% Three-pulse echo with phasecycling (saffron)
%==========================================================================
% creates a three pulse echo, with phase cycles on the first and third pulse
% that keep only the stimulated echo. The four-step cycle on the first pulse
% is needed because the simulation has no inhomogeneous broadening, so the
% unrefocused pathway does not decay as in a real sample.

clear

% Spin system
Sys.g = [2.00906 2.0062 2.0023];
Sys.A = [11.5 11.5 95]; % MHz
Sys.Nucs = '14N';

% Pulse definition
p90.tp = 0.020; % µs
p90.Flip = pi/2; % rad

tau = 0.5; % µs
T = 1; % µs

% Experiment
Exp.Sequence = {p90 tau p90 T p90 tau};
Exp.mwFreq = 9.1; % GHz
Exp.Field = 324.9; % mT
Exp.DetWindow = [-0.1 0.1]; % detection window in µs

% Phase cycle: one row per step, with the phases of the three pulses and
% the receiver phase (rad). The first pulse goes through +x, +y, -x, -y, the
% third pulse through +x, -x, with all combinations (8 steps).
[ph1,ph3] = ndgrid((0:3)*pi/2,[0 pi]);
Exp.PhaseCycle = [ph1(:) 0*ph1(:) ph3(:) -ph1(:)-ph3(:)];

saffron(Sys,Exp);
