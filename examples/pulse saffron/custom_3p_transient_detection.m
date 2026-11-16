% Three-pulse ESEEM with integrated echo transients (saffron)
%==========================================================================
% This example defines a three-pulse ESEEM experiment with rectangular
% pulses (pulse sequence in Exp.Sequence, incrementation scheme in
% Exp.nPoints and Exp.Dim1). For each point of the indirect dimension, the
% echo transient is detected in a window around the echo and integrated.
% For this, saffron uses the slower full-space method.

clear

Sys.Nucs = '1H';
Sys.A_ = [7 3]; % MHz

tau = 0.2; % µs

p90.Flip = pi/2; % flip angle, rad
p90.tp = 0.02; % pulse length, µs

Exp.Field = 350; % mT
Exp.Sequence = {p90 tau p90 0.1 p90 tau};
Exp.mwFreq = 9.81; % GHz

Exp.DetWindow = [-0.05 0.05]; % µs
Exp.DetIntegrate = true;

Exp.nPoints = 100; % 100 data points in indirect dimension
Exp.Dim1 = {'d2', 0.005}; % increment of second delay by 5 ns per step

% Phase cycles that keep only the stimulated echo. The four-step cycle on the
% first pulse is needed because the simulation has no inhomogeneous
% broadening, so the unrefocused pathway does not decay as in a real sample.
% Each row of Exp.PhaseCycle is one step: the phases of the three pulses
% and the receiver phase (rad).
[ph1,ph3] = ndgrid((0:3)*pi/2,[0 pi]);  % first pulse +x,+y,-x,-y; third pulse +x,-x
Exp.PhaseCycle = [ph1(:) 0*ph1(:) ph3(:) -ph1(:)-ph3(:)];

Opt.GridSize = 10;
Opt.Verbosity = 1;

saffron(Sys,Exp,Opt);
