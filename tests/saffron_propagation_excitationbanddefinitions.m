function ok = test() %#ok<FNDEF>

% Both ways of giving the excitation band of a chirp pulse give the same
% result: relative to Exp.mwFreq, or absolute with Exp.mwFreq = 0. With a
% single pulse: in later pulses, the phases differ, since the carrier at
% Exp.mwFreq is phase-continuous, while frequencies given relative to
% Exp.mwFreq = 0 start at the beginning of each pulse. Single orientation,
% Sz and S+ detected at the end of the sequence.

% System ------------------------
Sys.S = 1/2;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.tp = 0.1;
Pulse.Flip = pi;

% First method ----------------------
Exp.mwFreq = 33.5;
Exp.Field = mhz2mt(33500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Pulse.Frequency = 1000* [-0.100 0.100];

Exp.Sequence = {Pulse 0.5};

% Detection -------------------------
Exp.DetOperator = {'z1','+1'};
Exp.DetFreq = 33.5;  % GHz, also for the second method (Exp.mwFreq = 0)
Exp.DetWindow = 0;

% Function Call -----------------------------

signal1 = saffron(Sys,Exp);

% Second Method -------------------------
Exp.mwFreq = 0;
Pulse.Frequency = 1000* [33.400 33.600];
Exp.Sequence = {Pulse 0.5};

signal2 = saffron(Sys,Exp);

ok = areequal(signal1,signal2,1e-4,'abs');
