function ok = test() %#ok<FNDEF>

% An indirect dimension over one edge of the frequency band of a chirp
% pulse: increments of both edges, of one edge, and a list of values give
% the same result (single orientation, Sz detected over the whole sequence)

% System ------------------------
Sys.S = 1/2;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.Qcrit = 7;
Pulse.Frequency = 1000* [-0.1 0.1];
Pulse.tp = 0.1;

Exp.Sequence = {Pulse 0.5 Pulse};
Exp.mwFreq = 33.5;
Exp.Field = mhz2mt(33500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';

% Options ---------------------------
Exp.DetOperator = {'z1'};

% To test I --------------------------
Exp.nPoints = 3;
Exp.Dim1 = {'p1.Frequency' [-0.05 0]};

[~, signal1] = saffron(Sys,Exp);

% To test II --------------------------
Exp.Dim1 = {'p1.Frequency(1)' -0.05};

[~, signal2] = saffron(Sys,Exp);

% To test III -------------------------
Exp.Dim1 = {'p1.Frequency(1)' [0 -0.05 -0.1]};

[~, signal3] = saffron(Sys,Exp);

ok = [areequal(signal1,signal2,1e-4,'abs') areequal(signal1,signal3,1e-4,'abs')];
