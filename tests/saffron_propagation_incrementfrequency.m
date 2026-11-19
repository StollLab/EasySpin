function ok = test() %#ok<FNDEF>

% An indirect dimension over the frequency band of a chirp pulse: a table
% of bands gives the same result as a linear increment of both edges
% (single orientation, Sz detected over the whole sequence)

% System ------------------------
Sys.S = 1/2;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.Qcrit = 7;
Pulse.tp = 0.2;
Pulse.Frequency = 1000* [-0.1 0.1];

Exp.Sequence = {Pulse 0.5 Pulse};
Exp.mwFreq = 33.5;
Exp.Field = mhz2mt(33500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';

% To test I --------------------------
Exp.nPoints = 3;
Exp.Dim1 = {'p1.Frequency' [0 0; -0.05 0.05; -0.1 0.1]};

% Options ---------------------------
Exp.DetOperator = {'z1'};
% simulation frame at the microwave frequency: long time steps
Opt.SimFreq = 33.5;

[~, signal1] = saffron(Sys,Exp,Opt);


% To test II --------------------------
Exp.Dim1 = {'p1.Frequency' [-0.05 0.05]};

[~, signal2] = saffron(Sys,Exp,Opt);

ok = areequal(signal1,signal2,1e-4,'abs');
