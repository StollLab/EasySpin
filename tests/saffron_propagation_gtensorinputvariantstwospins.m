function ok = test() %#ok<FNDEF>

% The input forms of isotropic g tensors (one value per electron, axial,
% rhombic, full matrices) give the same signal for two electrons (single
% orientation, detection over the whole sequence)

% System ------------------------
Sys.S = [1/2 1/2];
Sys.J = 0;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.tp = 0.1;
Pulse.Frequency = 1000* [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5};
Exp.Field = 1195;
Exp.mwFreq = 33.5;
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';

% Options ---------------------------
Opt.SimFreq = 32;

% Isotropic
Sys.g = [gfree gfree];
[~, signal1] = saffron(Sys,Exp,Opt);

% Axial
Sys.g = [gfree gfree; gfree gfree];
[~, signal2] = saffron(Sys,Exp,Opt);

% Orthorombic
Sys.g = [gfree gfree gfree; gfree gfree gfree];
[~, signal3] = saffron(Sys,Exp,Opt);

% full tensor
Sys.g = diag([gfree gfree gfree]);
Sys.g = repmat(Sys.g, [2 1]);
[~, signal4] = saffron(Sys,Exp,Opt);

ok = [isequal(signal1,signal2) isequal(signal1,signal3) isequal(signal1,signal4)];
