function ok = test() %#ok<FNDEF>

% The input forms of an isotropic g tensor (one value, axial, rhombic, full
% matrix, six symmetric elements) give the same signal for one electron
% (single orientation, detection over the whole sequence)

% System ------------------------
Sys.S = 1/2;

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
Sys.g = gfree;
[~, signal1] = saffron(Sys,Exp,Opt);

% Axial
Sys.g = [gfree gfree];
[~, signal2] = saffron(Sys,Exp,Opt);

% Orthorombic
Sys.g = [gfree gfree gfree];
[~, signal3] = saffron(Sys,Exp,Opt);

% full tensor
Sys.g = diag([gfree gfree gfree]);
[~, signal4] = saffron(Sys,Exp,Opt);

% symmetric tensor [xx yy zz xy xz yz]
Sys.g = [gfree gfree gfree 0 0 0];
[~, signal5] = saffron(Sys,Exp,Opt);

ok = [isequal(signal1,signal2) isequal(signal1,signal3) isequal(signal1,signal4) isequal(signal1,signal5)];
