function ok = test() %#ok<FNDEF>

% A single detection operator can be given as a string, as a cell array
% with a string, or as a matrix (single orientation, detection during the
% pulse)

% System ------------------------
Sys.S = 1/2;

% Experiment -------------------
Pulse.tp = 0.03;
Pulse.Flip = pi;

Exp.Sequence = {Pulse};
Exp.mwFreq = 33.5;
Exp.Field = mhz2mt(33500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';

Exp.DetOperator = 'z1';
sig1 = saffron(Sys,Exp);

Exp.DetOperator = {'z1'};
sig2 = saffron(Sys,Exp);

Exp.DetOperator = [1/2 0; 0 -1/2];
sig3 = saffron(Sys,Exp);

ok = [isequal(sig1,sig2) isequal(sig1,sig3)];
