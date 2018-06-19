function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = [1/2 1/2];
Sys.J = 0;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.tp = 0.1;
Pulse.Frequency = [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5};
Exp.Field = 1195; 
Exp.TimeStep = 0.0001; % us

Exp.mwFreq = 33.5;
Exp.DetSequence = [1 1]; 

% Options ---------------------------
Opt.SimFrequency = 32;

% Isotropic
Sys.g = [gfree gfree];
[~, signal1] = spidyan(Sys,Exp,Opt);

% Axial
Sys.g = [gfree gfree; gfree gfree];
[~, signal2] = spidyan(Sys,Exp,Opt);

% Orthorombic
Sys.g = [gfree gfree gfree; gfree gfree gfree];
[~, signal3] = spidyan(Sys,Exp,Opt);

% full tensor
Sys.g = diag([gfree gfree gfree]);
Sys.g = repmat(Sys.g, [2 1]);
[~, signal4] = spidyan(Sys,Exp,Opt);

if any([~isequal(signal1,signal2) ~isequal(signal1,signal3) ~isequal(signal1,signal4)])
  err = 1;
else
  err = 0;
end

data = [];

