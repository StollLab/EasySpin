function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = [1/2];

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us

Exp.t = [0.1 0.5];
Exp.Pulses = {Pulse 0};
Exp.Field = 1195; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.100 0.100];
Exp.Flip = [pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = [1 1]; 

% Options ---------------------------
Opt.FrameShift = 32;

% Isotropic
Sys.g = [gfree];
[~, signal1] = spidyan(Sys,Exp,Opt);

% Axial
Sys.g = [gfree gfree];
[~, signal2] = spidyan(Sys,Exp,Opt);

% Orthorombic
Sys.g = [gfree gfree gfree];
[~, signal3] = spidyan(Sys,Exp,Opt);

% full tensor
Sys.g = diag([gfree gfree gfree]);
[~, signal4] = spidyan(Sys,Exp,Opt);

if any([~isequal(signal1,signal2) ~isequal(signal1,signal3) ~isequal(signal1,signal4)])
  err = 1;
else
  err = 0;
end

data = [];

