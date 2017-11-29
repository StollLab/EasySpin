function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us

Exp.t = [0.1 0.5 0.1];
Exp.Pulses = {Pulse 0 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.100 0.100];
Exp.Flip = [pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = [1 0]; 

% Options ---------------------------
Opt.FrameShift = 32;

% State Trajectories ---------------------------
Opt.StateTrajectories = [1 0 1];

[~, ~, ~, sigmas1] = spidyan(Sys,Exp,Opt);

data.sigmas1 = sigmas1;

% State Trajectories ---------------------------
Opt.StateTrajectories = [1 1 1];
[~, ~, ~, sigmas2] = spidyan(Sys,Exp,Opt);

data.sigmas2 = sigmas2;

if ~isempty(olddata)
  err = [~isequal(sigmas1,olddata.sigmas1) ~isequal(sigmas2,olddata.sigmas2)];
else
  err = [];
end
