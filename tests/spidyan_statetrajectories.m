function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500; % GHz

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us

Exp.t = [0.1 0.5 0.1]; % us
Exp.Pulses = {Pulse 0 Pulse};
Exp.Field = 1240; % mT
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.100 0.100];
Exp.Flip = [pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = [1 0]; 

% Options ---------------------------
Opt.FrameShift = 32; % GHz

% State Trajectories ---------------------------
Opt.StateTrajectories = [1 0 1];

[~, ~, out1] = spidyan(Sys,Exp,Opt);

data.sigmas1 = out1.StateTrajectories;

% State Trajectories ---------------------------
Opt.StateTrajectories = [1 1 1];
[~, ~, out2] = spidyan(Sys,Exp,Opt);

data.sigmas2 = out2.StateTrajectories;

if ~isempty(olddata)
  err(1) = ~isequal(out1.StateTrajectories,olddata.sigmas1);
  err(2) = ~isequal(out2.StateTrajectories,olddata.sigmas2);
else
  err = [];
end
