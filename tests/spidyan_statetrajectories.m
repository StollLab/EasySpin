function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500; % GHz

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.tp = 0.1;
Pulse.Frequency = [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240; % mT
Exp.TimeStep = 0.0001; % us

Exp.mwFreq = 33.5;
Exp.DetSequence = [1 0 0]; 

% Options ---------------------------
Opt.FrameShift = 32; % GHz
Opt.SimulationMode = 'FrameShift';

% State Trajectories ---------------------------
Opt.StateTrajectories = [1 0 1];

[~, ~, out1] = spidyan(Sys,Exp,Opt);

data.sigmas1 = out1.StateTrajectories;

% State Trajectories ---------------------------
Opt.StateTrajectories = [1 1 1];
[~, ~, out2] = spidyan(Sys,Exp,Opt);

data.sigmas2 = out2.StateTrajectories;

% Three elements in the StateTrajectories are tested, the 1st, last and
% ntest:
ntest = round(length(data.sigmas1));

if ~isempty(olddata)
  err(1) = any(~areequal(out1.StateTrajectories{1},olddata.sigmas1{1},1e-10) && ~areequal(out1.StateTrajectories{ntest},olddata.sigmas1{ntest},1e-10) && ~areequal(out1.StateTrajectories{end},olddata.sigmas1{end},1e-10));
  err(2) = any(~areequal(out2.StateTrajectories{1},olddata.sigmas2{1},1e-10) && ~areequal(out2.StateTrajectories{ntest},olddata.sigmas2{ntest},1e-10) && ~areequal(out2.StateTrajectories{end},olddata.sigmas2{end},1e-10));
else
  err = [];
end
