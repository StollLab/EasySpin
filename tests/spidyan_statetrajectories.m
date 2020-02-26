function [ok,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500; % GHz

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.tp = 0.1;
Pulse.Frequency = 1000* [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240; % mT
Exp.TimeStep = 0.0001; % us

Exp.mwFreq = 33.5;
Exp.DetSequence = [1 0 0]; 

% Options ---------------------------
Opt.SimFreq = 32; % GHz

% State Trajectories ---------------------------
Opt.StateTrajectories = [1 0 1];

[~, ~, out1] = spidyan(Sys,Exp,Opt);

data.sigmas1 = out1.StateTrajectories;

% State Trajectories ---------------------------
Opt.StateTrajectories = [1 1 1];
[x, y, out2] = spidyan(Sys,Exp,Opt);

data.sigmas2 = out2.StateTrajectories;

if ~isempty(olddata)
  % Three elements in the StateTrajectories are tested.
  N = numel(olddata.sigmas1);
  idx = [1 round(N/2) N];
  ok(1) = true;
  ok(2) = true;
  for k = 1:numel(idx)
    ok(1) = ok(1) && areequal(out1.StateTrajectories{idx(k)},olddata.sigmas1{idx(k)},1e-10,'abs');
    ok(2) = ok(2) && areequal(out2.StateTrajectories{idx(k)},olddata.sigmas2{idx(k)},1e-10,'abs');
  end
else
  ok = [];
end
