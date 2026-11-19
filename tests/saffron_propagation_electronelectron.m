function [ok,data] = test(opt,olddata) %#ok<FNDEF>

% Two electron spins, chirp pulse, detection over the whole sequence
% (single orientation)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = [1/2 1/2];
Sys.g = [gfree gfree];
Sys.J = 0;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.tp = 0.1;
Pulse.Frequency = 1000* [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5};
Exp.Field = 1195;
Opt.IntTimeStep = 0.0001; % µs
Exp.mwFreq = 33.5;
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';
Exp.DetPhase = 0;

% Options ---------------------------
Opt.SimFreq = 32;

[t1, signal1] = saffron(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

if ~isempty(olddata)
  ok = areequal(signal1(:),olddata.signal1(:),1e-4,'abs');
else
  ok = [];
end

warning(orig_state);
