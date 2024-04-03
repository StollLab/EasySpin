function [ok,data] = test(opt,olddata)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = [1/2];
Sys.g = gfree;
Sys = nucspinadd(Sys,'1H',[8 45 45]);

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.tp = 0.1;
Pulse.Frequency = 1000* [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5};
Exp.TimeStep = 0.0001; % µs

Exp.Field = 1195; 
Opt.IntTimeStep = 0.0001; % µs

Exp.mwFreq = 33.5;
Exp.DetSequence = [1 1]; 
Exp.DetPhase = 0;

% Options ---------------------------
Exp.DetOperator = {'z1'};
Exp.DetFreq = [0]; 
Opt.SimFreq = 32;

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

if ~isempty(olddata)
  ok = areequal(signal1,olddata.signal1,1e-4,'abs');
else
  ok = [];
end

warning(orig_state);