function [ok,data] = test(opt,olddata) %#ok<FNDEF>

% S = 1/2 with two chirp pulses, Sz and S+ detected over the whole
% sequence (single orientation). The last point of the transient equals
% the single point detected at the end of the sequence.

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = 1/2;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.tp = 0.1;
Pulse.Flip = pi;
Pulse.Frequency = 1000* [-0.100 0.100];

Exp.Sequence = {Pulse 0.5 Pulse};
Opt.IntTimeStep = 0.0001; % µs
Exp.mwFreq = 33.5;
Exp.Field = mhz2mt(33500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';

% Detection -------------------------
Exp.DetOperator = {'z1','+1'};

% Options ---------------------------
Opt.SimFreq = 32;

% Function Call -----------------------------

[t, signal] = saffron(Sys,Exp,Opt);

data.t = t;
data.signal = signal;

Exp.DetWindow = 0;
signalEnd = saffron(Sys,Exp,Opt);

if ~isempty(olddata)
  ok = [areequal(signal(end,:),signalEnd,1e-10,'abs') areequal(signal(:),olddata.signal(:),1e-4,'abs')];
else
  ok = [];
end

warning(orig_state);
