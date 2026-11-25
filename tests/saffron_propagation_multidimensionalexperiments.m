function [ok,data] = test(opt,olddata) %#ok<FNDEF>

% Two indirect dimensions over pulse shape parameters, Sz detected over the
% whole sequence (single orientation)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = 1/2;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.Qcrit = 7;
Pulse.tp = 0.2;  % µs
Pulse.Frequency = 1000* [-0.1 0.1];  % MHz

Exp.Sequence = {Pulse 0.5 Pulse};
Opt.IntTimeStep = 0.0001; % µs
Exp.mwFreq = 33.5;  % GHz
Exp.Field = mhz2mt(33500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';
Exp.DetPhase = 0;

Exp.nPoints = [3 3];
Exp.Dim1 = {'p1.trise' 0.05};
Exp.Dim2 = {'p2.Qcrit' 0.05};

% Options ---------------------------
Exp.DetOperator = {'z1'};
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
