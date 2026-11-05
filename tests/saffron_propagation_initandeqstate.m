function [ok,data] = test(opt,olddata) %#ok<FNDEF>

% User-defined initial state (Sys.initState), with relaxation times in Sys
% but relaxation switched off, Sz detected over the whole sequence (single
% orientation)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = 1/2;
Sys.T1 = 1;
Sys.T2 = 0.4;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.tp = 0.1;
Pulse.Flip = pi;
Pulse.Frequency = 1000* [-0.100 0.100];

Exp.Sequence = {Pulse 0.5};
Opt.IntTimeStep = 0.0001; % µs
Exp.mwFreq = 33.5;
Exp.Field = mhz2mt(33500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';
Exp.DetPhase = 0;

% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.SimFreq = 32;
Opt.Relaxation = false;

% Define initial state
Sys.initState = -sop(Sys.S,'z');

[~, signal1] = saffron(Sys,Exp,Opt);

data.signal1 = signal1;

if ~isempty(olddata)
  ok = areequal(signal1(:),olddata.signal1(:),1e-4,'abs');
else
  ok = [];
end

warning(orig_state);
