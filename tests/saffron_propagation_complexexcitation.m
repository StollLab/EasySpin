function [ok,data] = test(opt,olddata) %#ok<FNDEF>

% Chirp pulse with circular polarization, Sz detected during the pulse
% (single orientation)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = 1/2;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.tp = 0.1;
Pulse.Frequency = 1000* [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse};
Opt.IntTimeStep = 0.0001; % µs

Exp.mwFreq = 33.5;
Exp.Field = mhz2mt(33500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';
Exp.DetPhase = 0;

% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.SimFreq = 32;
Exp.mwPolarization = 'circular';

% Function Call -----------------------------

[t1, signal1] = saffron(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

if ~isempty(olddata)
  ok = areequal(signal1(:),olddata.signal1(:),1e-4,'abs');
else
  ok = [];
end

warning(orig_state);
