function [ok,data] = test(opt,olddata) %#ok<FNDEF>

% One indirect dimension over a flip angle, with relaxation, Sz detected
% over the whole sequence (single orientation)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = 1/2;

% Experiment -------------------
Pulse.Type = 'rectangular';
Pulse.tp = 0.1;
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5 Pulse};
Opt.IntTimeStep = 0.0001; % µs
Exp.mwFreq = 33.5;
Exp.Field = mhz2mt(33500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';
Exp.DetPhase = 0;

% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.SimFreq = 32;

% To test --------------------------
Exp.nPoints = 3;
Exp.Dim1 = {'p1.Flip' 0.05};
Opt.Relaxation = true;
Sys.T1 = 1;
Sys.T2 = 0.5;

% Run saffron ---------------------
[t, signal] = saffron(Sys,Exp,Opt);

data.t = t;
data.signal = signal;

if ~isempty(olddata)
  ok = areequal(signal(:),olddata.signal(:),1e-4,'abs');
else
  ok = [];
end

warning(orig_state);
