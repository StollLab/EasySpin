function [ok,data] = test(opt,olddata) %#ok<FNDEF>

% A pulse moved along an indirect dimension (Position), Sz detected over
% the whole sequence (single orientation)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = 1/2;

% Experiment -------------------
Pulse.Type = 'rectangular';
Pulse.Flip = pi;
Pulse.tp = 0.1;

Exp.Sequence = {Pulse 0.5 Pulse 0.5 Pulse 0.5};
Opt.IntTimeStep = 0.0001; % µs
Exp.mwFreq = 33.5;
Exp.Field = mhz2mt(33500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';
Exp.DetPhase = 0;

Exp.nPoints = 4;
Exp.Dim1 = {'p2.Position' 0.25};

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
