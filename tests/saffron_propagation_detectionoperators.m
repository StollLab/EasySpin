function [ok,data] = test(opt,olddata) %#ok<FNDEF>

% Detection operators given as strings or as matrices give the same
% signals (single orientation, detection during the first pulse)

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

Exp.Sequence = {Pulse};  % later events don't affect the detected signal
Opt.IntTimeStep = 0.0001; % µs

Exp.mwFreq = 33.5;
Exp.Field = mhz2mt(33500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';

% Options ---------------------------
Opt.SimFreq = 32;
% Detection -------------------------
Exp.DetOperator = {'z1','+1'};

% Function Call -----------------------------
[t1, signal1] = saffron(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

% Test Detection operators I ---------------------------
Exp.DetOperator = {'z1',[0 1; 0 0]};

[~, signal2] = saffron(Sys,Exp,Opt);

% Test Detection operators II ---------------------------
Exp.DetOperator = {[1/2 0; 0 -1/2],[0 1; 0 0]};

[~, signal3] = saffron(Sys,Exp,Opt);

if ~isempty(olddata)
  ok = [areequal(signal1(:),olddata.signal1(:),1e-4,'abs') ...
        areequal(signal2(:),olddata.signal1(:),1e-4,'abs') ...
        areequal(signal3(:),olddata.signal1(:),1e-4,'abs')];
else
  ok = [];
end

warning(orig_state);
