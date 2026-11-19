function [ok,data] = test(opt,olddata) %#ok<FNDEF>

% Chirp pulse with circular polarization and relaxation, Sz detected over
% the whole sequence (single orientation). The last point of the transient
% equals the single point detected at the end of the sequence.

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
Exp.DetPhase = 0;

% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.SimFreq = 32;
Opt.Relaxation = 1;

% Testing Complex excitation

Exp.mwPolarization = 'circular';
Exp.DetWindow = 0;
signalEnd = saffron(Sys,Exp,Opt);

Exp.DetWindow = 'sequence';
[~, signalcomplex] = saffron(Sys,Exp,Opt);

data.signalcomplex = signalcomplex;

if ~isempty(olddata)
  ok = [areequal(signalcomplex(end),signalEnd,1e-10,'abs') ...
    areequal(signalcomplex(:),olddata.signalcomplex(:),1e-4,'abs')];
else
  ok = [];
end

warning(orig_state);
