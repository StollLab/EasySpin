function [ok,data] = test(opt,olddata)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.tp = 0.1;
Pulse.Frequency = 1000* [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5};
Exp.Field = 1240; 
Opt.IntTimeStep = 0.0001; % µs

Exp.mwFreq = 33.5;
Exp.DetSequence = [1 0]; 
Exp.DetPhase = 0;

% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.SimFreq = 32;
Exp.mwPolarization = 'circular';

% Function Call -----------------------------

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

% A custom excitation operator equal to the default one gives the same
% result, the operator for the imaginary part (Sy) is derived from it
Opt.ExcOperator = {sop(Sys.S,'x(1|2)')};

[~, signal2] = spidyan(Sys,Exp,Opt);

% Non-Hermitian excitation operators are rejected
Opt.ExcOperator = {sop(Sys.S,'+')};
try
  spidyan(Sys,Exp,Opt);
  nonHermitianRejected = false;
catch
  nonHermitianRejected = true;
end

if ~isempty(olddata)
  ok = [areequal(signal1,olddata.signal1,1e-4,'abs') areequal(signal1,signal2,1e-12,'abs') nonHermitianRejected];
else
  ok = [];
end

warning(orig_state);
