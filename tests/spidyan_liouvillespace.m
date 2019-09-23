function [err,data] = test(opt,olddata)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;
Sys.T1 = 1;
Sys.T2 = 0.4;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.tp = 0.1;
Pulse.Flip = pi;
Pulse.Frequency = 1000* [-0.100 0.100];

Exp.Sequence = {Pulse 0.5};
Exp.Field = 1240; 
Opt.IntTimeStep = 0.0001; % us
Exp.mwFreq = 33.5;
Exp.DetSequence = [1 1]; 
Exp.DetPhase = 0;

% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.SimFreq = 32;
Opt.Relaxation = 1;

% Function Call -----------------------------
[~, signal1] = spidyan(Sys,Exp,Opt);

data.signal1 = signal1;

% Test without T1 -----------------
Sys1 = rmfield(Sys,'T1');
[~, signal2] = spidyan(Sys1,Exp,Opt);

data.signal2 = signal2;

% Test without T2 -----------------
Sys2 = rmfield(Sys,'T2');
Opt.Relaxation = [1 1];

[~, signal3] = spidyan(Sys2,Exp,Opt);

data.signal3 = signal3;

if ~isempty(olddata)
  err = ~areequal(signal1,olddata.signal1,1e-4,'abs') || ...
    ~areequal(signal2,olddata.signal2,1e-4,'abs') || ...
    ~areequal(signal3,olddata.signal3,1e-4,'abs');
else
  err = [];
end

warning(orig_state);