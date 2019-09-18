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
Exp.DetSequence = 1; 
Exp.DetPhase = 0;

% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.SimFreq = 32;
Opt.Relaxation = false;

% Define initial state
Sys.initState = -sop(Sys.S,'z');

[~, signal1] = spidyan(Sys,Exp,Opt);

data.signal1 = signal1;

% Define initial state
Sys.eqState = sop(Sys.S,'z');
Opt.Relaxation = 1;

[~, signal2] = spidyan(Sys,Exp,Opt);

data.signal2 = signal2;

if ~isempty(olddata)
  err = [~areequal(signal1,olddata.signal1,1e-4) ...
    ~areequal(signal2,olddata.signal2,1e-4)];
else
  err = [];
end

warning(orig_state);