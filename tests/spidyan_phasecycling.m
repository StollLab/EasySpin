function [err,data] = test(opt,olddata)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = [1/2];
Sys.ZeemanFreq = [33.500];

% Experiment -------------------
Pulse.Type = 'rectangular';
Pulse.tp = 0.1;
Pulse.Flip = pi;

Pulse2 = Pulse;
Pulse2.Phase = pi;

PC = [0, 1; pi, -1];

Exp.Sequence = {Pulse Pulse2 Pulse};
Exp.Field = 1240; 
Opt.TimeStep = 0.0001; % us
Exp.mwFreq = 33.5;
Exp.DetSequence = 1; 
Exp.PhaseCycle = {PC [] PC};
Exp.DetPhase = 0;

% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.SimFreq = 32;

% Function Call -----------------------------

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

if ~isempty(olddata)
  err = ~areequal(signal1,olddata.signal1,1e-4);
else
  err = [];
end

warning(orig_state);