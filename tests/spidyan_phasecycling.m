function [err,data] = test(opt,olddata)

clear Sys Exp Vary Opt Pulse sigmas Det

% System ------------------------
Sys.S = [1/2];
Sys.ZeemanFreq = [33.500];

% Experiment -------------------
Pulse.Type = 'rectangular';

PC = [0, 1; pi, -1];

Exp.t = [0.1 0.1 0.1];
Exp.Pulses = {Pulse Pulse Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [0 0; 0 0; 0 0];
Exp.Flip = [pi pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = 1; 
Exp.PhaseCycle = {PC [] PC};
Exp.Phase = [0 pi];

% Options ---------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;

% Function Call -----------------------------

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

if ~isempty(olddata)
  err = ~areequal(signal1,olddata.signal1,1e-4);
else
  err = [];
end

