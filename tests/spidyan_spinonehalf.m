function [err,data] = test(opt,olddata)

clear Sys Exp Vary Opt Pulse sigmas Det

% System ------------------------
Sys.S = [1/2];
Sys.ZeemanFreq = [33.500];


% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.tp = 0.1;
Pulse.Flip = pi;
Pulse.Frequency = [-0.100 0.100];

Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.mwFreq = 33.5;
Exp.DetSequence = [1 1 1]; 

% Detection -------------------------
Exp.DetOperator = {'z1','+1'};
Exp.DetFrequency = [0 33.5]; 

% Options ---------------------------
Opt.SimFrequency = 32;

% Function Call -----------------------------

[t, signal, out1] = spidyan(Sys,Exp,Opt);

data.t = t;
data.signal = signal;

Exp.DetEvents = [0 0 0]; 

[~, ~, out2] = spidyan(Sys,Exp,Opt);

if ~isempty(olddata)
  err = [~areequal(out1.FinalState,out2.FinalState,1e-4) ~areequal(signal,olddata.signal,1e-4)];
else
  err = [];
end

