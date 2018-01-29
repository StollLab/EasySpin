function [err,data] = test(opt,olddata)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = [1/2];
Sys.ZeemanFreq = [33.500];


% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us


Exp.t = [0.1 0.5 0.1];
Exp.Pulses = {Pulse 0 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.100 0.100];
Exp.Flip = [pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = [1 1 1]; 

% Detection -------------------------
Opt.DetOperator = {'z1','+1'};
Opt.FreqTranslation = [0 -33.5]; 

% Options ---------------------------
Opt.FrameShift = 32;

% Function Call -----------------------------

[~, ~, out1] = spidyan(Sys,Exp,Opt);

Opt.FrameShift = -32;

[~, ~, out2] = spidyan(Sys,Exp,Opt);

if ~areequal(out1.FinalState,out2.FinalState,1e-4)
  err = 1;
else
  err = 0;
end

data = [];

warning(orig_state);

