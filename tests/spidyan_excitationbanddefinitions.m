function [err,data] = test(opt,olddata)

% Tests if both ways of inputting excitation bands work/give the same
% result

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us

Exp.t = [0.1 0.5 0.1];
Exp.Pulses = {Pulse 0 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.00001; % us
Exp.Flip = [pi pi];

% First method ----------------------
Exp.Frequency = [-0.100 0.100];
Exp.mwFreq = 33.5;

% Detection -------------------------
Opt.DetOperator = {'z1','+1'};
Opt.FreqTranslation = [0 -33.5]; 


% Function Call -----------------------------

[~, ~, out1] = spidyan(Sys,Exp,Opt);

% Second Method -------------------------
Exp = rmfield(Exp,'mwFreq');
Exp.Frequency = [33.400 33.600];

[~, ~, out2] = spidyan(Sys,Exp,Opt);

if ~areequal(out1.FinalState,out2.FinalState,1e-4)
  err = 1;
else
  err = 0;
end

data = [];

