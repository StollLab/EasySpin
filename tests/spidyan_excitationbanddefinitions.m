function [err,data] = test(opt,olddata)

% Tests if both ways of inputting excitation bands work/give the same
% result

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.tp = 0.1;
Pulse.Flip = pi;

% First method ----------------------
Exp.mwFreq = 33.5;
Pulse.Frequency = [-0.100 0.100];

Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.00001; % us


% Detection -------------------------
Exp.DetOperator = {'z1','+1'};
Exp.DetFrequency = [0 33.5]; 


% Function Call -----------------------------

[~, ~, out1] = spidyan(Sys,Exp);

% Second Method -------------------------
Exp = rmfield(Exp,'mwFreq');
Pulse.Frequency = [33.400 33.600];
Exp.Sequence = {Pulse 0.5 Pulse};

[~, ~, out2] = spidyan(Sys,Exp);

if ~areequal(out1.FinalState,out2.FinalState,1e-4)
  err = 1;
else
  err = 0;
end

data = [];

