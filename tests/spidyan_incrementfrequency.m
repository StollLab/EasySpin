function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.Qcrit = 7;
Pulse.tp = 0.2;
Pulse.Frequency = [-0.1 0.1];

Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.mwFreq = 33.5;
Exp.DetSequence = 1; 

% To test I --------------------------
Exp.nPoints = [3];
Exp.Dim1 = {'p1.Frequency' [0 0; -0.05 0.05; -0.1 0.1]};

% Options ---------------------------
Exp.DetOperator = {'z1'};

[~, signal1] = spidyan(Sys,Exp);


% To test II --------------------------
Exp.Dim1 = {'p1.Frequency' [-0.05 0.05]};

[~, signal2] = spidyan(Sys,Exp);

if ~areequal(signal1,signal2,1e-4)
  err = 1;
else
  err = 0;
end

data = [];