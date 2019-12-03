function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.Qcrit = 7;
Pulse.Frequency = 1000* [-0.1 0.1];
Pulse.tp = 0.1;

Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.mwFreq = 33.5;
Exp.DetSequence = 1; 

% Options ---------------------------
Exp.DetOperator = {'z1'};

% To test I --------------------------
Exp.nPoints = [3];
Exp.Dim1 = {'p1.Frequency' [-0.05 0]};

[~, signal1] = spidyan(Sys,Exp);

% To test II --------------------------
Exp.Dim1 = {'p1.Frequency(1)' -0.05};

[~, signal2] = spidyan(Sys,Exp);

% To test III -------------------------
Exp.Dim1 = {'p1.Frequency(1)' [0 -0.05 -0.1]};

[~, signal3] = spidyan(Sys,Exp);

err = ~areequal(signal1,signal2,1e-4,'abs') || ~areequal(signal1,signal3,1e-4,'abs');

data = [];