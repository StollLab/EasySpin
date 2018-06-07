function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.Qcrit = 7;

Exp.t = [0.2 0.5 0.2];
Exp.Pulses = {Pulse 0 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.1 0.1];
Exp.mwFreq = 33.5;
Exp.DetEvents = 1; 

% To test I --------------------------
Exp.nPoints = [3];
Exp.Dim1 = {'p1.Frequency' [0 0; -0.05 0.05; -0.1 0.1]};

% Options ---------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;

[~, signal1] = spidyan(Sys,Exp,Opt);


% To test II --------------------------
Exp.Dim1 = {'p1.Frequency' [-0.05 0.05]};

[~, signal2] = spidyan(Sys,Exp,Opt);

if ~areequal(signal1,signal2,1e-4)
  err = 1;
else
  err = 0;
end

data = [];