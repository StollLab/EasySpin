function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;


% Experiment -------------------
Pulse.Type = 'rectangular';

Exp.t = [0.1 0.5 0.1];
Exp.Pulses = {Pulse 0 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = 0;
Exp.Flip = [pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = 1; 

Exp.nPoints = 3;
Exp.Dim = {'p1.Flip' 0.05};

% Options ---------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

if ~isempty(olddata)
  err = ~areequal(signal1,olddata.signal1,1e-4);
else
  err = [];
end

