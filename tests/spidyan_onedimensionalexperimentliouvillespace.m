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


% Options ---------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;
Opt.SimulationMode = 'FrameShift';

% To test --------------------------
Exp.nPoints = 3;
Exp.Dim1 = {'p1.Flip' 0.05};
Opt.Relaxation = true;
Sys.T1 = 1;
Sys.T2 = 0.5;

% Run spidyan ---------------------
[t, signal] = spidyan(Sys,Exp,Opt);

data.t = t;
data.signal = signal;

if ~isempty(olddata)
  err = ~areequal(signal,olddata.signal,1e-4);
else
  err = [];
end

