function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'rectangular';
Pulse.tp = 0.1;
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.mwFreq = 33.5;
Exp.DetSequence = 1; 


% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.SimFreq = 32;

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

