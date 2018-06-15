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

Exp.nPoints = 3;
Exp.Dim1 = {'p1.Flip' 0.05};

% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.FrameShift = 32;
Opt.SimulationMode = 'FrameShift';

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

if ~isempty(olddata)
  err = ~areequal(signal1,olddata.signal1,1e-4);
else
  err = [];
end

