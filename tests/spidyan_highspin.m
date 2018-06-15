function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = [3/2];
Sys.D = 200;
Sys.g = gfree;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.tp = 0.1;
Pulse.Frequency = [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5};
Exp.Field = 1195; 
Exp.TimeStep = 0.0001; % us

Exp.mwFreq = 33.5;
Exp.DetSequence = [1 1]; 

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

