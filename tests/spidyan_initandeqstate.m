function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;
Sys.T1 = 1;
Sys.T2 = 0.4;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us

Exp.t = [0.1 0.5];
Exp.Pulses = {Pulse 0};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.100 0.100];
Exp.Flip = [pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = 1; 

% Options ---------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;
Opt.SimulationMode = 'FrameShift';

% Define initial state
Sys.initState = -sop(Sys.S,'z');

[~, signal1] = spidyan(Sys,Exp,Opt);

data.signal1 = signal1;

% Define initial state
Sys.eqState = sop(Sys.S,'z');
Opt.Relaxation = 1;

[~, signal2] = spidyan(Sys,Exp,Opt);

data.signal2 = signal2;

if ~isempty(olddata)
  err = [~areequal(signal1,olddata.signal1,1e-4) ...
    ~areequal(signal2,olddata.signal2,1e-4)];
else
  err = [];
end

