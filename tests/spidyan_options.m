function [err,data] = test(opt,olddata)

clear Sys Exp Vary Opt Pulse sigmas Det

% System ------------------------
Sys.S = [1/2];
Sys.ZeemanFreq = [33.500];
Sys.T1 = 0.5;
Sys.T2 = 0.2;


% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us


Exp.t = [0.1 0.5 0.1];
Exp.Pulses = {Pulse 0 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.100 0.100];
Exp.Flip = [pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = [1 1 1]; 

% Detection -------------------------
Opt.DetOperator = {'z1','+1'};
Opt.FreqTranslation = [0 -33.5]; 

% Options ---------------------------
Opt.FrameShift = 32;

% Function Call -----------------------------
Opt.Relaxation = [1 1];
[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

Opt.ExcOperator = {[0 1; 0 0] '+1'};
Opt.StateTrajectories = [1 0 1];
Opt.DetOperator = {'z1',[0 1; 0 0]};

[t2, signal2, state, sigmas] = spidyan(Sys,Exp,Opt);

data.t2 = t2;
data.signal2 = signal2;
data.state = state;
data.sigmas = sigmas;

if ~isempty(olddata)
  err = [~areequal(signal1,olddata.signal1,1e-4) ~areequal(signal2,olddata.signal2,1e-4) ~areequal(signal2,olddata.signal2,1e-4)];
else
  err = [];
end

