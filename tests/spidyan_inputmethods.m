function [err,data] = test(opt,olddata)

clear Sys Exp Vary Opt Pulse sigmas Det

% System ------------------------
Sys.S = [1/2];
Sys.ZeemanFreq = [33.500];

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us


Exp.t = [0.1 0.5];
Exp.Pulses = {Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.100 0.100];
Exp.Flip = [pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = [1 0]; 

% Detection -------------------------
Opt.DetOperator = {'z1','+1'};
Opt.FreqTranslation = [0 -33.5]; 

% Options ---------------------------
Opt.FrameShift = 32;

% Function Call -----------------------------

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

% Additional Options to test ---------------------------
Opt.ExcOperator = {[0 1/2; 1/2 0] 'x1'};
Opt.StateTrajectories = [1 0 1];
Opt.DetOperator = {'z1',[0 1; 0 0]};

[t2, signal2, state, sigmas] = spidyan(Sys,Exp,Opt);

data.t2 = t2;
data.signal2 = signal2;
data.state = state;
data.sigmas = sigmas;

Opt.StateTrajectories = [1 1 1];
[~, ~, ~, sigmas2] = spidyan(Sys,Exp,Opt);

data.sigmas2 = sigmas2;

Exp = rmfield(Exp,'mwFreq');
Exp.Frequency = [33.400 33.600];
Exp.TimeStep = 0.00001; % us

[t3, signal3] = spidyan(Sys,Exp,Opt);

data.signal3 = signal3;
data.t3 = t3;

if ~isempty(olddata)
  err = [~isequal(signal1,signal2) ~areequal(signal3,olddata.signal3,1e-4) ~isequal(sigmas2,olddata.sigmas2)];
  err = [err ~areequal(signal1,olddata.signal1,1e-4) ~areequal(signal2,olddata.signal2,1e-4) ~isequal(sigmas,olddata.sigmas)];
else
  err = [];
end

