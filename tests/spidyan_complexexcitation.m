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

% Options ---------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;
Opt.ComplexExcitation = [1];

% Function Call -----------------------------

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

% Using a custom complex excitation operator

Opt.ExcOperator = {sop(Sys.S,'x(1|2)')+sop(Sys.S,'y(1|2)')};
Opt.ComplexExcitation = [1 1];

[~, signal2] = spidyan(Sys,Exp,Opt);

if ~isempty(olddata)
  err = [~areequal(signal1,olddata.signal1,1e-4) ~areequal(signal1,signal2,1e-4)];
else
  err = [];
end

