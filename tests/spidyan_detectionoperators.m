function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

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
Exp.DetEvents = [1 0 0]; 

% Options ---------------------------
Opt.FrameShift = 32;
% Detection -------------------------
Opt.DetOperator = {'z1','+1'};
Opt.FreqTranslation = [0 -33.5]; 

% Function Call -----------------------------
[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

% Test Detection operators I ---------------------------
Opt.DetOperator = {'z1',[0 1; 0 0]};

[~, signal2] = spidyan(Sys,Exp,Opt);

% Test Detection operators II ---------------------------
Opt.DetOperator = {[1/2 0; 0 -1/2],[0 1; 0 0]};

[~, signal3] = spidyan(Sys,Exp,Opt);

if ~isempty(olddata)
  err = [~areequal(signal1,olddata.signal1,1e-4) ... 
         ~areequal(signal2,olddata.signal1,1e-4) ...
         ~areequal(signal3,olddata.signal1,1e-4)];
else
  err = [];
end

