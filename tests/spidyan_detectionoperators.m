function [err,data] = test(opt,olddata)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.tp = 0.1;
Pulse.Frequency = 1000* [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240; 
Opt.IntTimeStep = 0.0001; % us

Exp.mwFreq = 33.5;
Exp.DetSequence = [1 0 0]; 

% Options ---------------------------
Opt.SimFreq = 32;
% Detection -------------------------
Exp.DetOperator = {'z1','+1'};
Exp.DetFreq = [0 33.5]; 

% Function Call -----------------------------
[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

% Test Detection operators I ---------------------------
Exp.DetOperator = {'z1',[0 1; 0 0]};

[~, signal2] = spidyan(Sys,Exp,Opt);

% Test Detection operators II ---------------------------
Exp.DetOperator = {[1/2 0; 0 -1/2],[0 1; 0 0]};

[~, signal3] = spidyan(Sys,Exp,Opt);

if ~isempty(olddata)
  err = [~areequal(signal1,olddata.signal1,1e-4) ... 
         ~areequal(signal2,olddata.signal1,1e-4) ...
         ~areequal(signal3,olddata.signal1,1e-4)];
else
  err = [];
end

warning(orig_state);