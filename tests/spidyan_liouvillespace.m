function [err,data] = test(opt,olddata)

clear Sys Exp Vary Opt Pulse sigmas Det

% System ------------------------
Sys.S = [1/2];
Sys.ZeemanFreq = [33.500];
Sys.T1 = 1;
Sys.T2 = 0.4;

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
Exp.DetEvents = [1 1]; 

% Options ---------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;
Opt.Relaxation = 1;

% Function Call -----------------------------

[~, signal, ~, state1] = spidyan(Sys,Exp,Opt);

data.signal = signal;

% Test without T1 -----------------
Sys1 = rmfield(Sys,'T1');
[~, signal2] = spidyan(Sys1,Exp,Opt);

data.signal2 = signal2;

% Test without T2 -----------------
Sys2 = rmfield(Sys,'T2');
Sys2.initState = -sop(Sys2.S,'z');
Sys2.eqState = -sop(Sys2.S,'z');
Opt.Relaxation = [1 1];

[~, signal3] = spidyan(Sys2,Exp,Opt);

data.signal3 = signal3;


% Additional Options to test ---------------------------
Exp.DetEvents = [0 0]; 
[~, ~, ~, state2] = spidyan(Sys,Exp,Opt);

Opt.ComplexExcitation = 1;
Exp.DetEvents = [1 1];
[~, ~, ~, statecomplex1] = spidyan(Sys,Exp,Opt);

Exp.DetEvents = [0 0];
[~, signalcomplex, ~, statecomplex2] = spidyan(Sys,Exp,Opt);

data.signalcomplex = signalcomplex;

if ~isempty(olddata)
  err = [~areequal(signal,olddata.signal,1e-4) ~areequal(state1,state2,1e-4)];
  err = [err ~areequal(signal2,olddata.signal2,1e-4) ~areequal(signal3,olddata.signal3,1e-4)];
  err = [err ~areequal(statecomplex1,statecomplex2,1e-4) ~areequal(signalcomplex,olddata.signalcomplex,1e-4)];
else
  err = [];
end

