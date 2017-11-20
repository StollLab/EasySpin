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
Exp.Pulses = {Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.100 0.100];
Exp.Flip = [pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = [0 0]; 

% Options ---------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;
Opt.Relaxation = 1;

% Testing Complex excitation

Opt.ComplexExcitation = 1;
[~, ~, ~, statecomplex1] = spidyan(Sys,Exp,Opt);

Exp.DetEvents = [1 1];
[~, signalcomplex, ~, statecomplex2] = spidyan(Sys,Exp,Opt);

data.signalcomplex = signalcomplex;

if ~isempty(olddata)
  err = [~areequal(statecomplex1,statecomplex2,1e-4) ...
    ~areequal(signalcomplex,olddata.signalcomplex,1e-4)];
else
  err = [];
end

