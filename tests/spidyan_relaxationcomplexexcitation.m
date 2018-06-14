function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;
Sys.T1 = 1;
Sys.T2 = 0.4;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.tp = 0.1;
Pulse.Flip = pi;
Pulse.Frequency = [-0.100 0.100];

Exp.Sequence = {Pulse 0.5};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us

Exp.mwFreq = 33.5;
Exp.DetEvents = [0 0]; 

% Options ---------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;
Opt.Relaxation = 1;
Opt.SimulationMode = 'FrameShift';

% Testing Complex excitation

Opt.ComplexExcitation = 1;
[~, ~, out1] = spidyan(Sys,Exp,Opt);

Exp.DetEvents = [1 1];
[~, signalcomplex, out2] = spidyan(Sys,Exp,Opt);

data.signalcomplex = signalcomplex;

if ~isempty(olddata)
  err = [~areequal(out1.FinalState,out2.FinalState,1e-4) ...
    ~areequal(signalcomplex,olddata.signalcomplex,1e-4)];
else
  err = [];
end

