function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.tp = 0.1;
Pulse.Frequency = [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us

Exp.mwFreq = 33.5;
Exp.DetSequence = [1 0]; 

% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.FrameShift = 32;
Opt.SimulationMode = 'FrameShift';
Opt.ComplexExcitation = 1;

% Function Call -----------------------------

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

% Using a custom complex excitation operator ----------------------
Opt.ExcOperator = {sop(Sys.S,'x(1|2)')+sop(Sys.S,'y(1|2)')};
Opt.ComplexExcitation = [1 1];

[~, signal2] = spidyan(Sys,Exp,Opt);

if ~isempty(olddata)
  err = [~areequal(signal1,olddata.signal1,1e-4) ~areequal(signal1,signal2,1e-4)];
else
  err = [];
end

