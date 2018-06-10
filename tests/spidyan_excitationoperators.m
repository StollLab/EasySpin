function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us

Exp.t = [0.1 0.1];
Exp.Pulses = {Pulse Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.100 0.100];
Exp.Flip = [pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = [1 0]; 

% Options ---------------------------
Opt.DetOperator = {'+','+(1|2)'};
Opt.FreqTranslation = [-33.5 -33.5];
Opt.FrameShift = 32;
Opt.SimulationMode = 'FrameShift';

% Test custom exc operator syntax
Opt.ExcOperator = {sop(Sys.S,'x(1|2)')};

[t, signal1] = spidyan(Sys,Exp,Opt);

data.t = t;
data.signal1 = signal1;

Opt.ExcOperator = {'x(1|2)' []};
[~, signal2] = spidyan(Sys,Exp,Opt);

% Test custom exc operator input syntax
Opt.ExcOperator = {'x(1|2)' 'x(1|3)'};
[~, signal3] = spidyan(Sys,Exp,Opt);

data.t = t;
data.signal3 = signal3;

if ~isempty(olddata)
  err = [~areequal(signal1,olddata.signal1,1e-4) ...
    ~areequal(signal3,olddata.signal3,1e-4) ~isequal(signal1,signal2)];
else
  err = [];
end

