function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.tp = 0.1;
Pulse.Frequency = 1000* [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse Pulse};
Exp.Field = 1240; 
Opt.TimeStep = 0.0001; % us

Exp.mwFreq = 33.5;
Exp.DetSequence = [1 0]; 

% Options ---------------------------
Exp.DetOperator = {'+','+(1|2)'};
Exp.DetFreq = [33.5 33.5];
Opt.SimFreq = 32;

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

