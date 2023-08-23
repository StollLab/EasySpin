function [ok,data] = test(opt,olddata)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = [1/2];
Sys.ZeemanFreq = [33.400];

% Experiment -------------------
Pulse.Type = 'rectangular';
Pulse.tp = 0.1;
Pulse.Flip = pi;


Exp.Sequence = {Pulse};
Exp.Field = 1240; 
Opt.IntTimeStep = 0.0001; % µs
Exp.mwFreq = 33.5;
Exp.DetSequence = 1; 
Exp.DetPhase = 0;

Exp.ResonatorFrequency = 33.5;
Exp.ResonatorQL = 300;

% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.SimFreq = 32;

% Function Call -----------------------------

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

Exp.Sequence = {Pulse 0.5};
Exp.ResonatorMode = 'compensate';

[t2, signal2] = spidyan(Sys,Exp,Opt);

data.t2 = t2;
data.signal2 = signal2;

Exp = rmfield(Exp,'ResonatorMode');
Exp = rmfield(Exp,'ResonatorFrequency');
Exp = rmfield(Exp,'ResonatorQL');

FrequencyResponse(1,:) = 33:0.0001:34;
FrequencyResponse(2,:) = ones(1,length(FrequencyResponse));
Exp.FrequencyResponse = FrequencyResponse;

[t3, signal3] = spidyan(Sys,Exp,Opt);

data.t3 = t3;
data.signal3 = signal3;

if ~isempty(olddata)
  ok = [areequal(signal1,olddata.signal1,1e-4,'abs') ...
    areequal(signal2,olddata.signal2,1e-4,'abs') areequal(signal3,olddata.signal3,1e-4,'abs')];
else
  ok = [];
end

warning(orig_state);
