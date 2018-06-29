function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'rectangular';
Pulse.tp = 0.1;
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.mwFreq = 33.5;
Exp.DetSequence = 1; 

% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.SimFreq = 32;

% To test --------------------------
Exp.Resonator.ResonatorFrequency = 33.5;
Exp.Resonator.ResonatorQL = 300;

Exp.nPoints = 3;
Exp.Dim1 = {'d1' 0.05;
            'p1.t,p2.tp' 0.03};
          
[t, signal] = spidyan(Sys,Exp,Opt);

data.t = t;
data.signal = signal;

if ~isempty(olddata)
  for i = 1 : length(signal)
    err(i) = ~areequal(signal{i},olddata.signal{i},1e-4);
  end
else
  err = [];
end

