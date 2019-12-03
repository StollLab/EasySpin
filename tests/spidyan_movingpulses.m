function [err,data] = test(opt,olddata)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'rectangular';
Pulse.Flip = pi;
Pulse.tp = 0.1;

Exp.Sequence = {Pulse 0.5 Pulse 0.5 Pulse 0.5};
Exp.Field = 1240; 
Opt.IntTimeStep = 0.0001; % us
Exp.mwFreq = 33.5;
Exp.DetSequence = 1;
Exp.DetPhase = 0;

Exp.nPoints = 4;
Exp.Dim1 = {'p2.Position' 0.25};

% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.SimFreq = 32;

[t1, signal1] = spidyan(Sys,Exp,Opt);     

data.t1 = t1;
data.signal1 = signal1;

if ~isempty(olddata)
  err = ~areequal(signal1,olddata.signal1,1e-4,'abs');
else
  err = [];
end

warning(orig_state);