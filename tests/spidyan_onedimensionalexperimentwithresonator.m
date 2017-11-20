function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'rectangular';

Exp.t = [0.1 0.5 0.1];
Exp.Pulses = {Pulse 0 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = 0;
Exp.Flip = [pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = 1; 

% Options ---------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;

% To test --------------------------
Exp.Resonator.nu0 = 33.5;
Exp.Resonator.QL = 300;

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

