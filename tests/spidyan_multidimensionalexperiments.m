function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.Qcrit = 7;

Exp.t = [0.2 0.5 0.2];
Exp.Pulses = {Pulse 0 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.1 0.1];
Exp.mwFreq = 33.5;
Exp.DetEvents = 1; 

% To test I --------------------------
Exp.nPoints = [3 3];
Exp.Dim1 = {'p1.trise' 0.05};
Exp.Dim2 = {'p2.Qcrit' 0.05};

% Options ---------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;
Opt.SimulationMode = 'FrameShift';

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

% To test II --------------------------
Exp.Dim1 = {'d1' [0 0.2 0.8]
            'p1.tp' 0.05};
Exp.Dim2 = {'p1.Qcrit' -1};

[t2, signal2] = spidyan(Sys,Exp,Opt);

data.t2 = t2;
data.signal2 = signal2;

if ~isempty(olddata)
  err = ~areequal(signal1,olddata.signal1,1e-4);
  for i = 1 : length(signal2)
    err(i+1) = ~areequal(signal2{i},olddata.signal2{i},1e-4);
  end
else
  err = [];
end