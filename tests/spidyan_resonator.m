function [err,data] = test(opt,olddata)

clear Sys Exp Vary Opt Pulse sigmas Det

% System ------------------------
Sys.S = [1/2];
Sys.ZeemanFreq = [33.400];

% Experiment -------------------
Pulse.Type = 'rectangular';
Pulse.tp = 0.1;
Pulse.Flip = pi;


Exp.Sequence = {Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.mwFreq = 33.5;
Exp.DetEvents = 1; 

Exp.Resonator.nu0 = 33.5;
Exp.Resonator.QL = 300;

% Options ---------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;
Opt.SimulationMode = 'FrameShift';

% Function Call -----------------------------

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

Exp.Sequence = {Pulse 0.5};
Exp.Resonator.Mode = 'compensate';

[t2, signal2] = spidyan(Sys,Exp,Opt);

data.t2 = t2;
data.signal2 = signal2;

Exp = rmfield(Exp,'Resonator');
Exp.Resonator.nu = 33:0.0001:34;
Exp.Resonator.TransferFunction = ones(1,length(Exp.Resonator.nu));

[t3, signal3] = spidyan(Sys,Exp,Opt);

data.t3 = t3;
data.signal3 = signal3;

if ~isempty(olddata)
  err = [~areequal(signal1,olddata.signal1,1e-4) ~areequal(signal2,olddata.signal2,1e-4) ~areequal(signal3,olddata.signal3,1e-4)];
else
  err = [];
end

