function [err,data] = test(opt,olddata)

clear Sys Exp Vary Opt Pulse sigmas Det

% System ------------------------
Sys.S = [1/2];
Sys.ZeemanFreq = [33.400];

% Experiment -------------------
Pulse.Type = 'rectangular';

Exp.t = [0.1];
Exp.Pulses = {Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = 0;
Exp.Flip = pi;
Exp.mwFreq = 33.5;
Exp.DetEvents = 1; 

Exp.Resonator.nu0 = 33.5;
Exp.Resonator.QL = 300;

% Options ---------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;

% Function Call -----------------------------

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

Exp.t = [0.1 0.5];
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

