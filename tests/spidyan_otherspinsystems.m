function [err,data] = test(opt,olddata)

clear Sys Exp Vary Opt Pulse sigmas Det

Sys.S = [3/2];
Sys.D = 200;
Sys.g = gfree;

Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us


Exp.t = [0.1 0.5];
Exp.Pulses = {Pulse};
Exp.Field = 1195; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.100 0.100];
Exp.Flip = [pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = [1 1]; 

Opt.DetOperator = {'z1'};
Opt.FreqTranslation = [0 -33.5]; 
Opt.FrameShift = 32;

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

clear Sys
Sys.S = [1/2];
Sys.g = gfree;
Sys = nucspinadd(Sys,'1H',[8 45 45]);

[t2, signal2] = spidyan(Sys,Exp,Opt);

data.t2 = t2;
data.signal2 = signal2;

if ~isempty(olddata)
  err = [~areequal(signal1,olddata.signal1,1e-4) ~areequal(signal2,olddata.signal2,1e-4)];
else
  err = [];
end

