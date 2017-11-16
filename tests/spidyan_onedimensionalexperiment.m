function [err,data] = test(opt,olddata)

clear Sys Exp Vary Opt Pulse sigmas Det

Sys.S = [1/2];
Sys.ZeemanFreq = [33.500];

Pulse.Type = 'rectangular';

Exp.t = [0.1 0.5 0.1];
Exp.Pulses = {Pulse 0 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = 0;
Exp.Flip = [pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = 1; 

Exp.nPoints = 3;
Exp.Dim = {'p1.Flip' 0.05};

Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;

[t1, signal1] = spidyan(Sys,Exp,Opt);

data.t1 = t1;
data.signal1 = signal1;

Opt.Relaxation = true;
Sys.T1 = 1;
Sys.T2 = 0.5;

[t3, signal3] = spidyan(Sys,Exp,Opt);

data.t3 = t3;
data.signal3 = signal3;

Opt.Relaxation = false;

Exp = rmfield(Exp,'Dim');

Exp.Resonator.nu0 = 33.5;
Exp.Resonator.QL = 300;

Exp.Dim1 = {'d1' 0.05;
            'p1.t,p2.tp' 0.03};
          
[t2, signal2] = spidyan(Sys,Exp,Opt);

data.t2 = t2;
data.signal2 = signal2;


if ~isempty(olddata)
  err = [~areequal(signal1,olddata.signal1,1e-4) ~areequal(signal3,olddata.signal3,1e-4)];
  a = length(err);
  for i = 1 : length(signal2)
    err(i+a) = ~areequal(signal2{i},olddata.signal2{i},1e-4);
  end
else
  err = [];
end

