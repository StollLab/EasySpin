function [err,data] = test(opt,olddata)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = [1/2];
Sys.ZeemanFreq = [33.600];


% Pulse -------------------
Pulse.Type = 'rectangular';
Pulse.TimeStep = 0.0001;
Pulse.Flip = pi;
Pulse.Frequency = 1000* 0.100;
Pulse.tp = 0.1;

% Experiment with pulses internally defined -------------------
Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240; 
Opt.IntTimeStep = 0.0001; % us
Exp.mwFreq = 33.5;
Exp.DetSequence = [1 1 1]; 

% Detection -------------------------
Exp.DetOperator = {'z1'};

% Function Call -----------------------------
[signal1] = spidyan(Sys,Exp,Opt);

% Build custom IQ
[t,IQ] = pulse(Pulse);
Opt.dt = Opt.IntTimeStep/10;
[t, IQ] = rfmixer(t,IQ,Exp.mwFreq,'IQshift',Opt);


%  Exp with userdefined IQ
Exp2 = Exp;
IQPulse.IQ = IQ;
IQPulse.t = t;
Exp2.Sequence{1} = IQPulse;

[signal2] = spidyan(Sys,Exp2,Opt);

Exp3 = Exp;
IQPulse2.IQ = {IQ};
IQPulse2.t = {t};

Exp3.Sequence{1} = IQPulse2;

[signal3] = spidyan(Sys,Exp3,Opt);

if any([~areequal(signal1,signal2,1e-4) ~areequal(signal1,signal3,1e-4)])
  err = 1;
else
  err = 0;
end

data = [];

warning(orig_state);
