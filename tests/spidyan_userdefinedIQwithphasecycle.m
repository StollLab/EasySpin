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
Pulse.Frequency = 0.100;
Pulse.tp = 0.1;

PC = [0, 1; pi, -1];

% Experiment with pulses internally defined -------------------
Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.mwFreq = 33.5;
Exp.DetEvents = [1 1 1]; 
Exp.PhaseCycle = {PC []};

% Detection -------------------------
Opt.DetOperator = {'z1'};
Opt.FreqTranslation = [0]; 

% Options ---------------------------
Opt.FrameShift = 32;

% Function Call -----------------------------

[~,signal1] = spidyan(Sys,Exp,Opt);

% Build custom IQ
Pulse.Phase = 0;
Pulse.Frequency = Pulse.Frequency*1000; % GHz to MHz
[t,IQ1] = pulse(Pulse);
Opt.dt = Exp.TimeStep/10;
[~, IQ1] = rfmixer(t,IQ1,Exp.mwFreq,'IQshift',Opt);

Pulse.Phase = pi;
[t,IQ2] = pulse(Pulse);
Opt.dt = Exp.TimeStep/10;
[t, IQ2] = rfmixer(t,IQ2,Exp.mwFreq,'IQshift',Opt);

IQPulse.IQ = [IQ1; IQ2];
IQPulse.t = t;

%  Exp with userdefined IQ
Exp2 = Exp;
Exp2.Sequence{1}= IQPulse;

[~, signal2] = spidyan(Sys,Exp2,Opt);

Exp3 = Exp;
IQPulse2.IQ = {[IQ1; IQ2]};
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

