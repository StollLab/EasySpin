function ok = test()

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
Opt.IntTimeStep = 0.0001; % µs
Exp.mwFreq = 33.5;
Exp.DetSequence = [1 1 1];
Exp.nPoints = 2;

% Detection -------------------------
Exp.DetOperator = {'z1'};
Exp.DetFreq = 0;

Exp1 = Exp;
Exp1.Dim1 = {'p1.Phase', pi};


% Options ---------------------------
Opt.SimFreq = 32;

% Function Call -----------------------------

[signal1] = spidyan(Sys,Exp1,Opt);

% Build custom IQ
Pulse.Phase = 0;
[t,IQ1] = pulse(Pulse);
Opt.dt = Opt.IntTimeStep;
[t1, IQ1] = rfmixer(t,IQ1,Exp.mwFreq,'IQshift',Opt);

Pulse.Phase = pi;
[t,IQ2] = pulse(Pulse);
Opt.dt = Opt.IntTimeStep;
[t2, IQ2] = rfmixer(t,IQ2,Exp.mwFreq,'IQshift',Opt);

PulseIQ.IQ = {IQ1 IQ2};
PulseIQ.t = {t1 t2};

%  Exp with userdefined IQ
Exp2 = Exp;
Exp2.Sequence{1} = PulseIQ;
Exp2.Dim1 = {'p1.IQ' []};

[signal2] = spidyan(Sys,Exp2,Opt);

Exp3 = Exp;
Exp3.nPoints = [2 5];
Exp3.Dim1 = {'p1.Phase', pi};
Exp3.Dim2 = {'d1' 0.05};

[signal3] = spidyan(Sys,Exp3,Opt);

Exp4 = Exp;
PulseIQ.IQ = {IQ1; IQ2};
PulseIQ.t = {t1; t2};

Exp4.Sequence{1} = PulseIQ;
Exp4.Dim1 = {'p1.IQ' []};
Exp4.Dim2 = {'d1' 0.05};
Exp4.nPoints = [2 5];

[signal4] = spidyan(Sys,Exp4,Opt);

Exp5 = Exp4;
PulseIQ.IQ = {IQ1 IQ2};
PulseIQ.t = {t1; t2};

Exp5.Sequence{1} = PulseIQ;

[signal5] = spidyan(Sys,Exp5,Opt);

ok = areequal(signal1,signal2,1e-4,'abs') && ...
  areequal(signal3{2,4},signal4{2,4},1e-4,'abs') && ...
  isequal(signal5,signal4);

warning(orig_state);
