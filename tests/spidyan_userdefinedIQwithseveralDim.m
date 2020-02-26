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
Exp.t = [0.1 0.5 0.1];
Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240; 
Opt.IntTimeStep = 0.0001; % us
Exp.mwFreq = 33.5;
Exp.DetSequence = [1 1 1];

Exp.nPoints = [2 3];

% Detection -------------------------
Exp.DetOperator = {'z1'};
Exp.DetFreq = [0]; 

Exp1 = Exp;
Exp1.Dim1 = {'p2.Phase', pi};
Exp1.Dim2 = {'p2.t', 0.1};

% Options ---------------------------
Opt.SimFreq = 32;

% Function Call -----------------------------

[signal1] = spidyan(Sys,Exp1,Opt);

% Build custom IQ
Pulse.Phase = 0;
[t,IQ11] = pulse(Pulse);
Opt.dt = Opt.IntTimeStep;
[t11, IQ11] = rfmixer(t,IQ11,Exp.mwFreq,'IQshift',Opt);

Pulse.Phase = pi;
[t,IQ21] = pulse(Pulse);
Opt.dt = Opt.IntTimeStep;
[t21, IQ21] = rfmixer(t,IQ21,Exp.mwFreq,'IQshift',Opt);

Pulse.tp = 0.2;

Pulse.Phase = 0;
[t,IQ12] = pulse(Pulse);
Opt.dt = Opt.IntTimeStep;
[t12, IQ12] = rfmixer(t,IQ12,Exp.mwFreq,'IQshift',Opt);

Pulse.Phase = pi;
[t,IQ22] = pulse(Pulse);
Opt.dt = Opt.IntTimeStep;
[t22, IQ22] = rfmixer(t,IQ22,Exp.mwFreq,'IQshift',Opt);

Pulse.tp = 0.3;

Pulse.Phase = 0;
[t,IQ13] = pulse(Pulse);
Opt.dt = Opt.IntTimeStep;
[t13, IQ13] = rfmixer(t,IQ13,Exp.mwFreq,'IQshift',Opt);

Pulse.Phase = pi;
[t,IQ23] = pulse(Pulse);
Opt.dt = Opt.IntTimeStep;
[t23, IQ23] = rfmixer(t,IQ23,Exp.mwFreq,'IQshift',Opt);

PulseIQ.IQ = {IQ11 IQ12 IQ13; IQ21 IQ22 IQ23};
PulseIQ.t = {t11 t12 t13; t21 t22 t23};

%  Exp with userdefined IQ
Exp2 = Exp;
Exp2.Sequence{3} = PulseIQ;
Exp2.Dim1 = {'p2.IQ' []};
Exp2.Dim2 = {'p2.IQ' []};

[signal2] = spidyan(Sys,Exp2,Opt);

Exp3 = Exp;
Exp3.nPoints = [4 2 3];
Exp3.Dim1 = {'d1' 0.05};
Exp3.Dim2 = {'p2.Phase', pi};
Exp3.Dim3 = {'p2.t', 0.1};

[signal3] = spidyan(Sys,Exp3,Opt);


Exp4 = Exp;

Exp4.nPoints = [4 2 3];
Exp4.Dim1 = {'d1' 0.05};
Exp4.Dim2 = {'p2.IQ' []};
Exp4.Dim3 = {'p2.IQ' []};

PulseIQ.IQ = cell(1,2,3);
PulseIQ.IQ(1,:,:) = {IQ11 IQ12 IQ13; IQ21 IQ22 IQ23};

PulseIQ.t = cell(1,2,3);
PulseIQ.t(1,:,:) = {t11 t12 t13; t21 t22 t23};

Exp4.Sequence{3} = PulseIQ;

[signal4] = spidyan(Sys,Exp4,Opt);

ok = areequal(signal1{2},signal2{2},1e-4,'abs') && areequal(signal1{5},signal2{5},1e-4,'abs') &&...
     areequal(signal3{2},signal4{2},1e-4,'abs') && areequal(signal3{5},signal4{5},1e-4,'abs');

warning(orig_state);
