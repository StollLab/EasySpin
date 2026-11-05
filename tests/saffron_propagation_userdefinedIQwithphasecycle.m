function ok = test() %#ok<FNDEF>

% User-defined IQ waveforms with one row per phase cycle step (array or
% cell) give the same signal as the built pulse with a phase cycle (single
% orientation, detection over the whole sequence)

orig_state = warning;
warning('off','all');

% System ------------------------
Sys.S = 1/2;


% Pulse -------------------
Pulse.Type = 'rectangular';
Pulse.TimeStep = 0.0001;
Pulse.Flip = pi;
Pulse.Frequency = 1000* 0.100;
Pulse.tp = 0.1;

PC = [0, 1; pi, -1];

% Experiment with pulses internally defined -------------------
Exp.Sequence = {Pulse 0.5 Pulse};
Opt.IntTimeStep = 0.0001; % µs
Exp.mwFreq = 33.5;
Exp.Field = mhz2mt(33600);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';
Exp.PhaseCycle = {PC []};

% Options ---------------------------
Opt.SimFreq = 32;

% Function Call -----------------------------

[~,signal1] = saffron(Sys,Exp,Opt);

% Build custom IQ
Pulse.Phase = 0;
[t,IQ1] = pulse(Pulse);
Opt.dt = Opt.IntTimeStep/10;
[~, IQ1] = rfmixer(t,IQ1,Exp.mwFreq,'IQshift',Opt);

Pulse.Phase = pi;
[t,IQ2] = pulse(Pulse);
Opt.dt = Opt.IntTimeStep/10;
[t, IQ2] = rfmixer(t,IQ2,Exp.mwFreq,'IQshift',Opt);

IQPulse.IQ = [IQ1; IQ2];
IQPulse.t = t;

%  Exp with userdefined IQ
Exp2 = Exp;
Exp2.Sequence{1}= IQPulse;

[~, signal2] = saffron(Sys,Exp2,Opt);

Exp3 = Exp;
IQPulse2.IQ = {[IQ1; IQ2]};
IQPulse2.t = {t};
Exp3.Sequence{1} = IQPulse2;

signal3 = saffron(Sys,Exp3,Opt);

ok = [areequal(signal1,signal2,1e-4,'abs') areequal(signal1,signal3,1e-4,'abs')];

warning(orig_state);
