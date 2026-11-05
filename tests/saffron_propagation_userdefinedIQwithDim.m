function ok = test() %#ok<FNDEF>

% User-defined IQ waveforms along indirect dimensions (Dim with 'p1.IQ')
% give the same signals as the built pulses with incremented phases, for
% one and for two indirect dimensions, and for row and column cell arrays
% of waveforms (single orientation, Sz detected in the lab frame)

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

% Experiment with pulses internally defined -------------------
Exp.Sequence = {Pulse 0.5 Pulse};
Opt.IntTimeStep = 0.0001; % µs
Exp.mwFreq = 33.5;
Exp.Field = mhz2mt(33600);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';
Exp.nPoints = 2;

% Detection -------------------------
Exp.DetOperator = {'z1'};
Exp.DetFreq = 0;

Exp1 = Exp;
Exp1.Dim1 = {'p1.Phase', pi};


% Options ---------------------------
Opt.SimFreq = 32;

% Function Call -----------------------------

signal1 = saffron(Sys,Exp1,Opt);

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

signal2 = saffron(Sys,Exp2,Opt);

% Second dimension over the delay: detection after the last pulse, so that
% the detected transients have the same length for all points
Exp3 = Exp;
Exp3.DetWindow = [0 0.05];  % µs
Exp3.nPoints = [2 5];
Exp3.Dim1 = {'p1.Phase', pi};
Exp3.Dim2 = {'d1' 0.05};

signal3 = saffron(Sys,Exp3,Opt);

Exp4 = Exp3;
PulseIQ.IQ = {IQ1; IQ2};
PulseIQ.t = {t1; t2};

Exp4.Sequence{1} = PulseIQ;
Exp4.Dim1 = {'p1.IQ' []};

signal4 = saffron(Sys,Exp4,Opt);

Exp5 = Exp4;
PulseIQ.IQ = {IQ1 IQ2};
PulseIQ.t = {t1; t2};

Exp5.Sequence{1} = PulseIQ;

signal5 = saffron(Sys,Exp5,Opt);

ok = [areequal(signal1,signal2,1e-4,'abs') areequal(signal3,signal4,1e-4,'abs') ...
  isequal(signal5,signal4)];

warning(orig_state);
