function ok = test() %#ok<FNDEF>

% A user-defined IQ waveform (row or column vector, or in a cell array)
% gives the same signal as the built pulse with a phase cycle: the phase
% cycle shifts the phase of the waveform (single orientation, detection
% over the whole sequence). Several waveforms in Pulse.IQ give an error.

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
Exp.PhaseCycle = [0 0 0; pi 0 pi];  % first pulse +x, -x

% Function Call -----------------------------

[~,signal1] = saffron(Sys,Exp,Opt);

% Build custom IQ, upconverted with the carrier
Pulse.Phase = 0;
[t,IQ] = pulse(Pulse);
Opt.dt = Opt.IntTimeStep/10;
[t,IQ] = rfmixer(t,IQ,Exp.mwFreq,'IQshift',Opt);

IQPulse.IQ = IQ;
IQPulse.t = t;

% Exp with user-defined IQ
Exp2 = Exp;
Exp2.Sequence{1} = IQPulse;
[~,signal2] = saffron(Sys,Exp2,Opt);

% column vector
Exp3 = Exp2;
Exp3.Sequence{1}.IQ = IQ(:);
[~,signal3] = saffron(Sys,Exp3,Opt);

% cell arrays
Exp4 = Exp;
IQPulse2.IQ = {IQ};
IQPulse2.t = {t};
Exp4.Sequence{1} = IQPulse2;
[~,signal4] = saffron(Sys,Exp4,Opt);

ok = [areequal(signal1,signal2,1e-4,'abs') areequal(signal1,signal3,1e-4,'abs') ...
  areequal(signal1,signal4,1e-4,'abs')];

% several waveforms in Pulse.IQ
Exp5 = Exp2;
Exp5.Sequence{1}.IQ = [IQ; -IQ];
try
  saffron(Sys,Exp5,Opt);
  ok(end+1) = false;
catch err
  ok(end+1) = contains(err.message,'single waveform');
end

warning(orig_state);
