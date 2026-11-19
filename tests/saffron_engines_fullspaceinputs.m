function ok = test() %#ok<FNDEF>

% Sequences that only the full-space engine can represent are handed to it
% with Opt.Method = 'auto' (and rejected with Opt.Method = 'subspace'),
% rather than causing an error while the engine is chosen: pulses given
% as IQ waveforms without a flip angle, and indirect dimensions over
% frequency bands of pulses

Sys.Nucs = '1H';
Sys.A = [2 3 9];  % MHz

Pulse.Type = 'rectangular';
Pulse.tp = 0.02;  % µs
Pulse.Flip = pi/2;

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0 0];  % single orientation
Opt.IntTimeStep = 1e-4;  % µs

% Pulse given as an IQ waveform (upconverted to the carrier), without Flip
OptIQ = Opt;
OptIQ.IntTimeStep = 1e-6;  % µs, resolves the lab-frame waveform
[t,IQ] = pulse(setfield(Pulse,'TimeStep',OptIQ.IntTimeStep)); %#ok<SFLD>
[t,IQ] = rfmixer(t,IQ,Exp.mwFreq,'IQshift',struct('dt',OptIQ.IntTimeStep));
IQPulse.IQ = IQ;
IQPulse.t = t;
E = Exp;
E.Sequence = {IQPulse 0.1};
[~,~,info] = saffron(Sys,E,OptIQ);
ok(1) = strcmp(info.Method,'fullspace') && contains(info.MethodReason,'IQ');
ok(2) = throwsWith(@() saffron(Sys,E,setfield(OptIQ,'Method','subspace')),'IQ');

% Indirect dimension over the frequency band of a chirp: one row per
% point, or one row of band-edge increments
Chirp.Type = 'rectangular/linear';
Chirp.tp = 0.1;  % µs
Chirp.Frequency = [-20 20];  % MHz
Chirp.Flip = pi/2;
OptC = Opt;
OptC.IntTimeStep = 2e-6;  % µs, resolves the chirp
E = Exp;
E.Sequence = {Chirp 0.1};
E.nPoints = 3;
E.Dim1 = {'p1.Frequency' [0 0; -5 5; -10 10]};  % MHz
[~,~,info] = saffron(Sys,E,OptC);
ok(3) = strcmp(info.Method,'fullspace') && contains(info.MethodReason,'frequency bands');
E.Dim1 = {'p1.Frequency' [-5 5]};  % MHz
[~,~,info] = saffron(Sys,E,OptC);
ok(4) = strcmp(info.Method,'fullspace') && contains(info.MethodReason,'frequency bands');
ok(5) = throwsWith(@() saffron(Sys,E,setfield(OptC,'Method','subspace')),'frequency bands');

% A rectangular pulse without Flip, Amplitude or IQ is still an error
E = Exp;
E.Sequence = {rmfield(Pulse,'Flip') 0.1};
ok(6) = throwsWith(@() saffron(Sys,E,Opt),'Flip');

end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch err
  ok = contains(err.message,str);
end
end
